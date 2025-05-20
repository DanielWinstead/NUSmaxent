# This script is designed to be able to run as a batch job or on RStudio for Parallel Processing
# Some notes and code used here are from the ENMevaluation vignette (https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0-vignette.html)
# Inputs: csv of species names of interest, Training rasters (Available in GitHub Repository), GBIF login info, ENMeval settings
# Ouptuts: ENMeval object(rds); predictive model outputs (rds) for 2011-40 and 2041-70;
#          Maps of optimal maxent suitability model, 2011-40 and 2041-70 predictions,
#          and map of difference in suitabilities(geoTIFF); csv of summary statistics
#          which include average change in global suitability, change in suitability regionally
#          root mean squares difference between suitability for predictions.
# Environment: conda environment (provided in GitHub Repository).
# "###" indicates inputs needed
# "#" indicates notes

library(ENMeval)
library(dplyr)
library(raster)
library(CoordinateCleaner)
library(rgbif)
library(rasterVis)
library(sf)
library(doParallel)
library(readr)
library(foreach)
library(fs)

### Enter your working directory
setwd("")

### enter your GBIF login info
Sys.setenv(GBIF_USER = "",
           GBIF_PWD  = "",
           GBIF_EMAIL = "")

envs.files <- c(list.files(path="EnvironmentFiles/allignedRasters_Control", full.names=TRUE))

# Find the descriptions of the bioclimatic variables here: 
# https://www.worldclim.org/data/bioclim.html

envs.rasters <- lapply(envs.files, terra::rast)
envs <- terra::rast(envs.rasters)

# create function to be used in loop

process_species <- function(species_name) {
  print(paste("Processing species:", species_name))
  
  # Replace spaces with underscores for directory and file naming
  species_name_underscore <- gsub(" ", "_", species_name)
  output_dir <- paste0("/storage/group/mgj2/default/Species_Outputs/", species_name_underscore)
  
  # Check if the output directory exists, if not, create it
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    print(paste("Directory created at:", output_dir))
  } else {
    print(paste("Directory already exists at:", output_dir))
  }
  
  # Define the path to the .rds file
  rds_file_path <- paste0(output_dir, "/", species_name_underscore, "_e_mx_soil.rds")
  
  # Check if the .rds file already exists
  if (file.exists(rds_file_path)) {
    print(paste("Skipping species:", species_name, "- .rds file already exists"))
    return() # Skip processing this species
  }
  
  ########################### GBIF DOWNLOAD
  download_key <- occ_download(pred("taxonKey", name_backbone(species_name)$usageKey), # Retrieve taxonKey dynamically
                               pred("hasCoordinate", TRUE),  # Only records with coordinates
                               pred_in("basisOfRecord", c("HUMAN_OBSERVATION","LIVING_SPECIMEN","OCCURRENCE")), 
                               format = "SIMPLE_CSV" # Output as a CSV
  )
  # Check status in a loop (polling)
  status <- "PREPARING"
  while (status %in% c("PREPARING", "RUNNING")) {
    Sys.sleep(60)  # Wait 60 seconds before checking again (adjust as needed)
    status <- occ_download_meta(download_key)$status
    cat("Download status:", status, "\n")
  }
  
  # Once complete, download the file
  if (status == "SUCCEEDED") {
    occ_download_get(download_key, path = output_dir)
    cat("Download completed successfully!\n")
  } else {
    cat("Download failed with status:", status, "\n")
  }
  meta <- occ_download_meta(download_key)
  meta_text <- capture.output(meta)
  writeLines(meta_text, file.path(output_dir, paste0(download_key, "_citation.txt")))

  zip_file <- file.path(output_dir, paste0(download_key, ".zip"))
  unzip(zip_file, exdir = output_dir)
  csv_file <- file.path(output_dir, paste0(download_key, ".csv"))
  data <- read.csv(csv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Check that data is populated with occurrences
  if (!is.null(data) && nrow(data) > 0) {
    occsdf <- as_tibble(data)
    
    if ("issue" %in% colnames(occsdf)) {
      occsdf <- subset(occsdf, issue != "zerocd" & issue != "cum" & issue != "bri")
      occsdf <- occsdf[!is.na(occsdf$decimalLatitude) & !is.na(occsdf$decimalLongitude), ]
    } else {
      print(paste("Warning: 'issue' column missing for species:", species_name))
    }
    
    # Continue processing only if occsdf is not empty
    if (nrow(occsdf) > 0) {
      occs.clean <- clean_coordinates(occsdf, lat = "decimalLatitude", lon = "decimalLongitude", 
                                      outliers_method = "distance")
      
      occs <- filter(occs.clean, .summary == "TRUE" & decimalLatitude > -60)[,22:23]
      
      # Remove occurrences that have the same coordinates
      occs <- occs[!duplicated(occs),] %>% as.data.frame()
      occs <- occs %>% dplyr::select(2,1)
      
      occs.sp <- terra::vect(occs, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(envs[[24]]))
      occs.cells <- terra::cellFromXY(envs[[24]], terra::crds(occs.sp))
      occs.cellDups <- duplicated(occs.cells)
      occs <- occs[!occs.cellDups,]
      if (count(occs)<10){
        print("Species has less than 10 occurrence points.")
        return()
      }
      
     
      occs.sf <- sf::st_as_sf(occs, coords = c("decimalLongitude","decimalLatitude"), crs = raster::crs(envs))
      
      # Project point data to an equal-area projection. This is ideal for buffering.
      # We use the typical Eckert IV projection.
      eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
      occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)
      
      occs.buf <- sf::st_buffer(occs.sf, dist = 500000) %>% 
        sf::st_union() %>% 
        sf::st_sf() %>%
        sf::st_transform(crs = raster::crs(envs))
      
      # Crop environmental rasters to match the study extent
      envs.bg <- raster::crop(envs, occs.buf)
      # Next, mask the rasters to the shape of the buffers
      envs.bg <- raster::mask(envs.bg, occs.buf)
      
      # Randomly sample 10,000 background points from one background extent raster 
      # (only one per cell without replacement).
      raster_layer <- raster(envs.bg[[9]])
      bg <- dismo::randomPoints(raster_layer, n = 10000) %>% as.data.frame()
      colnames(occs) <- c("longitude", "latitude")
      colnames(bg) <- colnames(occs)
      
      ############ RUNNING ENMeval
      
      # ENMevaluation can only accept rasterstack
      envs.stack <- raster::stack(envs)
      
      ### Setup cluster for parallel processing
      cl <- makeCluster(9, outfile = paste0(species_name,"_parallel_log.txt"))
      registerDoParallel(cl)
      # make safe function
      tune_args <- list(
        fc = c("L", "LQ", "LQH", "H"), 
        rm = 1:5
      )
      # Safe wrapper function with automatic retry if Hinge fails due to too few occurrences
      safe_ENMevaluate <- function(occs, envs, bg, tune_args, ...) {
        tryCatch(
          {
            # Try running with Hinge first
            result <- ENMevaluate(occs = occs, envs = envs, bg = bg, tune.args = tune_args, ...)
            return(result)  # If successful, return result
          },
          error = function(e) {
            message("Error encountered: ", conditionMessage(e))
            
            # If error is due to Hinge, retry without Hinge
            if (grepl("non-conformable arguments", conditionMessage(e))) {
              message("Retrying without Hinge (fc = H)...")
              
              new_tune_args <- list(
                fc = c("L", "LQ", "LQH"), # Remove "H"
                rm = 1:5
              )
              
              return(tryCatch(
                ENMevaluate(occs = occs, envs = envs, bg = bg, tune.args = new_tune_args, ...),
                error = function(e2) {
                  message("Retry without Hinge also failed: ", conditionMessage(e2))
                  return(NA)
                }
              ))
            }
            
            return(NA)  # If error is unrelated, return NA
          }
        )
      }
      
      ### Define ENMevaluate parameters
      eval_params <- list(
        algorithm = 'maxnet', 
        partitions = 'randomkfold',
        partition.settings = list(kfolds = 10), 
        parallel = TRUE, 
        updateProgress = FALSE, 
        numCores = 9, 
        parallelType = "doParallel", 
        quiet = F
      )
      
      
      # Run the evaluation
      
      e.mx_soil <- do.call(safe_ENMevaluate, c(list(occs = occs, envs = envs.stack, bg = bg, tune_args = tune_args), eval_params))
      
      stopCluster(cl)
      write_rds(e.mx_soil, paste0(output_dir,"/",species_name_underscore, "_e_mx_soil.rds"))
  
      
      # Select the optimal model based on delta.AICc
      res <- eval.results(e.mx_soil)
      opt.aicc <- res %>% filter(delta.AICc == 0)
      if (nrow(opt.aicc) == 0) {
        # If no model has delta AICc of 0, select the one with the lowest AICc
        opt.aicc <- res %>% filter(AICc == min(AICc))
      }
      # Tie breaker is highest test AUC
      if (nrow(opt.aicc) > 1) {
        opt.aicc <- opt.aicc %>% filter(auc.test == max(auc.test))
      }
      print(paste("Optimal model for species", species_name, "based on AICc:", opt.aicc$tune.args))
      
      
      # use optimal model to create geotiff and use for predictions
      envs.bg <- raster(envs.bg)
      optimal_model <- e.mx_soil@models[[opt.aicc$tune.args]]
      optimal_Ras <- eval.predictions(e.mx_soil)[[opt.aicc$tune.args]]
      optimal_Ras <- crop(optimal_Ras, extent(envs.bg))
      optimal_Ras <- mask(optimal_Ras, envs.bg)
      metadata(optimal_Ras)$units <- "Suitability Index (cloglog)"
      metadata(optimal_Ras)$title <- paste0(species_name," ",opt.aicc$tune.args)
      writeRaster(optimal_Ras, paste0(output_dir, "/Suitability_current_",species_name_underscore, ".tif"), 
                  
                  format = "GTIFF",
                  overwrite = TRUE)
      # Get predictions for the current time period (2011-2040)
      PredEnvs_11_40 <- raster::stack(c(list.files(path = "EnvironmentFiles/allignedRasters_11-40", full.names = TRUE)))
      names(PredEnvs_11_40) <- names(envs)
      
      pred_11_40 <- dismo::predict(PredEnvs_11_40, optimal_model, clamp = TRUE, type = "cloglog")
      write_rds(pred_11_40, paste0(output_dir, "/Pred_", species_name, "_11-40.rds"))
      metadata(pred_11_40)$units <- "Suitability Index (cloglog)"
      metadata(pred_11_40)$title <- paste0(species_name," Prediction for 2011-2040 using GFDL ESM4 ssp370")
      pred_11_40 <- crop(pred_11_40, extent(envs.bg))
      pred_11_40 <- mask(pred_11_40, envs.bg)
      writeRaster(pred_11_40, paste0(output_dir, "/Suitability_11_40_", species_name_underscore, ".tif"), 
                  
                  format = "GTiff",
                  overwrite = TRUE)
      
      
      # Get predictions for the future time period (2041-2070)
      PredEnvs_41_70 <- raster::stack(c(list.files(path = "EnvironmentFiles/allignedRasters_41-70", full.names = TRUE)))
      names(PredEnvs_41_70) <- names(envs)
      
      pred_41_70 <- dismo::predict(PredEnvs_41_70, optimal_model, clamp = TRUE, type = "cloglog")
      
      write_rds(pred_41_70, paste0(output_dir, "/Pred_", species_name, "_41-70.rds"))
      metadata(pred_41_70)$units <- "Suitability Index (cloglog)"
      metadata(pred_41_70)$title <- paste0(species_name," Prediction for 2041-2070 using GFDL ESM4 ssp370")
      pred_41_70 <- crop(pred_41_70, extent(envs.bg))
      pred_41_70 <- mask(pred_41_70, envs.bg)
      writeRaster(pred_41_70, paste0(output_dir, "/Suitability_41_70_", species_name_underscore, ".tif"), 
                  
                  format = "GTiff",
                  overwrite = TRUE)
      
      # Calculate global max from all three rasters
      max_current <- cellStats(optimal_Ras, 'max', na.rm = TRUE)
      max_11 <- cellStats(pred_11_40, 'max', na.rm = TRUE)
      max_41 <- cellStats(pred_41_70, 'max', na.rm = TRUE)
      global_max <- max(c(max_current, max_11, max_41), na.rm = TRUE)
      
      # Function to normalize raster and save
      normalize_and_save <- function(r, out_path, global_max) {
        r_norm <- r / global_max
        r_norm[r_norm > 1] <- 1
        writeRaster(r_norm, out_path, format = "GTiff", overwrite = TRUE)
      }
      
      # Build normalized filenames
      norm_path_current <- paste0(output_dir, "/Suitability_current_", species_name_underscore, "_normalized.tif")
      norm_path_11 <- paste0(output_dir, "/Suitability_11_40_", species_name_underscore, "_normalized.tif")
      norm_path_41 <- paste0(output_dir, "/Suitability_41_70_", species_name_underscore, "_normalized.tif")
      
      # Normalize and save
      normalize_and_save(optimal_Ras, norm_path_current, global_max)
      normalize_and_save(pred_11_40, norm_path_11, global_max)
      normalize_and_save(pred_41_70, norm_path_41, global_max)
      
      
      # Calculate the difference in suitability between the two time periods and save it
      # need to add units and titles for ID
      suitability_diff_11 <- pred_11_40 - optimal_Ras
      metadata(suitability_diff_11)$units <- "Change in Suitability Index (cloglog)"
      metadata(suitability_diff_11)$title <- paste0(species_name," Difference from currenct conditions to 2011-2040 using GFDL ESM4 ssp370")
      writeRaster(suitability_diff_11, paste0(output_dir, "/Suitability_diff_11_", species_name_underscore, ".tif"), 
                  
                  format = "GTiff",
                  overwrite = TRUE)
      suitability_diff_41 <- pred_41_70 - optimal_Ras
      metadata(suitability_diff_41)$units <- "Change in Suitability Index (cloglog)"
      metadata(suitability_diff_41)$title <- paste0(species_name," Difference from currenct conditions to 2041-2070 using GFDL ESM4 ssp370")
      writeRaster(suitability_diff_11, paste0(output_dir, "/Suitability_diff_41_", species_name_underscore, ".tif"), 
                  
                  format = "GTiff",
                  overwrite = TRUE)
    
      # Extract numeric values from the raster layers
      suitability_diff_11_vals <- values(suitability_diff_11)
      suitability_diff_41_vals <- values(suitability_diff_41)
      optimal_Ras_vals <- values(optimal_Ras)
      pred_11_40_vals <- values(pred_11_40)
      pred_41_70_vals <- values(pred_41_70)
      optimal_Ras_rspl_vals <- values(optimal_Ras_rspl)
      # Components of Stability Index SI
      # Get single mean values
      # Compute mean suitability per scenario
      mean_current <- mean(optimal_Ras_vals, na.rm = TRUE)
      mean_future_1 <- mean(pred_11_40_vals, na.rm = TRUE)
      mean_future_2 <- mean(pred_41_70_vals, na.rm = TRUE)
      
      # Calculate range between scenario means
      scenario_means <- c(mean_current, mean_future_1, mean_future_2)
      range_between_scenarios <- max(scenario_means, na.rm = TRUE) - min(scenario_means, na.rm = TRUE)
      
      #Calculate and save summary statistics
      summary_stats <- data.frame(
        species = species_name,
        avg_change_11 = mean(suitability_diff_11_vals, na.rm = TRUE),
        avg_change_41 = mean(suitability_diff_41_vals, na.rm = TRUE),
        rms_c_11 = sqrt(mean((optimal_Ras_vals - pred_11_40_vals)^2, na.rm = TRUE)),
        rms_11_41 = sqrt(mean((pred_41_70_vals - pred_11_40_vals)^2, na.rm = TRUE)),
        rms_c_41 = sqrt(mean((pred_41_70_vals - optimal_Ras_vals)^2, na.rm = TRUE)),
        SI_C_Future = mean(scenario_means, na.rm = TRUE) / (range_between_scenarios + 1),
      )
      write.csv(summary_stats, paste0(output_dir, "/SummaryStats_", species_name_underscore, ".csv"))
      
      ###################################################################################
      
      
      #GENERATING METADATA
      # Generate a rangeModelMetadata object based on the information stored in the 
      # ENMevaluate object.
      # Generate a rangeModelMetadata object based on the ENMevaluate object
      rmm <- eval.rmm(e.mx_soil)
      
      ### Author info
      
      rmm$authorship$rmmName <- ""
      rmm$authorship$contact <- ""
      
      # Populate model selection rules
      rmm$model$selectionRules <- "lowest AICc, then highest AUC.test if equal AICc"
      
      # Add details about the optimal model
      rmm$model$finalModelSettings <- paste(opt.aicc$tune.args, collapse = ", ")
      rmm$model$evaluationMetrics <- list(
        AICc = opt.aicc$AICc,
        deltaAICc = opt.aicc$delta.AICc,
        validationAUC = opt.aicc$Mean.AUC
      )
      
      # Prediction metadata
      rmm$prediction$current <- list(
        extent = paste(extent(optimal_Ras), collapse = ", "),
        resolution = paste(res(optimal_Ras), collapse = "x"),
        minVal = cellStats(optimal_Ras, "min"),
        maxVal = cellStats(optimal_Ras, "max"),
        units = "Suitability Index (cloglog)"
      )
      
      rmm$prediction$future <- list(
        timePeriods = c("2011-2040", "2041-2070"),
        extent = paste(extent(pred_11_40), collapse = ", "),
        resolution = paste(res(pred_11_40), collapse = "x"),
        units = "Suitability Index (cloglog)"
      )
      
      rmm$prediction$nuclear <- list(
        timePeriods = c("years 01-05","years 06-10","years 11-15","years 16-20"),
        extent = paste(extent(pred_150_1), collapse = ", "),
        resolution = paste(res(pred_150_1), collapse = "x"),
        units = "Suitability Index (cloglog)"
      )
      
      # Input data metadata
      rmm$data$occurrences <- list(
        source = "GBIF via spocc",
        totalPoints = nrow(occs),
        filteringSteps = "Removed duplicates, cell duplicates, and spatial outliers. Use CleanCoordinates() from CoordinateCleaner with outlier method = distance. All points below 60°S were removed."
      )
      rmm$data$environment <- list(
        source = "WorldClim/CHELSA/GFDL ESM4 SSP370/SoilGrids/ORNL DACC/ Coupe et al. 2019",
        variables = names(envs),
        resolution = paste(res(envs), collapse = "x"),
        units = c("Temperature (°C)", "Precipitation (mm)", "Slope (degrees)", "Elevation (Meters)", "Soil composition (%)", "Soil acidity (pH)", "Depth (Meters)")
      )
      
      # Processing metadata
      rmm$process$software <- list(
        R = R.version.string,
        ENMeval = packageVersion("ENMeval"),
        MaxEnt = "MaxEnt (maxnet R package implementation)"
      )
      rmm$process$computing <- list(
        environment = "HPC cluster",
        numCores = 9,
        parallelLibrary = "doParallel"
      )
      
      # Save the metadata to a CSV file
      rmm_file <- paste0(output_dir, "/Metadata_", species_name, ".csv")
      rangeModelMetadata::rmmToCSV(rmm, rmm_file)
      print(paste("Metadata saved to:", rmm_file))
      
      gc()
    } else {
      cat("No valid occurrences for species:", species_name, "\n")
      return(NULL)
    }
  } else {
    cat("No data found for species:", species_name, "\n")
    return(NULL)
  }
}

### create loop and read in species name list

NUSdf <- read.csv("")

NUS <- paste0(NUSdf$Genus," ", NUSdf$Species)

#run function in loop through species names

for (species in NUS) {
  process_species(species)  
  
}
