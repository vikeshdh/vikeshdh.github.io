#' Create resistance stack (VERSION 2)
#' Liv Hemond, based off original code by Sam Muir
#' Last updated: 2024-11-29
#' 
#' This function creates a raster of resistance values to be used for circuitscape/ least cost pathway modeling, specifically for elephants in the Kavango-Zambezi Trans-frontier Conservation Area (KAZA TFCA) in Africa. The resistance raster is created by stacking multiple rasters and assigning resistance values based on the resistance rank of each cell. The resistance rank is determined by the maximum value of the stack, with additional rules for slope, & river-road crossings. These ranks are then converted to resistance values. The ranks and resistance values are based on Van de Perre 2014 (https://doi.org/10.1111/aje.12139) and additional literature review. 
#' 
#' @param landcover land cover raster with resistance rank value (1-6)
#' @param fences fence data raster with rank 6
#' @param roads road data raster
#' @param slope slope data raster in degrees
#' @param wet_season_rivers wet season rivers raster 
#' @param dry_season_rivers dry season rivers raster
#' @param river_barriers rivers with high flow that act as a barrier for movement
#' @param seasonal_rivers_buff rivers buffers that act as an attractant
#' @param ewh_02_04 ephemeral watering holes with rank 1
#' @param ewh_04_08 ephemeral watering holes with rank 2
#' @param ewh_08_1 ephemeral watering holes with rank 3
#' @param season "wet" or "dry"
#'
#' @return raster with assigned resistance values
#' 
#' @references https://doi.org/10.1111/aje.12139
#'

create_resistance_stack <- function(landcover, fences, roads, slope, wet_season_rivers, dry_season_rivers, river_barriers, river_buffers, ewh_02_04, ewh_04_08, ewh_08_1, season) {
  
  fences <- resample(fences, landcover)
  roads <- resample(roads, landcover)
  slope <- resample(slope, landcover)
  wet_season_rivers <- resample(wet_season_rivers, landcover)
  dry_season_rivers <- resample(dry_season_rivers, landcover)
  river_barriers <- resample(river_barriers, landcover)
  river_buffers <- resample(river_buffers, landcover)
  ewh_02_04 <- resample(ewh_02_04, landcover)
  ewh_04_08 <- resample(ewh_04_08, landcover)
  ewh_08_1 <- resample(ewh_08_1, landcover)
  
  if (season == "wet") { # define rules for wet season
    
    # modify river rasters where they intersect with roads
    # if a road is present, the river rank is set to the road rank; and if not, just return the river rank
    wet_season_rivers <- ifel((roads > 0) & (wet_season_rivers > 0), roads, wet_season_rivers)
    
    # stack land cover and roads
    landcover_stack <- c(landcover, roads)
    
    # find max across stack
    landcover_stack_max <- app(landcover_stack, max, na.rm = TRUE)
    
    # add rank/rank increase based on slope thresholds
    landcover_stack_max <- ifel(slope >= 15 & slope <= 30, landcover_stack_max + 1,
                                ifel(slope > 30, 6, landcover_stack_max))
    
    # add rank/rank decrease based on river & river buffers
    landcover_stack_max <- ifel(wet_season_rivers > 0, landcover_stack_max - 1, landcover_stack_max)
    landcover_stack_max <- ifel(river_buffers > 0, landcover_stack_max - 1, landcover_stack_max)
    landcover_stack_max <- ifel(river_barriers > 0, 6, landcover_stack_max)
    
    # decrease rank based on ephemeral watering holes 
    # landcover_stack_max <- ifel(ewh_02_04 == 1, landcover_stack_max - 1,
    #                             ifel(ewh_04_08 == 2, landcover_stack_max - 1,
    #                                  ifel(ewh_08_1 == 3, landcover_stack_max - 1, landcover_stack_max)))
    
    # include fences
    # this is last in case there are overlapping fence & road/river/etc cells
    landcover_stack_max <- ifel(fences > 0, 6, landcover_stack_max)
    max_resistance_rank <- app(landcover_stack_max, max, na.rm = TRUE)
    
    # round to make sure we get whole number values
    max_resistance_rank[max_resistance_rank < 1] <- 1
    max_resistance_rank <- round(max_resistance_rank)
    
  } else if (season == "dry") { # define rules for dry season
    
    # modify river rasters where they intersect with roads
    # if a road is present, the river rank is set to the road rank; and if not, just return the river rank
    dry_season_rivers <- ifel((roads > 0) & (dry_season_rivers > 0), roads, dry_season_rivers)
    
    # stack land cover and roads
    landcover_stack <- c(landcover, roads)
    
    # find max across stack
    landcover_stack_max <- app(landcover_stack, max, na.rm = TRUE)
    
    # add rank/rank increase based on slope thresholds
    landcover_stack_max <- ifel(slope >= 15 & slope <= 30, landcover_stack_max + 1,
                                ifel(slope > 30, 6, landcover_stack_max))
    
    # add rank/rank decrease based on river & river buffers
    landcover_stack_max <- ifel(dry_season_rivers > 0, landcover_stack_max - 2, landcover_stack_max)
    landcover_stack_max <- ifel(river_buffers > 0, landcover_stack_max - 2, landcover_stack_max)
    landcover_stack_max <- ifel(river_barriers > 0, 6, landcover_stack_max)
    
    
    # decrease rank based on ephemeral watering holes 
    # landcover_stack_max <- ifel(ewh_04_08 == 2, landcover_stack_max - 1,
    #                           ifel(ewh_08_1 == 3, landcover_stack_max - 2, landcover_stack_max))
    
    # include fences
    # this is last in case there are overlapping fence & road/river/etc cells
    landcover_stack_max <- ifel(fences > 0, 6, landcover_stack_max)
    max_resistance_rank <- app(landcover_stack_max, max, na.rm = TRUE)
    
    # round to make sure we get whole number values
    max_resistance_rank[max_resistance_rank < 1] <- 1
    max_resistance_rank <- round(max_resistance_rank)
    
  } else {
    stop("season must be either 'wet' or 'dry'") # do not run if season != 'wet' or 'dry'
  }
  
  # assign resistance based on rank (see reference doi)
  max_resistance <- ifel(max_resistance_rank == 1, 1,
                         ifel(max_resistance_rank == 2, 5,
                              ifel(max_resistance_rank == 3, 20,
                                   ifel(max_resistance_rank == 4, 100,
                                        ifel(max_resistance_rank == 5, 500,
                                             ifel(max_resistance_rank == 6, 2000,
                                                  ifel(max_resistance_rank == 0, .5,
                                                       ifel(max_resistance_rank == -1, 0.1,
                                                            landcover_stack_max))))))))
  
  return(max_resistance)
}
