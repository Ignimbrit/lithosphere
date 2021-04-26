#' R6 Class representing the fundamental framework for building a geological 3D-model in Lithosphere 
#' 
#' @description 
#' A Lithosphere-Object has, in xy-direction, a minimum and a maximum coordinate and a resolution as well as a predefined stratigraphic sequence
#' 
#' @details 
#' Once the model is initialized with the basic grid constraints, provide well markers indicating lithology changes and (optionally) a dem and start the interpolation process.
#' 
#' @export
#' 
#' @examples 
#' 
#' data("synthetic_dem_1") 
#' data("synthetic_welldata_1")
#' 
#' mymodel <- Lithosphere$new(
#' xmin = 0, xmax = 1000, xres = 10, 
#' ymin = 0, ymax = 1000, yres = 10, 
#' stratigraphy = c("upper_sand", "intermediate_clay", "lower_marl")
#' )
#' 
#' mymodel$set_DEM(synthetic_dem_1)
#' mymodel$set_markers(xyzg = synthetic_welldata_1[, c("x", "y", "z", "stratigraphy")])
#' mymodel$interpolate_surfaces()
#' 
#' mymodel$export_to_data_frame() 
#' 
#' 

Lithosphere <- R6::R6Class(
  "Lithosphere",
  public = list(
    
    #' @description 
    #' Create a new lithosphere model
    #' @param xmin minimum x coordinate (left border). Will be truncated to integer.
    #' @param xmax maximum x coordinate (right border). Will be truncated to integer.
    #' @param ymin minimum y coordinate (lower border). Will be truncated to integer.
    #' @param ymax maximum y coordinate (upper border). Will be truncated to integer.
    #' @param xres numeric vector of length 1 to set the horizontal resolution. Will be truncated to integer.
    #' @param yres numeric vector of length 1 to set the vertical resolution. Will be truncated to integer.
    #' @param stratigraphy character vector with names of the stratigraphic sequence (the layers) that shall be modeled sorted from the uppermost to the lowermost unit. Getting the correct order here ist really important! 
    #' @param verbose logical. Print out what you ordered?
    initialize = function(
      xmin, xmax, xres, ymin, ymax, yres, stratigraphy, verbose = TRUE 
    ){
      private$xmin <- as.integer(xmin)
      private$xmax <- as.integer(xmax)
      private$xres <- as.integer(xres)
      private$ymin <- as.integer(ymin)
      private$ymax <- as.integer(ymax)
      private$yres <- as.integer(yres)
      private$stratigraphy <- as.character(stratigraphy)
      private$has_DEM <- FALSE
      
      # build initial grid
      X = seq(from = private$xmin, to = private$xmax, by = private$xres)
      Y = seq(from = private$ymin, to = private$ymax, by = private$yres)
      
      if(verbose){
        cat(
          paste0(
            "creating grid with ", length(X), " columns and ", length(Y), " rows\n",
            "ranging from ", private$xmin, " to ", private$xmax, " by ", private$xres, " in x-direction and\n",
            "ranging from ", private$ymin, " to ", private$ymax, " by ", private$yres, " in y-direction, respectively."
          )
        )
      }
      
      private$grid <- expand.grid(X = X, Y = Y)
      
      invisible(self)
      
    },
    
    #' @description 
    #' The models printing method.
    print = function(){
      print(paste0('A "Lithosphere" geological model with '), length(private$stratigraphy), " layers:\n")
      print(paste0(private$stratigraphy))
    },
    
    #' @description 
    #' Set the models uppermost border. Usually thought to be a digital elevation model (DEM) representing the ground level.
    #' @param xyz matrix or data.frame with exactly three columns: x, y and z coordinates. 
    #' 
    set_DEM = function(xyz){
      private$DEM <- raster::extract(
        raster::rasterFromXYZ(
          xyz = xyz,
          res = c(private$xres, private$yres)
          ), 
        private$grid
        )
      private$has_DEM <- TRUE
      
      invisible(self)
    },
    
    #' @description 
    #' Add well markers
    #' @param xyzg a \code{data.frame} with four columns: the x, the y and the z coordinate of a marker point as well as a label corresponding to the stratigraphy as defined at model initialization.
    set_markers = function(xyzg){
      
      stopifnot(is.data.frame(xyzg))
      stopifnot(all(xyzg[[4]] %in% private$stratigraphy))
      
      private$marker <- data.frame(
        X = xyzg[[1]],
        Y = xyzg[[2]],
        Z = xyzg[[3]],
        G = xyzg[[4]]
      )
      
      invisible(self)
      
    },
    
    #' @description 
    #' Start the interpolation and possibly the subsequent erosion of all stratigraphic layers within a model
    #' @param engine a character vector of length 1 stating the interpolation algorithm used to transform markers into continuous surfaces. Currently available are "idw" (inverse distance weighted), "ok" (ordinary kriging), and "uk" (universal kriging - the default).
    #' @param erosion a character vector of length 1 stating the type of erosion applied to each layer (with the exception of the DEM) after interpolation. Available options are "topdown" (the default), "bottomup" and "none" (which will not cut any intersecting surfaces at all).
    interpolate_surfaces = function(engine = "uk", erosion = "topdown"){
      horizon_marker <- purrr::map(
        private$stratigraphy,
        function(x){
          private$marker[private$marker[["G"]] == x, ]
        }
      )
      
      interpolator <- private$get_interpolator(which = engine)
      
      horizons_raw <- purrr::map_dfc(
        horizon_marker, interpolator
      )
      
      private$surface_df <- cbind(private$grid, horizons_raw)
      
      private$erode(direction = erosion)
      
      invisible(self)
    },
    
    #' @description 
    #' Return a data.frame of interpolated surface coordinates
    #' @param shape a character vector of length 1: either "long" (the default) or "wide".
    #' 
    #' @return 
    #' A \code{data.frame} formated either in line with \code{\link[tidyr]{pivot_wider}} or \code{\link[tidyr]{pivot_longer}}, depending on which format was specified via the \code{shape} argument.
    export_to_data_frame = function(shape = "long"){
      
      if(private$has_DEM){
        
        raw <- cbind(
          cbind(
            private$grid,
            DEM = private$DEM
          ),
          private$surface_df[, private$stratigraphy]
        )
        
      } else {
        raw <- private$surface_df
      }
      
      if(shape == "wide"){
        raw
      } else if(shape == "long"){
        longdf <- tidyr::pivot_longer(
          raw, !c(X, Y), names_to = "strat", values_to = "Z"
        )
        
        longdf[, c("strat", "X", "Y", "Z")]
      } else{
        stop('shape formats "long" and "wide" supported only')
      }
    },
    
    #' @description 
    #' Plot a 3D-representation of the model
    #' @param surface_colors A character vector of color names as accepted by R. Ideally, the length of the vector should correspond to the number of stratigraphic units (+DEM) of the model but selecting any other number of colors will be dealt with, too.
    #' @param exag_fct numeric specifying the exaggeration factor of the visualization in z-direction. Defaults to 1 but maybe try 3.
    plot3d = function(surface_colors = NA, exag_fct = 1){
      raw <- self$export_to_data_frame(shape = "wide")
      
      ncolors = length(private$stratigraphy) + as.integer(private$has_DEM)
      
      if(all(is.na(surface_colors))){
        surface_colors <- base::sample(
          grDevices::colors(distinct = TRUE),
          ncolors
        )
      } else if(length(surface_colors) < ncolors){
        base::message("Too few colors provided. Padding with random colors")
        
        ncolorsmissing <- ncolors - length(surface_colors)
        surface_colors <- c(
          surface_colors, 
          base::sample(
            grDevices::colors(distinct = TRUE),
            ncolorsmissing
          )
        )
      } else if(length(surface_colors) > ncolors){
        base::message("Too many colors provided. Truncating...")
        surface_colors <- surface_colors[1:ncolors]
      }
      
      for(i in seq_along(surface_colors)){
        z <- raster::as.matrix(raster::rasterFromXYZ(raw[, c(1, 2, i+2)]))
        
        rgl::surface3d(
          x = private$xres*(1:nrow(z)),
          y = private$yres*(1:ncol(z)),
          z = z*exag_fct,
          color = surface_colors[i]
        )
      }
    },
    
    #' @field idp inverse distance power 
    idp = 0.5
  ),
  private = list(
    xmin = NA_integer_,
    xmax = NA_integer_,
    xres = NA_integer_,
    ymin = NA_integer_,
    ymax = NA_integer_,
    yres = NA_integer_,
    stratigraphy = NA_character_,
    grid = NA,
    DEM = NA,
    has_DEM = FALSE,
    marker = NA,
    surface_df = NA,
    get_interpolator = function(which){
      if(which == "idw"){
        function(hmarkerframe){
          setNames(
            data.frame(
              res = predict(
                gstat::gstat(formula = Z~1, locations = ~X+Y, data = hmarkerframe, set = list(idp = self$idp)),
                newdata = private$grid
              )[["var1.pred"]]
            ),
            hmarkerframe[["G"]][1]
          )
        }
      } else if(which == "ok"){
        function(hmarkerframe){
          setNames(
            data.frame(
              res = predict(
                gstat::gstat(
                  formula = Z~1, locations = ~X+Y, data = hmarkerframe,
                  model = automap::autofitVariogram(
                    formula = Z~1,
                    input_data = sf::as_Spatial(
                      sf::st_as_sf(hmarkerframe, coords = c("X", "Y"))
                    )
                  )$var_model
                ),
                newdata = private$grid
              )[["var1.pred"]]
            ),
            hmarkerframe[["G"]][1]
          )
        }
      } else if(which == "uk"){
        function(hmarkerframe){
          setNames(
            data.frame(
              res = predict(
                gstat::gstat(
                  formula = Z~X+Y, locations = ~X+Y, data = hmarkerframe,
                  model = automap::autofitVariogram(
                    formula = Z~coords.x1+coords.x2,
                    input_data = sf::as_Spatial(
                      sf::st_as_sf(hmarkerframe, coords = c("X", "Y"))
                    )
                  )$var_model
                ),
                newdata = private$grid
              )[["var1.pred"]]
            ),
            hmarkerframe[["G"]][1]
          )
        }
      }
    },
    erode = function(direction = "topdown"){
      
      if(private$has_DEM){
        surf_pre <- cbind(
          DEM = private$DEM,
          private$surface_df[, private$stratigraphy]
        )
      } else { 
        surf_pre <- private$surface_df[, private$stratigraphy]
      }
      
      if(direction == "topdown"){
        surf_post <- private$erode_topdown(surf_pre)
        private$surface_df <- cbind(private$surface_df[, c("X", "Y")], surf_post)
      } else if(direction == "bottomup"){
        surf_post <- private$erode_bottomup(surf_pre)
        private$surface_df <- cbind(private$surface_df[, c("X", "Y")], surf_post)
      } else if(direction == "none"){
        private$surface_df <- private$surface_df
      } else {
        stop('erosion must be either, "topdown", "bottomup" or "none"')
      }
      
    },
    erode_topdown = function(pre_erosion_surface_frame){
      
      wdf <- pre_erosion_surface_frame
      
      for(i in 2:ncol(pre_erosion_surface_frame)){
        upper <- wdf[[i-1]]
        lower <- wdf[[i]]
        
        wdf[[i]] <- dplyr::case_when(
          upper < lower ~ upper,
          upper >= lower ~ lower
        )
      }
      
      if(private$has_DEM){
        wdf$DEM <- NULL
      }
      
      wdf
      
    },
    erode_bottomup = function(pre_erosion_surface_frame){
      
      wdf <- pre_erosion_surface_frame
      
      # DEM always wins the erosion game 
      if(private$has_DEM){
        wdf[[2]] <- dplyr::case_when(
          wdf[[1]] < wdf[[2]] ~ wdf[[1]],
          wdf[[1]] >= wdf[[2]] ~ wdf[[2]]
        )
      }
      
      for(i in rev(2:(ncol(pre_erosion_surface_frame)-1))){
        upper <- wdf[[i]]
        lower <- wdf[[i+1]]
        
        wdf[[i]] <- dplyr::case_when(
          upper < lower ~ lower,
          upper >= lower ~ upper,
        )
      }
      
      if(private$has_DEM){
        wdf$DEM <- NULL
      }
      
      wdf
    }
  )
)