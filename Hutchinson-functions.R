# Function 'hutchinson', 'hutchinson_cat'
### These functions work in combination with the 'e_space' and 'e_space_cat' functions available in the E_space_functions script
### It allows the user to create transects in any environmental space and transfer them to the geographic space 
### or viceversa. It can be used to further explore and select areas on the environmental dimensions used to calibrate  
### the model or other environmental dimensions representing knowledge of the species niche. 

#
#' @param EtoG logical. If TRUE, E-space is plotted first to select transects and then,
#'   transects are drawn in G-space. If FALSE, G-space is plotted first.
#' @param data dataframe with the following columns: longitude, latitude, 
#' and environmental variables (one per column) for 'hutchinson' as it is obtained with 
#' 'e_space' function, or a dataframe with longitude, latitude, suitability category (from a thresholded raster)
#' and environmental variables as columns (one per column) for 'hutchinson_cat',
#' as it is obtained from applying the function 'e_space_cat'
#' @param calls vector of length two that indicates the columns that contain the
#' environmental variables to be used in the plots.
#' @param plyg polygon that delimits the geographical space of interest.
#' @param ntr number of transects to be drawn.
#' @param col.use vector of lenght two with the colors to be used in plots
#' 
#' @return
#' \code {hutchinso} returns a database with coordinates, values, 
#' and the selected tracks
#' \code{hutchinson_cat} returns a database with coordinates, values, categories,
#' and the selected tracks.
#' 
#' @describeIn hutchinson plots the E-space/G-space and allows the user
#' to select transects that are then represented in the G-space/E-space.
#' 
# CODE hutchinson -----------
# Dependencies: maptools, raster, rgdal, sp, gatepoints

hutchinson = function (EtoG, data, calls, plyg, ntr){
  # create object to save information about the transects created
  transects <- vector("list",length=ntr)
  if(EtoG==T){ # from E to G
    # Plot 1: E-space
    dev.new()
    plot(data[,calls], col = 'grey', main="E-space", pch = 3,
         cex = 0.5, xlab = colnames(data[calls[1]]), ylab = colnames(data[calls[2]]))
        # user defines transects
    for(j in 1:ntr){
      if(j == 1){
        print("Instructions:")
        print("(1) Select vertices of transect by clicking on the map")
        print("(2) Right-click and select 'stop' to finish drawinng a transect")
        print("(3) Repeat (1) and (2) for eah transect")
      }
      # left click to create points, right click to finish
      ss = fhs(data.frame(data[,calls]),pch=1+j)
      # subsetting the entire dataset according to the transect in which they were selected 
      transects[[j]] = cbind(as.matrix(data[ss,]),rep(j,length(ss)))
    }
    all.tr = data.frame(do.call(rbind,transects)) #I left it as a dataframe, originally matrix  #POTENTIALLY ELIMINATE DUPLICATES 
    colnames(all.tr) = c(colnames(data), 'transects') 
    # Plot 2: G-space
    dev.new()
    # add region limits
    plot(plyg, main = "G-space")
    points (all.tr[,1:2], col = c(1+all.tr[,ncol(all.tr)]), pch = 1+all.tr[,ncol(all.tr)], cex = 0.8)
    tr.names <- paste("Transect",1:ntr)
    legend('bottomleft', legend=tr.names, pch = 1+(1:ntr), col = c(1+(1:ntr)))
    # save new matrix with selected points
    return(all.tr)
  }
  if(EtoG==F){ # from G to E
    # Convert geographic coordinates into a Spatial object
    data("wrld_simpl", package = "maptools")
    sp_ob = sp::SpatialPointsDataFrame(data[,calls], data, proj4string = crs(wrld_simpl))
    #dim(sp_ob)
    # Plot 1: G-space
    dev.new()
    # plot points that cover the environmental region and are identified with
    # a suitability category (different colors)
    plot(sp_ob@data[,1:2], main = "G-space", pch = 3,
         col = 'grey', cex = 0.8)
    # add region limits
    plot (plyg, add = T)
    for(j in 1:ntr){
      if(j == 1){
        print("Instructions:")
        print("(1) Select vertices of transect by clicking on the map")
        print("(2) Right-click and select 'stop' to finish drawinng a transect")
        print("(3) Repeat (1) and (2) for eah transect")
      }
      # left click to create points, right click to finish
      ss = fhs(data.frame(sp_ob@data[,1], sp_ob@data[,2]),pch=1+j)
      # subsetting the entire dataset according to the transect in which they were selected 
      transects[[j]] = cbind(as.matrix(data[ss,]),rep(j,length(ss)))
    }
    all.tr = data.frame (do.call(rbind,transects))
    colnames(all.tr) <- c(colnames(data),"transects")
    # Plot 2: E-space
    dev.new()
    plot(all.tr[,calls[1]], all.tr[,calls[2]], cex = 0.8, main="E-space",
         col = c(1+all.tr[,dim(all.tr)[2]]), pch = 1+all.tr[,dim(all.tr)[2]],
         xlab = paste ('Env_var1:', colnames(all.tr[calls[1]])), 
         ylab = paste ('Env_var2:', colnames(all.tr[calls[2]])))
    tr.names = paste ('Transect', 1:ntr)
    legend ('bottomleft', legend = tr.names, pch = 1+(1:ntr), col =  c(1+(1:ntr)))
    return(all.tr)
  }
}
#
#' @describeIn hutchinson_cat plots the E-space/G-space from a categorized raster
#' and allows the user to select transects that are then represented in the G-space/E-space.
#' 
# CODE hutchinson_cat ---------
# Dependencies: maptools, sp, raster, rgdal, gatepoints
hutchinson_cat <- function(EtoG, data, calls, plyg, ntr, col.use = NULL){
  if(is.null(col.use)){
    print("Please define 'col.use' using two colors")
  } else{
    # defining color ramp for plotting: 
    pal5 = colorRampPalette(col.use)
    # create object to save information about the transects created
    transects <- vector("list",length=ntr)
    if(EtoG==T){ # from E to G
      # determine the number of categories
      catnum = length(unique(data[,3]))
      # Plot 1: E-space
      dev.new()
      plot(data[,calls], col = pal5(catnum)[data[,3]], main="E-space", pch = 15,
           cex = 0.7, xlab = colnames(data[calls[1]]), ylab = colnames(data[calls[2]]))
      suit_class = paste("Suitability",unique(data[,3]))
      legend('topleft', legend = suit_class, pch = 15, cex = 1, col = pal5(catnum))
      # user defines transects
      for(j in 1:ntr){
        if(j == 1){
          print("Instructions:")
          print("(1) Select vertices of transect by clicking on the map")
          print("(2) Right-click and select 'stop' to finish drawinng a transect")
          print("(3) Repeat (1) and (2) for eah transect")
        }
        # left click to create points, right click to finish
        ss = fhs(data.frame(data[,calls]),pch=1+j)
        # subsetting the entire dataset according to the transect in which they were selected 
        transects[[j]] = cbind(as.matrix(data[ss,]),rep(j,length(ss)))
      }
      all.tr = data.frame(do.call(rbind,transects))  
      colnames(all.tr) <- c(colnames(data),"transects")
      # Plot 2: G-space
      dev.new()
      # add region limits
      plot(plyg, main = "G-space")
      points (all.tr[,1:2], col = pal5(catnum)[all.tr[,3]], pch = 1+all.tr[,ncol(all.tr)], cex = 0.8)
      tr.names <- paste("Transect",1:ntr)
      legend('bottomleft', legend=tr.names, pch = 1+(1:ntr))
      # save new matrix with selected points
      return(all.tr)
    }
    if(EtoG==F){ # from G to E
      # determine the number of categories
      catnum = length(unique(data[,3]))
      # Convert geographic coordinates into a Spatial object
      data("wrld_simpl", package = "maptools")
      sp_ob = sp::SpatialPointsDataFrame(data[,calls], data, proj4string = crs(wrld_simpl))
      dim(sp_ob)
      # Plot 1: G-space
      dev.new()
      # plot points that cover the environmental region and are identified with
      # a suitability category (different colors)
      plot(sp_ob@data[,1:2], main = "G-space", pch = 15,
            col = pal5(length(unique(sp_ob@data[,3])))[sp_ob@data[,3]], cex = 0.8)
      # add region limits
      plot (plyg, add = T)
      suit_class = paste("Suitability",unique(sp_ob@data[,3]))
      legend('bottomleft', legend=suit_class, pch = 15,
             col= pal5(length(unique(sp_ob@data[,3]))), cex=0.9)
       for(j in 1:ntr){
        if(j == 1){
          print("Instructions:")
          print("(1) Select vertices of transect by clicking on the map")
          print("(2) Right-click and select 'stop' to finish drawinng a transect")
          print("(3) Repeat (1) and (2) for eah transect")
        }
        # left click to create points, right click to finish
        ss = fhs(data.frame(sp_ob@data[,1], sp_ob@data[,2]),pch=1+j)
        # subsetting the entire dataset according to the transect in which they were selected 
        transects[[j]] = cbind(as.matrix(data[ss,]),rep(j,length(ss)))
      }
      all.tr = data.frame (do.call(rbind,transects))
      colnames(all.tr) <- c(colnames(data),"transects")
      # Plot 2: E-space
      dev.new()
      plot(all.tr[,calls[1]], all.tr[,calls[2]], cex = 0.8, main="E-space",
           col = pal5(catnum)[all.tr[,3]], pch = 1+all.tr[,dim(all.tr)[2]],
           xlab = paste ('Env_var1:', colnames(all.tr[calls[1]])), 
           ylab = paste ('Env_var2:', colnames(all.tr[calls[2]])))
      tr.names = paste ('Transect', 1:ntr)
      legend ('bottomleft', legend = tr.names, pch = 1+(1:ntr))
      return(all.tr)
    }
  }
}


# END
# Daniel Romero-Alvarez & Laura Jimenez, 2020