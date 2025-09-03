###################### LOAD LIBRARIES ##########################################

library(png)
library(spatstat)
library(imager)

###################### CONSTRUCTION TEMPERATURE COVARIATE ######################

# Definition of observation window in pixel units 
window = owin(xrange = c(-0.501, 861.501), yrange = c(-0.501, 741.501))

temp_june <- list()

for(aa in 2004:2024){
  
  # create spatstat object
  spatstat_mean <- im(mat = matrix(0, 741, 861), xcol = seq(0, 861, length.out = 861),
                     yrow = seq(0, 741, length.out = 741))
  # compute monthly averages
  for(day in 1:30){
    # download images from Canadian wildfire control website: change 06 in range 01 to 12 to choose month 
    url <- paste("https://cwfis.cfs.nrcan.gc.ca/data/maps/fwi_fbp/", aa,"/temp", aa, "06", sprintf("%02d", day), ".png", sep = "")
    destfile <- paste("temp", aa, "06", sprintf("%02d", day), ".png", sep = "")
    tryCatch({
      download.file(url, destfile, mode = "wb")
      img <-  readPNG(destfile)
      
      # resize for small images
      if(dim(img)[1] == 618){
        dim(img) <- c(dim(img)[1:2], 1, dim(img)[3])
        img_cimg <- as.cimg(img)
        img <- resize(img_cimg, size_x = 760, size_y = 880)
      }
      # Average RGB channels (grayscale conversion from colored image to heatmap)
      if (length(dim(img)) == 3) {
        img_gray <- round(0.2989 * img[,,1] + 0.5870 * img[,,2] + 0.1140 * img[,,3],2)
      } else {
        img_gray <- round(0.2989 * img[,,,1] + 0.5870 * img[,,,2] + 0.1140 * img[,,,3],2)
      }
      
      # remove margins 
      img_gray <- img_gray[10:750, 10:870]
      
      # intensity values corresponding to reported temperatures
      keep = c(0.3, 0.11,  0.63, 0.29, 0.52, 0.89)
      # intensity values corresponding to background information
      background = c(0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.77, 0.78, 0.79,
                     0.76, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85)
      # maintaining of all information not in the background 
      img_gray[img_gray %in% background] <- NA
      # img_gray[!(img_gray %in% keep)] <- NA
      
      # dropping parts of the plot with legends, titles, ecc...
      img_gray[1:90, 1:329] <- NA               # top left
      img_gray[1:60, 1:399] <- NA               # top left ext
      img_gray[1:330, 610:861] <- NA            # top right
      img_gray[1:130, 510:861] <- NA            # top right ext
      img_gray[630:741, 1:391] <- NA            # bottom left
      img_gray[640:741, 680:861] <- NA          # bottom right
      
      img_gray[731:741, 1:861] <- NA            # bottom 
      img_gray[1:741, 850:861] <- NA            # right
      img_gray[1:12, 1:861] <- NA               # top
      img_gray[1:741, 1:15] <- NA               # left
      
      
      # match the heat scale with the temperature values  
      img_gray[img_gray == 0.3] <- 35
      img_gray[img_gray == 0.11] <- 0
      img_gray[img_gray == 0.63] <- 22.5
      img_gray[img_gray == 0.29] <- 5
      img_gray[img_gray == 0.52] <- 15
      img_gray[img_gray == 0.89] <- 27.5
      
      xrange <- c(0, ncol(img_gray))
      yrange <- c(0, nrow(img_gray))
      
      # for border pixels, match the value to the neighboring pixels
      borders = setdiff(unique(as.vector(img_gray)),keep)[-1]
      indices <- which(matrix(img_gray %in% borders, nrow = 741), arr.ind = TRUE)
      img_gray[indices] <- apply(indices, 1, function(x){
            # browser()
            xi <- pmax(pmin(c(-8:8) + x[1], 741),1)
            yi <- pmax(pmin(c(-8:8) + x[2], 861),1)
            return(max(0, median(img_gray[xi,yi], na.rm = TRUE),na.rm = TRUE))
          })
      
      # store the information in a pixel image
      heat_img <- im(mat = img_gray[nrow(img_gray):1,], xcol = seq(xrange[1], xrange[2], length.out = ncol(img_gray)),
                          yrow = seq(yrange[1], yrange[2], length.out = nrow(img_gray)))
      
      # average daily temperatures over month
      spatstat_mean <- spatstat_mean + 1/30 * heat_img
      }, error = function(e) {
      # Handle error: print message and do something else
      message("this image was not available:", aa, day, e$message)
    })
  }
  print(aa)
  temp_june[[paste0(aa)]] <- spatstat_mean
}

lapply(temp_june, plot)
###################### IMPUTE MISSING VALUES TO UNIFY DOMAIN ###################
# mask of non NA values in either of the years 
l1_binary <- lapply(temp_june, function(x){im(!is.na(x$v), xcol = heat_img$xcol, yrow = heat_img$yrow)})
binary1 <- Reduce("|", l1_binary)
binary1 <- binary1 & !is.na(temp_june$`2004`$v)
plot(binary1)
mask_canada <- binary1
# function for imputing NA values
impute <- function(img){
  # set isolated pixels to NA
  # for(i in 2:(nrow(img)-1)){
  #   for(j in 2:(ncol(img)-1)){
  #     if(!is.na(img$v[i,j]) & sum(img$v[(i-1):(i+1), (j-1):(j+1)], na.rm = TRUE) == img$v[i,j]){
  #       img$v[i,j] = NA
  #     }
  #   }
  # }
  Img_vals <- img$v
  # set all values outside the mask to NA
  img$v[!mask_canada$v] <- NA
  # indexes where current image is NA, but we want to impute 
  # because other years have information
  mask_to_impute <- mask_canada$v & is.na(Img_vals) 
  indices <- which(mask_to_impute == TRUE, arr.ind = TRUE)
  img$v[indices] <- apply(indices, 1, function(x){
    # browser()
    xi <- pmax(pmin(c(-8:8) + x[1], 741),1)
    yi <- pmax(pmin(c(-8:8) + x[2], 861),1)
    return(max(0, median(img$v[xi,yi], na.rm = TRUE),na.rm = TRUE))
  })
  return(img)
}

temp_june <- lapply(temp_june, impute)

###################### PARAMETERS FOR HUMIDITY COVARIATE #######################

# keep = c(0.3, 0.11, 0.89, 0.63, 0.52, 0.66)
# img_gray[img_gray == 0.66] <- 1
# img_gray[img_gray == 0.3] <- 1
# img_gray[img_gray == 0.52] <- 0.3
# img_gray[img_gray == 0.89] <- 0.55
# img_gray[img_gray == 0.63] <- 0.75

###################### CONSTRUCTION PRECIPITAION COVARIATE #####################

precipitation_june <- list() 

for(aa in 2004:2025){
  spatstat_mean <- im(mat = matrix(0, 741, 861), xcol = seq(0, 861, length.out = 861),
                     yrow = seq(0, 741, length.out = 741))
  for(day in 1:30){
    url <- paste("https://cwfis.cfs.nrcan.gc.ca/data/maps/fwi_fbp/", aa,"/precip", aa, "06", sprintf("%02d", day), ".png", sep = "")
    destfile <- paste("precip", aa, "06", sprintf("%02d", day), ".png", sep = "")
    tryCatch({
      download.file(url, destfile, mode = "wb")
      img <-  readPNG(destfile)
      if(dim(img)[1] == 618){
        dim(img) <- c(dim(img)[1:2], 1, dim(img)[3])
        img_cimg <- as.cimg(img)
        img <- resize(img_cimg, size_x = 760, size_y = 880)
      }
      if (length(dim(img)) == 3) {
        # Average RGB channels
        img_gray <- round(0.2989 * img[,,1] + 0.5870 * img[,,2] + 0.1140 * img[,,3],2)
      } else {
        img_gray <- round(0.2989 * img[,,,1] + 0.5870 * img[,,,2] + 0.1140 * img[,,,3],2)
      }
      
      img_gray <- img_gray[10:750, 10:870]
      
      # intensity values corresponding to reported precipitation with color change in 2015
      if(aa > 2015){
        keep = c(0.52, 0.57, 0.11, 0.66, 0.49, 0.29)
      }else{
        keep = c(0.52, 0.57, 0.11, 0.63, 0.49, 0.29)

      }
      # intensity values corresponding to background information
      background = c(0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.77, 0.78, 0.79,
                     0.76, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85)
      # dropping of background pixels
      img_gray[img_gray %in% background] <- NA
    
      # dropping parts of the plot with legends, titles, ecc...
      img_gray[1:90, 1:329] <- NA               # top left
      img_gray[1:60, 1:399] <- NA               # top left ext
      img_gray[1:330, 610:861] <- NA            # top right
      img_gray[1:130, 510:861] <- NA            # top right ext
      img_gray[630:741, 1:391] <- NA            # bottom left
      img_gray[640:741, 680:861] <- NA          # bottom right
      
      img_gray[731:741, 1:861] <- NA            # bottom 
      img_gray[1:741, 850:861] <- NA            # right
      img_gray[1:12, 1:861] <- NA               # top
      img_gray[1:741, 1:15] <- NA               # left
      
      # for border pixels, match the value to the neighboring pixels
      borders = setdiff(unique(as.vector(img_gray)),keep)[-1]
      indices <- which(matrix(img_gray %in% borders, nrow = 741), arr.ind = TRUE)
      img_gray[indices] <- apply(indices, 1, function(x){
        xi <- pmax(pmin(c(-8:8) + x[1], 741),1)
        yi <- pmax(pmin(c(-8:8) + x[2], 861),1)
        return(max(0, median(img_gray[xi,yi], na.rm = TRUE),na.rm = TRUE))
      })
      
      # transform the heat scale, so to match the covariate values/colors to intensity 
      img_gray[img_gray == 0.66] <- 0 
      img_gray[img_gray == 0.63] <- 0 
      img_gray[img_gray == 0.57] <- 15 #0.6
      img_gray[img_gray == 0.52] <- 7.5 #0.4
      img_gray[img_gray == 0.29] <- 2.5 #0.2
      img_gray[img_gray == 0.11] <- 35 #1
      img_gray[img_gray == 0.49] <- 25 #0.8
      
      # transform in pixel image
      xrange <- c(0, ncol(img_gray))
      yrange <- c(0, nrow(img_gray))
      heat_img <- im(mat = img_gray[nrow(img_gray):1,], xcol = seq(xrange[1], xrange[2], length.out = ncol(img_gray)),
        yrow = seq(yrange[1], yrange[2], length.out = nrow(img_gray)))
      # plot(heat_img)
      
      # return monthly average
      spatstat_mean <- spatstat_mean + 1/30 * heat_img
      
      }, error = function(e) {
      # Handle error: print message and do something else
      message("this image was not available:", aa, day, e$message)
    })
  }
  print(aa)
  precipitation_june[[paste0(aa)]] <- spatstat_mean
}

lapply(precipitation_june, plot)

###################### IMPUTE MISSING VALUES TO UNIFY DOMAIN ###################
l2_binary <- lapply(precipitation_june, function(x){im(!is.na(x$v), xcol = heat_img$xcol, yrow = heat_img$yrow)})
binary2 <- Reduce("|", l2_binary)
binary2 <- binary2 & !is.na(precipitation_june$`2004`$v)
plot(binary2)
# the two mask are the same! 
sum(binary1 != binary2)

precipitation_june <- lapply(precipitation_june, impute)

###################### CONSTRUCTION WIND SPEED COVARIATE #######################

wind_speed_june <- list()

# function finding most represented value in a set, solving draws by randomization  
# like median, but avoiding values outside the color scheme to maintain 

median_in_set <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) return(NA)
  
  x_sorted <- sort(x)
  
  if (n %% 2 == 1) {
    # Odd length: take the middle value
    return(x_sorted[(n + 1) / 2])
  } else {
    # Even length: randomly pick one of the two middle values
    idx <- n / 2
    return(sample(x_sorted[c(idx, idx + 1)], 1))
  }
}

for(aa in 2004:2025){
  spatstat_mean <- im(mat = matrix(0, 741, 861), xcol = seq(0, 861, length.out = 861),
                      yrow = seq(0, 741, length.out = 741))
  for(day in 1:30){
    url <- paste("https://cwfis.cfs.nrcan.gc.ca/data/maps/fwi_fbp/", aa,"/ws", aa, "06", sprintf("%02d", day), ".png", sep = "")
    destfile <- paste("ws", aa, "06", sprintf("%02d", day), ".png", sep = "")
    tryCatch({
      download.file(url, destfile, mode = "wb")
      img <-  readPNG(destfile)
      
      if(dim(img)[1] == 618){
        dim(img) <- c(dim(img)[1:2], 1, dim(img)[3])
        img_cimg <- as.cimg(img)
        img <- resize(img_cimg, size_x = 760, size_y = 880)
      }
      
      if (length(dim(img)) == 3) {
        img_gray <- round(0.2989 * img[,,1] + 0.5870 * img[,,2] + 0.1140 * img[,,3],2)
      } else {
        img_gray <- round(0.2989 * img[,,,1] + 0.5870 * img[,,,2] + 0.1140 * img[,,,3],2)
      }
      
      img_gray <- img_gray[10:750, 10:870]
      
      # intensity values corresponding to reported wind speeds
      keep = c(0.3,0.11,0.63,0.29,0.52,0.89)
      # pixel values for background 
      background = c(0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.77, 0.78, 0.79,
                     0.76, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85)
      img_gray[img_gray %in% background] <- NA
      
      img_gray[1:90, 1:329] <- NA               # top left
      img_gray[1:60, 1:399] <- NA               # top left ext
      img_gray[1:330, 610:861] <- NA            # top right
      img_gray[1:130, 510:861] <- NA            # top right ext
      img_gray[630:741, 1:391] <- NA            # bottom left
      img_gray[640:741, 680:861] <- NA          # bottom right
      
      img_gray[731:741, 1:861] <- NA            # bottom 
      img_gray[1:741, 850:861] <- NA            # right
      img_gray[1:12, 1:861] <- NA               # top
      img_gray[1:741, 1:15] <- NA               # left
      
      # handling of all arrow and border pixels
      borders_arrows = setdiff(unique(as.vector(img_gray)),keep)[-1]
      indices <- which(matrix(img_gray %in% borders_arrows, nrow = 741), arr.ind = TRUE)
      # first set to NA, so they do not provide color information 
      img_gray[indices] <- NA
      # assign median neighbor color 
      img_gray[indices] <- apply(indices, 1, function(x){
        # browser()
        xi <- pmax(pmin(c(-6:6) + x[1], 741),1)
        yi <- pmax(pmin(c(-6:6) + x[2], 861),1)
        col = max(0, median_in_set(img_gray[xi,yi], na.rm = TRUE),na.rm = TRUE)
        return(col)
      })
      
      # transform the heat scale, so to match the covariate values/colors to intensity 
      img_gray[img_gray == 0.3] <- 30     #1
      img_gray[img_gray == 0.63] <- 22.5  #0.80
      img_gray[img_gray == 0.89] <- 17.5  #0.60
      img_gray[img_gray == 0.11] <- 2.5   #0
      img_gray[img_gray == 0.29] <- 7.5   #0.20
      img_gray[img_gray == 0.52] <- 12.5  #0.40

      heat_img <- im(mat = img_gray[nrow(img_gray):1,], xcol = seq(xrange[1], xrange[2], length.out = ncol(img_gray)),
                           yrow = seq(yrange[1], yrange[2], length.out = nrow(img_gray)))
      plot(heat_img)
      spatstat_mean <- spatstat_mean + 1/30 * heat_img

      
    }, error = function(e) {
      # Handle error: print message and do something else
      message("this image was not available:", aa, day, e$message)
    })
  }
  print(aa)
  wind_speed_june[[paste0(aa)]] <- spatstat_mean
}

lapply(wind_speed_june, function(x)plot(x, main = names(wind_speed_june)[which(sapply(wind_speed_june, identical, x))]))

###################### IMPUTE MISSING VALUES TO UNIFY DOMAIN ###################
l3_binary <- lapply(wind_speed_june, function(x){im(!is.na(x$v), xcol = heat_img$xcol, yrow = heat_img$yrow)})
binary3 <- Reduce("|", l3_binary)
binary3 <- binary3 & !is.na(wind_speed_june$`2004`$v)
plot(binary3)
# again mask is the same! 
sum(binary1 != binary3)

wind_speed_june <- lapply(wind_speed_june, impute)

###################### DEFINITION OF NON_BORDER PIXELS #########################
{
binary_image <- eval.im(!is.na(heat_img))
for(day in 1:30){
destfile <- paste("temp", 2024, "06", sprintf("%02d", day), ".png", sep = "")
img <-  readPNG(destfile)
img_gray <- round(0.2989 * img[,,1] + 0.5870 * img[,,2] + 0.1140 * img[,,3],2)
img_gray <- img_gray[10:750, 10:870]
keep = c(0.3, 0.11,  0.63, 0.29, 0.52, 0.89)
img_gray[!(img_gray %in% keep)] <- NA
img_gray[1:90, 1:329] <- NA               # top left
img_gray[1:60, 1:399] <- NA               # top left ext
img_gray[1:330, 610:861] <- NA            # top right
img_gray[1:130, 510:861] <- NA            # top right ext
img_gray[630:741, 1:391] <- NA            # bottom left
img_gray[640:741, 680:861] <- NA          # bottom right
img_gray[731:741, 1:861] <- NA            # bottom
img_gray[1:741, 850:861] <- NA            # right
img_gray[1:12, 1:861] <- NA               # top
img_gray[1:741, 1:15] <- NA               # left
xrange <- c(0, ncol(img_gray))
yrange <- c(0, nrow(img_gray))
heat_img <- im(mat = img_gray[nrow(img_gray):1,], xcol = seq(xrange[1], xrange[2], length.out = ncol(img_gray)),
               yrow = seq(yrange[1], yrange[2], length.out = nrow(img_gray)))
binary_image <- binary_image & eval.im(!is.na(heat_img))
}
}
plot(binary_image)
###################### CONSTRUCTION OF WILDFIRE POINT PROCESS LIST #############
ppp_june <- list()

for(aa in 2004:2024){
  points <- ppp(window = window)
  for(day in 1:30){
    url <- paste("https://cwfis.cfs.nrcan.gc.ca/data/maps/fireM3/", aa,"/tri", aa, "06", sprintf("%02d", day), ".png", sep = "")
    destfile <- paste("tri", aa, "06", sprintf("%02d", day), ".png", sep = "")
    tryCatch({
      # browser()
      download.file(url, destfile, mode = "wb")
      img <-  readPNG(destfile)
      img_gray <- round(0.2989 * img[,,1] + 0.5870 * img[,,2] + 0.1140 * img[,,3],2)
      img_gray <- img_gray[10:750, 10:870]
      
      # red and dark red pixels identify fire locations according to years
      keep = c(0.3, 0.36)
      img_gray[!img_gray %in% keep] <- NA
      
      # remove background points 
      img_gray[1:130, 1:231] <- NA
      img_gray[1:250, 590:861] <- NA
      img_gray[1:350, 650:861] <- NA
      img_gray[620:741, 1:281] <- NA
      img_gray[640:741, 700:861] <- NA
      img_gray[680:741, 1:400] <- NA
      img_gray[1:250, 590:861] <- NA
      
      # transforom the heatmap into pixel image 
      ppp_img <- im(mat = img_gray[nrow(img_gray):1,], xcol = seq(0, 861, length.out = 861),
                    yrow = seq(0, 741, length.out = 741))
      
      # exclude fires outside Canada and border pixels 
      ppp_img$v[!binary_image$v] <- NA
      
      # find indices of colored pixels (locations of fires)
      index <- which(!is.na(ppp_img$v), arr.ind = TRUE)
      
      # check how the image looks like in either of those locations 
      # {
      #   xi <- pmax(pmin(c(-8:8) + index[100,1], 741),1)
      #   yi <- pmax(pmin(c(-8:8) + index[100,2], 861),1)
      #   ppp_img$v[xi,yi]
      # }
      
      x_coords <- ppp_img$xcol[index[,2]]
      y_coords <- ppp_img$yrow[index[,1]]
      
      # set points in colored pixel locations
      points_ppp <- ppp(x = x_coords, y = y_coords, window = window)
      
      # each M3 fire location is identified by multiple pixels
      # identify close ppp objects so to maintain an observation each "square" 
      j = 1
      n = points_ppp$n
      dd <- pairdist(points_ppp)
      surv_points <- ppp(window = window)
      while(j<n){
        # browser()
        del <- which(dd[j,j:n]<9)
        surv_points <- superimpose(surv_points, points_ppp[j])
        j = j + length(del)
      }
      print(paste(day, aa, sep = "_"))
      points <- superimpose(points, surv_points)
      
    }, error = function(e) {
      # Handle error: print message and do something else
      message("this image was not available:", aa, day, e$message)
    })
  }
  ppp_june[[paste0(aa)]] <- points
}
# remove strange boundary points 
ppp_june <- lapply(ppp_june, function(x)return(x[coords(x)$y != 741]))

lapply(ppp_june, plot)


###################### DEFINE AND SUBSET ON THE BASIS OF REGIONS ###############

## division of Canada into regions by means of contour lines extrapolated from images ###
{
  mask_regions <- binary_image & mask_canada
  contour(mask_regions)
  # extract contour lines from 
  clines <- contourLines(mask_regions$xcol, mask_regions$yrow, t(mask_regions$v), levels = 0.5)
  # build list of polygons contouring each region
  poly_list = list()
  for(i in 1:length(clines)){
    coords <- clines[[i]]
    poly_list[[i]] <- tryCatch({
      owin(poly = list(x = rev(coords$x), y = rev(coords$y)))
    }, error = function(e) {
      # Handle error: print message and do something else
      message("Error caught: ", i, e$message)
      # Alternative execution
      owin(poly = list(x = (coords$x), y = (coords$y)))
    })
    plot(poly_list[[i]], lwd = 1, add = TRUE)
  }
  # merge all polygonal shapes
  merged_poly <- Reduce(union.owin, poly_list)
  merged_poly$xrange <- mask_regions$xrange
  merged_poly$yrange <- mask_regions$yrange
  plot(merged_poly, main = "Merged Polygon")
  
  # determine regions by means of enclosed area  
  areas <- sapply(poly_list, area.owin)
  tail(sort(areas),11)
  saska <- poly_list[[which(round(areas,2) == 14502.89)]]
  # 5389.02       10484.81        10497.84        14502.89          14562.04          14979.09    
  # LABRADOR      YUKON           NUNAVUT1        SASKATCHEWAN      MANITOBA          ALBERTA     
  # 18421.73      21086.41        21826.27        25205.75          32852.94
  # NUNAVUT2      BRIT COLUMBIA   ONTARIO         NW TERRITORIES    QUEBEC+ISL
  plot(saska)
  plot(temp_june$`2009`[saska, drop = FALSE])
  plot(ppp_june$`2009`[saska], add = TRUE)
  
  cov_saska <- temp_june$`2009`[saska, drop = FALSE]  
  saska_dom <- owin(saska$xrange, saska$yrange)
  cov_saska <- cov_saska[saska_dom]
  plot(cov_saska)
  cov_saska <- blur(cov_saska, sigma = 4, NA, bleed=FALSE, normalise = TRUE)
  plot(cov_saska)
}

###################### SMALL DOMAIN DATA RESTRICTION: saska ####################

saska_temp_june <- lapply(temp_june, function(x){y <- x[saska, drop = FALSE];
return(blur(y[saska_dom], sigma = 4, NA, bleed=FALSE, normalise = TRUE)) })
lapply(saska_temp_june, function(x)plot(x, main = names(saska_temp_june)[which(sapply(saska_temp_june, identical, x))]))
################################################################################
saska_precip_june <- lapply(precipitation_june, function(x){y <- x[saska, drop = FALSE];
return(blur(y[saska_dom], sigma = 4, NA, bleed=FALSE, normalise = TRUE)) })
lapply(saska_precip_june, function(x)plot(x, main = names(saska_precip_june)[which(sapply(saska_precip_june, identical, x))]))
################################################################################
saska_wind_june <- lapply(wind_speed_june, function(x){y <- x[saska, drop = FALSE];
return(blur(y[saska_dom], sigma = 4, NA, bleed=FALSE, normalise = TRUE)) })
lapply(saska_wind_june, function(x)plot(x, main = names(saska_wind_june)[which(sapply(saska_wind_june, identical, x))]))
################################################################################
saska_ppp_june <- lapply(ppp_june, function(x){x[saska]})
################################################################################

# save("areas", "saska_ppp_june", "saska_precip_june", "saska_temp_june" , 
#      "saska_wind_june" ,  "saska","precipitation_june", 
#      "wind_speed_june" ,"temp_june" , "ppp_june"  , "poly_list", 
#      "mask_canada", "window", "mask_regions", "saska_dom" ,file = "data_june_saska.RData")

save("saska_precip_june", 
     "saska_wind_june",
     "saska_temp_june",
     "saska_ppp_june",
     "mask_canada", "window", "saska", "saska_dom",
     file = "saska_data_june_orig_scale.RData")
