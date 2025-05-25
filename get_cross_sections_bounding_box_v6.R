library(lidR)
library(data.table)
library(png)
require(rgl)
require(sp)

# --- Slice Extraction ---
extract_slice <- function(las, center, angle_deg, thickness = 0.5) {
  angle_rad <- angle_deg * pi / 180
  dir_x <- cos(angle_rad)
  dir_y <- sin(angle_rad)
  
  dx <- las@data$X - center[1]
  dy <- las@data$Y - center[2]
  
  proj_along <- dx * dir_x + dy * dir_y  # x'
  proj_perp <- dx * (-dir_y) + dy * dir_x  # depth
  
  keep <- abs(proj_perp) <= (thickness / 2)
  depth_vals <- proj_perp[keep]
  depth_norm <- (depth_vals + thickness / 2) / thickness
  
  colors <- colorRampPalette(c("#f0f0f0", "#1f78b4"))(100)
  color_idx <- pmin(pmax(round(depth_norm * 99) + 1, 1), 100)
  
  data.frame(
    x_proj = proj_along[keep],
    z = las@data$Z[keep],
    X = las@data$X[keep],
    Y = las@data$Y[keep],
    Z = las@data$Z[keep],
    col = colors[color_idx]
  )
}




# --- Parameters ---

# select next tree (for each tree increase fnr by 1)
fnr=1

path_to_file <- list.files("D:/Tree_PCs", pattern=".laz", full.names = T)
filen <- basename(path_to_file[fnr])

las1 <- readLAS(path_to_file[fnr])
las <- filter_poi(las1, runif(npoints(las1)) < 0.25)

# check whether las file contains data
if (is.empty(las)) stop("LAS file is empty")

# calculate some base metrics
# center of point clouds
center <- c(mean(las@data$X), mean(las@data$Y))
# determine/define minimum and maximum height
zmin <- min(las@data$Z)
zmax <- zmin + 30

####
#### these settings matter
#### 

# define thickness of transect views
slice_thickness <- 1.25

# define size of and steps of transects cuts through point cloud (in degree steps)
# more angles means more slices to interpret per tree => will potentially be more accurate but also
# more work 
steps = 36
angle <- seq(0,360-steps,steps)  # degrees



# create empty list to store resulting digitized coordinates
coords_tree <- list()

# loop through angles and label data

for (i in 1:length(angle)){
  
  # --- Step 1: Create Slice and Save PNG ---
  slice_df <- extract_slice(las, center, angle[i], slice_thickness)
  xlim_range <- range(slice_df$x_proj)
  ylim_range <- c(zmin, zmax)
  
  # define folder to save temporary images  
  setwd("D:/Tree_PCs_hulls")
  
  png("slice_plot.png", width = 800, height = 1200)
  # Margins area
  #par(oma=c(0,0,0,0)) # all sides have 3 lines of space
  #par(mar=c(0,0,0,0))
  plot(slice_df$z ~ slice_df$x_proj,
       col = slice_df$col, pch = 16, cex = 0.3,
       xlab = "x'", ylab = "Z", xlim = xlim_range, ylim = ylim_range)
  abline(v=0)
  abline(v=5)
  abline(v=-5)
  
  dev.off()
  
  # --- Step 2: Display PNG in Correct Coordinate System ---
  img <- readPNG("slice_plot.png")
  
  plot(NA, xlim = xlim_range, ylim = ylim_range,
       xlab = "x'", ylab = "Z", main = paste("Click crown points - Angle", angle[i]))
  rasterImage(img, xlim_range[1], ylim_range[1], xlim_range[2], ylim_range[2])
  
  message("Click as many points as you want. ESC to stop.")
  clicks <- locator(type = "p", col = "red", pch = 16)
  
  # --- Step 3: Convert Clicks to 3D Coordinates ---
  if (!is.null(clicks)) {
    angle_rad <- angle[i] * pi / 180
    dir_x <- cos(angle_rad)
    dir_y <- sin(angle_rad)
    
    x_proj_clicked <- clicks$x
    z_clicked <- clicks$y
    
    coords_3D <- data.frame(
      X = center[1] + x_proj_clicked * dir_x,
      Y = center[2] + x_proj_clicked * dir_y,
      Z = z_clicked
    )
    
    print(coords_3D)
    coords_tree[[i]] <- coords_3D
    print(paste0("Completed ", i, "th transect view"))
  }
  
  
}


# merge all coordinates to one file
coords_fin <- do.call(rbind,coords_tree)
head(coords_fin)

# convert to las file
hull <- LAS(coords_fin)

# plot las f ile
plot(hull)

# save hull fil
writeLAS(hull, paste0("hull_" , filen))a



# create 2D convex hull
hull_coords <- hull@data[, c("X", "Y")]
hull_indices <- chull(hull_coords)
hull_coords2 <- hull_coords[c(hull_indices, hull_indices[1]), ]  # Close the polygon

# Create a SpatialPolygons object from the convex hull
hull_polygon <- SpatialPolygons(list(Polygons(list(Polygon(hull_coords2)), ID = "1")))
#hull_polygon2 <- gBuffer(hull_polygon, width = 0.1)

# Clip the first point cloud using the convex hull polygon
proj4string(hull_polygon) <- crs(las)

clipped <- clip_roi(las1, hull_polygon)

# write to file
writeLAS(clipped, paste0("extr_tree_", filen))

