###############################################################################
# plot_spatial_maps.R -- yearly spatial plug-in intensity maps with ONE colour scale per figure
#
# Redraws every spatial yearly figure of the paper (Ontario, main text and Supplement) and of
# Supplement S5 (Saskatchewan, British Columbia) with a common intensity scale across the years of
# each figure and a single legend, so that the panels of a figure can be compared directly.
# (plot_canada.R in ../Application uses a separate scale for each panel.)
#
# Reads only (no sampling, nothing written outside --out):
#   Ontario:                  ../Application/Results/Canada/Fits       (fits of the paper)
#   Saskatchewan, Br. Col.:   Results/Canada_Regions/Fits              (run_additional_regions.R)
#   data:                     ../Application/Data
# The posterior summaries are computed exactly as in plot_canada.R (plug-in posterior mean
# rho_hat(Z^(i)(x)) at every valid pixel, ../R/canada_posterior.R), with the same map layout,
# palette and geographical axes; only the colour scale is shared.
#
# Run from 'Data Analysis/Additional Applications':
#   Rscript plot_spatial_maps.R                              # all figures, Results/Spatial_Maps
#   Rscript plot_spatial_maps.R --regions=saskatchewan,british_columbia
#   Rscript plot_spatial_maps.R --out=Results/Spatial_Maps_test
# Figures of the 3D fits: the main-text years 2013, 2015, 2021 for every province (ontario_3d.pdf,
# saska_3d.pdf, BT_3d.pdf) and all years 2004-2022 in one figure (<same>_all.pdf).
# Output: <out>/<file>.pdf, named as in the manuscript's figures/ folder, and <out>/scales.csv
# (colour-scale range of each figure and the maximum of each yearly map).
###############################################################################
local({ need <- c("Matrix", "ggplot2", "patchwork")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(patchwork) })
if (!file.exists(file.path("..", "R", "canada_posterior.R")) || !dir.exists(file.path("..", "Application", "Data")))
  stop("Run this script from 'Data Analysis/Additional Applications'")
source(file.path("..", "R", "gp_sampler.R"))
source(file.path("..", "R", "canada_data.R"))
source(file.path("..", "R", "canada_posterior.R"))
CANADA_DATA_DIR <- file.path("..", "Application", "Data")          # read-only

FIT_DIRS <- c(ontario = file.path("..", "Application", "Results", "Canada", "Fits"),
              saskatchewan = file.path("Results", "Canada_Regions", "Fits"),
              british_columbia = file.path("Results", "Canada_Regions", "Fits"))

opt <- list(out = file.path("Results", "Spatial_Maps"), regions = "all")
for (a in commandArgs(TRUE)) {
  k <- sub("^--([a-z_]+)=.*$", "\\1", a)
  if (!grepl("^--[a-z_]+=", a) || !k %in% names(opt)) stop("Unknown option ", a)
  opt[[k]] <- sub("^--[a-z_]+=", "", a)
}
regions <- if (identical(opt$regions, "all")) names(FIT_DIRS) else strsplit(opt$regions, ",")[[1]]
if (!all(regions %in% names(FIT_DIRS))) stop("--regions must be a subset of ", paste(names(FIT_DIRS), collapse = ","))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(opt$out)
if (startsWith(out, normalizePath(file.path("..", "Application")))) stop("--out must not lie inside ../Application")

PAL <- c("darkblue", "deeppink3", "darkgoldenrod1")                # as in plot_canada.R

# ------------------------------------------------------------------ data and fits (cached)
.regions <- new.env(); .pools <- new.env()
get_region <- function(region) {
  if (is.null(.regions[[region]])) .regions[[region]] <- load_canada_region(region)
  .regions[[region]]
}
get_pool <- function(region, model) {
  key <- paste(region, model, sep = "__")
  if (!is.null(.pools[[key]])) return(.pools[[key]])
  fd <- file.path(FIT_DIRS[[region]], fit_folder(key)); st <- file.path(fd, "status.rds")
  if (!file.exists(st)) { message("Skipping ", region, ":", model, " (no assessed fit in ", fd, ")"); return(NULL) }
  s <- readRDS(st)
  if (!isTRUE(s$final)) { message("Skipping ", region, ":", model, " (fit not final)"); return(NULL) }
  if (!isTRUE(s$passed)) message("Note: ", region, ":", model, " stopped at the sweep cap without meeting all convergence criteria")
  .pools[[key]] <- load_pool(fd)
  .pools[[key]]
}
drop_pools <- function() { rm(list = ls(.pools), envir = .pools); gc(FALSE) }

# ------------------------------------------------------------------ one map (plot_canada.R layout)
map_plot <- function(R, values, year, limits) {
  geo <- CANADA_GEO[[R$region]]
  xb <- seq(geo$lon_px[1], geo$lon_px[2], length.out = geo$n_xbreaks)
  yb <- seq(geo$lat_px[1], geo$lat_px[2], length.out = 5)
  lon <- function(x) geo$lon[1] + (x - geo$lon_px[1]) / diff(geo$lon_px) * diff(geo$lon)
  lat <- function(y) geo$lat[1] + (y - geo$lat_px[1]) / diff(geo$lat_px) * diff(geo$lat)
  df <- data.frame(x = R$xcol[R$valid_col], y = R$yrow[R$valid_row], value = values)
  ev <- R$events[[as.character(year)]]
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile(width = diff(R$xcol[1:2]), height = diff(R$yrow[1:2])) +
    scale_fill_gradientn(colours = PAL, limits = limits, na.value = "white", name = expression(lambda(x))) +
    coord_fixed() + theme_minimal() +
    annotate("rect", xmin = min(R$xcol), xmax = max(R$xcol), ymin = min(R$yrow), ymax = max(R$yrow),
             fill = NA, colour = "black", linewidth = 0.3) +
    scale_x_continuous(breaks = xb, labels = round(lon(xb), 2), name = "Longitude") +
    scale_y_continuous(breaks = yb, labels = round(lat(yb), 2), name = "Latitude") +
    geom_point(data = ev, aes(x = x, y = y), inherit.aes = FALSE, shape = 16, colour = "turquoise1", size = 0.9) +
    ggtitle(as.character(year))
}

# ------------------------------------------------------------------ one figure, common scale
scales <- list()
spatial_figure <- function(region, model, years, name, ncol = 3) {
  pool <- get_pool(region, model); if (is.null(pool)) return(invisible(NULL))
  R <- get_region(region)
  vals <- lapply(years, function(y) plugin_intensity(pool, R, y))
  lim <- range(unlist(vals))
  plots <- Map(function(v, y) map_plot(R, v, y, lim), vals, years)
  nrow <- ceiling(length(plots) / ncol)
  f <- file.path(out, name)
  ggsave(f, wrap_plots(plots, ncol = ncol, guides = "collect"), width = 4.6 * ncol + 0.8, height = 4.2 * nrow)
  cat(sprintf("Written: %s  (scale %.4g - %.4g)\n", f, lim[1], lim[2])); flush.console()
  scales[[length(scales) + 1]] <<- data.frame(file = name, region = region, model = model,
    years = paste(years, collapse = " "), scale_min = lim[1], scale_max = lim[2],
    yearly_max = paste(signif(vapply(vals, max, 0), 3), collapse = " "))
}
run_fig <- function(expr) tryCatch(expr, error = function(e) message("Figure failed: ", conditionMessage(e)))

# ------------------------------------------------------------------ figures (as in plot_canada.R)
# Years shown for Saskatchewan and British Columbia: those of the main-text Ontario maps (Figs. 6
# and 8), in one row as there (previously nine years, 2006-2022, in a 3 x 3 grid).
main_years <- c(2013, 2015, 2021)
# All 19 years of the data, in one figure per province (3D fit, one common scale, 5 columns).
all_years <- 2004:2022
if ("ontario" %in% regions) {
  run_fig(spatial_figure("ontario", "temp", c(2013, 2015, 2021), "ontario_temp.pdf"))          # Fig. 6
  run_fig(spatial_figure("ontario", "full", c(2013, 2015, 2021), "ontario_3d.pdf"))            # Fig. 8
  run_fig(spatial_figure("ontario", "temp_prec", c(2013, 2015, 2021), "ontario_2d_spat.pdf"))  # Suppl. 2D maps
  run_fig(spatial_figure("ontario", "full", c(2006, 2008, 2010, 2016, 2018, 2022), "ontario_3d_add.pdf"))
  run_fig(spatial_figure("ontario", "full", 2004:2009, "ontario_3d_1.pdf"))
  run_fig(spatial_figure("ontario", "full", c(2010, 2011, 2012, 2014, 2016, 2017), "ontario_3d_2.pdf"))
  run_fig(spatial_figure("ontario", "full", c(2018, 2019, 2020, 2022), "ontario_3d_3.pdf", ncol = 2))
  run_fig(spatial_figure("ontario", "full", all_years, "ontario_3d_all.pdf", ncol = 5))
  drop_pools()
}
if ("saskatchewan" %in% regions) {
  run_fig(spatial_figure("saskatchewan", "full", main_years, "saska_3d.pdf"))
  run_fig(spatial_figure("saskatchewan", "full", all_years, "saska_3d_all.pdf", ncol = 5))
  drop_pools()
}
if ("british_columbia" %in% regions) {
  run_fig(spatial_figure("british_columbia", "full", main_years, "BT_3d.pdf"))
  run_fig(spatial_figure("british_columbia", "full", all_years, "BT_3d_all.pdf", ncol = 5))
  drop_pools()
}

if (length(scales)) {
  tab <- do.call(rbind, scales)
  sf <- file.path(out, "scales.csv")
  if (file.exists(sf)) { old <- read.csv(sf); tab <- rbind(old[!old$file %in% tab$file, ], tab) }
  write.csv(tab, sf, row.names = FALSE)
  cat("Written:", sf, "\n")
}
