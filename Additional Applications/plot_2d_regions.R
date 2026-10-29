###############################################################################
# plot_2d_regions.R -- figures of the 2D (temperature, precipitation) fits for Saskatchewan and
# British Columbia, as drawn for Ontario in ../Application/plot_canada.R.
#
# Reads only (no sampling, nothing written outside --out):
#   fits:  Results/Canada_Regions/Fits/<Province>_temp_prec   (run_additional_regions.R
#          --fits=saskatchewan:temp_prec,british_columbia:temp_prec)
#   data:  ../Application/Data
# Output (default Results/Canada_Regions/Figures_2D):
#   saska_2d_estim.pdf, bt_columbia_2d_estim.pdf  -- as ontario_2d_estim.pdf: (1) covariates at the
#       events with the convex hull of the observed (temperature, precipitation) pairs; (2) GP
#       posterior mean on the whole square, hull in white; (3) averaged CV kernel estimate inside
#       the hull (same code as plot_canada.R: ../R/competitors.R, kernel_competitor/kernel_on_grid).
#   saska_2d_spat.pdf, bt_columbia_2d_spat.pdf    -- plug-in spatial maps rho_hat(Z^(i)(x)) for 2013,
#       2015, 2021, one colour scale per figure (as plot_spatial_maps.R).
# The helper functions below are copied from ../Application/plot_canada.R and plot_spatial_maps.R
# (the shared code in ../R/ is sourced unchanged).
#
# Run from 'Data Analysis/Additional Applications':
#   Rscript plot_2d_regions.R
#   Rscript plot_2d_regions.R --regions=saskatchewan --kernel=FALSE
# From RStudio:  rs("plot_2d_regions.R")   (rs() as in README-ADDITIONAL-APPLICATIONS.md)
###############################################################################
local({ need <- c("Matrix", "posterior", "ggplot2", "patchwork")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(patchwork) })
if (!file.exists(file.path("..", "R", "canada_posterior.R")) || !dir.exists(file.path("..", "Application", "Data")))
  stop("Run this script from 'Data Analysis/Additional Applications'")
source(file.path("..", "R", "gp_sampler.R"))
source(file.path("..", "R", "canada_data.R"))
source(file.path("..", "R", "canada_posterior.R"))
source(file.path("..", "R", "competitors.R"))
if (!exists("rel_L2")) rel_L2 <- function(est, truth) sqrt(mean((est - truth)^2)) / sqrt(mean(truth^2))
CANADA_DATA_DIR <- file.path("..", "Application", "Data")          # read-only
FIT_DIR <- file.path("Results", "Canada_Regions", "Fits")

opt <- list(out = file.path("Results", "Canada_Regions", "Figures_2D"),
            regions = "saskatchewan,british_columbia", kernel = "TRUE")
for (a in commandArgs(TRUE)) {
  k <- sub("^--([a-z_]+)=.*$", "\\1", a)
  if (!grepl("^--[a-z_]+=", a) || !k %in% names(opt)) stop("Unknown option ", a)
  opt[[k]] <- sub("^--[a-z_]+=", "", a)
}
regions <- strsplit(opt$regions, ",")[[1]]
if (!all(regions %in% c("saskatchewan", "british_columbia")))
  stop("--regions must be a subset of saskatchewan,british_columbia")
use_kernel <- as.logical(opt$kernel)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(opt$out)
if (startsWith(out, normalizePath(file.path("..", "Application")))) stop("--out must not lie inside ../Application")

PAL <- c("darkblue", "deeppink3", "darkgoldenrod1")
FILE_PREFIX <- c(saskatchewan = "saska", british_columbia = "bt_columbia")
SPATIAL_YEARS <- c(2013, 2015, 2021)

# ------------------------------------------------------------------ helpers (from plot_canada.R)
.regions <- new.env(); .pools <- new.env(); .kern <- new.env()
get_region <- function(region) {
  if (is.null(.regions[[region]])) .regions[[region]] <- load_canada_region(region)
  .regions[[region]]
}
get_pool <- function(region, model) {
  key <- paste(region, model, sep = "__")
  if (!is.null(.pools[[key]])) return(.pools[[key]])
  fd <- file.path(FIT_DIR, fit_folder(key)); st <- file.path(fd, "status.rds")
  if (!file.exists(st)) { message("Skipping ", region, ":", model, " (no assessed fit in ", fd, ")"); return(NULL) }
  s <- readRDS(st)
  if (!isTRUE(s$final)) { message("Skipping ", region, ":", model, " (fit not final)"); return(NULL) }
  if (!isTRUE(s$passed)) message("Note: ", region, ":", model, " stopped at the sweep cap without meeting all convergence criteria")
  .pools[[key]] <- load_pool(fd)
  .pools[[key]]
}
run_fig <- function(expr) tryCatch(expr, error = function(e) message("Figure failed: ", conditionMessage(e)))
save_plot <- function(p, name, width, height) {
  f <- file.path(out, name); ggsave(f, p, width = width, height = height); cat("Written:", f, "\n"); flush.console()
}
axis_breaks <- c(0, 0.25, 0.5, 0.75, 1)
orig_labels <- function(R, cv) round(canada_backtransform(R, cv, axis_breaks), 2)
event_u <- function(R, cv) unlist(lapply(seq_along(R$years), function(i) R$U[[cv]][R$events[[i]]$pix, i]))
kernel_cv <- function(region, covs) {
  key <- paste(c(region, covs), collapse = "__")
  if (is.null(.kern[[key]])) {
    kc <- kernel_competitor(canada_reps(get_region(region), covs), length(covs),
                            rho0 = function(z) rep(1, nrow(z)))
    cat(sprintf("Kernel estimate %s (%s): CV bandwidth %s%s\n", region, paste(covs, collapse = ", "),
                paste(signif(kc$cv$h, 3), collapse = " x "),
                if (isTRUE(kc$cv$at_edge)) " (AT THE EDGE OF THE GRID)" else "")); flush.console()
    .kern[[key]] <- kc
  }
  .kern[[key]]
}
kernel_2d <- function(region, cv1, cv2, z) kernel_on_grid(kernel_cv(region, c(cv1, cv2)), z)[, 2]
inside_convex_hull <- function(P, z) {             # P: hull vertices in chull order
  n <- nrow(P); s <- matrix(0, nrow(z), n)
  for (k in seq_len(n)) {
    a <- P[k, ]; b <- P[if (k == n) 1 else k + 1, ]
    s[, k] <- (b[1] - a[1]) * (z[, 2] - a[2]) - (b[2] - a[2]) * (z[, 1] - a[1])
  }
  apply(s, 1, function(v) all(v <= 1e-12) || all(v >= -1e-12))
}

# ------------------------------------------------------------------ 2D estimate (plot_canada.R::plot_2d)
plot_2d <- function(region) {
  pool <- get_pool(region, "temp_prec"); if (is.null(pool)) return(invisible(NULL))
  R <- get_region(region); g <- seq(0, 1, length.out = 200)
  z <- as.matrix(expand.grid(z1 = g, z2 = g))
  P <- cbind(as.vector(R$U$temp), as.vector(R$U$prec)); P <- P[chull(P), ]
  inside <- inside_convex_hull(P, z)
  # GP posterior mean on the whole square (outside the hull it is driven by the prior and the
  # GP correlation, not by data); the kernel estimate exists only inside the hull.
  df <- data.frame(x = z[, 1], y = z[, 2], value = posterior_rho(pool, z)$mean)
  axes <- list(coord_fixed(), theme_minimal(), labs(x = CANADA_LABELS[["temp"]], y = CANADA_LABELS[["prec"]]),
               scale_x_continuous(breaks = axis_breaks, labels = orig_labels(R, "temp")),
               scale_y_continuous(breaks = axis_breaks, labels = orig_labels(R, "prec")))
  hull <- data.frame(x = c(P[, 1], P[1, 1]), y = c(P[, 2], P[1, 2]))
  hull_line <- geom_path(data = hull, aes(x = x, y = y), inherit.aes = FALSE, colour = "white", linewidth = 0.6)
  p1 <- ggplot() + geom_path(data = hull, aes(x = x, y = y), colour = "grey50", linewidth = 0.5) +
    geom_point(data = data.frame(x = event_u(R, "temp"), y = event_u(R, "prec")), aes(x = x, y = y),
               shape = 19, size = 1, colour = "black") + axes
  p2 <- ggplot(df, aes(x = x, y = y, fill = value)) + geom_raster() + hull_line +
    scale_fill_gradientn(colours = PAL, na.value = "white", name = expression(rho(z))) + axes
  plots <- list(p1, p2)
  if (use_kernel) {
    df$kern <- ifelse(inside, kernel_2d(region, "temp", "prec", z), NA)
    plots[[3]] <- ggplot(df, aes(x = x, y = y, fill = kern)) + geom_raster() + hull_line +
      scale_fill_gradientn(colours = PAL, na.value = "white", name = expression(rho(z))) + axes
  }
  save_plot(wrap_plots(plots, nrow = 1), sprintf("%s_2d_estim.pdf", FILE_PREFIX[[region]]), 4.4 * length(plots), 3.8)
}

# ------------------------------------------------------------------ spatial maps, common scale (plot_spatial_maps.R)
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
spatial_figure <- function(region, model, years, name, ncol = 3) {
  pool <- get_pool(region, model); if (is.null(pool)) return(invisible(NULL))
  R <- get_region(region)
  vals <- lapply(years, function(y) plugin_intensity(pool, R, y))
  lim <- range(unlist(vals))
  plots <- Map(function(v, y) map_plot(R, v, y, lim), vals, years)
  nrow <- ceiling(length(plots) / ncol)
  save_plot(wrap_plots(plots, ncol = ncol, guides = "collect"), name, 4.6 * ncol + 0.8, 4.2 * nrow)
  cat(sprintf("  colour scale %.4g - %.4g\n", lim[1], lim[2]))
}

# ------------------------------------------------------------------ figures
for (region in regions) {
  run_fig(plot_2d(region))
  run_fig(spatial_figure(region, "temp_prec", SPATIAL_YEARS, sprintf("%s_2d_spat.pdf", FILE_PREFIX[[region]])))
  rm(list = ls(.pools), envir = .pools); gc(FALSE)
}
