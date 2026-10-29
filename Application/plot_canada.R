###############################################################################
# plot_canada.R -- figures of the wildfire application from the revised fits.
# Run from this folder after run_canada.R:
#   Rscript plot_canada.R                      # uses Results/Canada
#   Rscript plot_canada.R --out=Results/Smoke --kernel=FALSE
# Output: <out>/Figures/, named as the corresponding files in the manuscript's figures/
# folder. Fits that are missing are skipped; fits that did not pass the convergence
# criteria are plotted with a warning.
#
# All GP summaries are computed exactly from the stored node draws (R/canada_posterior.R):
# no post-hoc smoothing (the old Smooth(), blur() and loess steps are no longer needed with
# linear interpolation). The kernel estimates (dashed curves in 1D, third panel in 2D) are the
# averaged ratio-form kernel estimator of the simulation study with leave-one-year-out
# likelihood-CV bandwidth, computed by the same function (../R/competitors.R,
# kernel_competitor, with its default extended bandwidth grid: 1000 bins and bandwidths in
# [0.003, 1] in 1D, 160 x 160 bins and [0.008, 2]^2 in 2D, on the ECDF scale): plain average of
# the yearly estimates, no further smoothing. (Changed 11 October 2026: the old
# script weighted the yearly rhohat estimates by their numbers of events, which biases the
# average upwards, and smoothed it with smooth.spline.)
###############################################################################

# Required packages (install once: see ../README-DATA-ANALYSIS.md, section Software)
local({ need <- c("Matrix", "posterior", "ggplot2", "patchwork")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "), ". Install with\n  install.packages(c(",
                         paste0('"', miss, '"', collapse = ", "), "))", call. = FALSE) })
suppressPackageStartupMessages({ library(Matrix); library(ggplot2); library(patchwork) })
source("../R/gp_sampler.R"); source("../R/canada_data.R"); source("../R/canada_posterior.R")

opt <- list(out = file.path("Results", "Canada"), kernel = "TRUE")
for (a in commandArgs(TRUE)) {
  k <- sub("^--([a-z_]+)=.*$", "\\1", a); if (!k %in% names(opt)) stop("Unknown option ", a)
  opt[[k]] <- sub("^--[a-z_]+=", "", a)
}
use_kernel <- as.logical(opt$kernel)
figdir <- file.path(opt$out, "Figures"); dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
PAL <- c("darkblue", "deeppink3", "darkgoldenrod1")
FILE_PREFIX <- c(ontario = "ontario", saskatchewan = "saska", british_columbia = "bt_columbia")
SPATIAL_3D_FILE <- c(ontario = "ontario_3d", saskatchewan = "saska_3d", british_columbia = "BT_3d")

# ------------------------------------------------------------------ helpers
.regions <- new.env()
get_region <- function(region) {
  if (is.null(.regions[[region]])) .regions[[region]] <- load_canada_region(region)
  .regions[[region]]
}
.pools <- new.env()
get_pool <- function(region, model) {
  key <- paste(region, model, sep = "__")
  if (!is.null(.pools[[key]])) return(.pools[[key]])
  fd <- file.path(opt$out, "Fits", fit_folder(key)); st <- file.path(fd, "status.rds")
  if (!file.exists(st)) { message("Skipping ", region, ":", model, " (no assessed fit)"); return(NULL) }
  s <- readRDS(st)
  if (!isTRUE(s$final)) { message("Skipping ", region, ":", model, " (still running or being extended)"); return(NULL) }
  if (!isTRUE(s$passed)) warning(region, ":", model, " has not passed the convergence criteria", call. = FALSE)
  .pools[[key]] <- load_pool(fd)
  .pools[[key]]
}
run_fig <- function(expr) tryCatch(expr, error = function(e) message("Figure failed: ", conditionMessage(e)))
save_plot <- function(p, name, width, height) {
  f <- file.path(figdir, name); ggsave(f, p, width = width, height = height); cat("Written:", f, "\n")
}
axis_breaks <- c(0, 0.25, 0.5, 0.75, 1)
orig_labels <- function(R, cv) round(canada_backtransform(R, cv, axis_breaks), 2)
event_covariate <- function(R, cv) unlist(lapply(seq_along(R$years), function(i) R$raw[[cv]][R$events[[i]]$pix, i]))
event_u <- function(R, cv) unlist(lapply(seq_along(R$years), function(i) R$U[[cv]][R$events[[i]]$pix, i]))

# Averaged ratio-form kernel estimator with leave-one-year-out likelihood-CV bandwidth: the
# simulation study's implementation, unchanged. kernel_competitor() also selects an oracle
# bandwidth from the true intensity, which is unknown here: a constant placeholder is passed and
# only the CV selection is used.
source(file.path("..", "R", "competitors.R"))
if (!exists("rel_L2")) rel_L2 <- function(est, truth) sqrt(mean((est - truth)^2)) / sqrt(mean(truth^2))
.kern <- new.env()
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
# CV kernel estimate at points u (1D) or z (2D): interpolation between bin centres (competitors.R)
kernel_1d <- function(region, cv, u) kernel_on_grid(kernel_cv(region, cv), matrix(u))[, 2]
kernel_2d <- function(region, cv1, cv2, z) kernel_on_grid(kernel_cv(region, c(cv1, cv2)), z)[, 2]
inside_convex_hull <- function(P, z) {             # P: hull vertices in chull order
  n <- nrow(P); s <- matrix(0, nrow(z), n)
  for (k in seq_len(n)) {
    a <- P[k, ]; b <- P[if (k == n) 1 else k + 1, ]
    s[, k] <- (b[1] - a[1]) * (z[, 2] - a[2]) - (b[2] - a[2]) * (z[, 1] - a[1])
  }
  apply(s, 1, function(v) all(v <= 1e-12) || all(v >= -1e-12))
}

# ------------------------------------------------------------------ 1D curves
plot_1d <- function(region) {
  R <- get_region(region); plots <- list()
  for (cv in c("temp", "prec", "wind")) {
    pool <- get_pool(region, cv); if (is.null(pool)) next
    u <- seq(0, 1, length.out = 201)
    pr <- posterior_rho(pool, matrix(u), probs = c(0.025, 0.975))
    df <- data.frame(x = canada_backtransform(R, cv, u), mean = pr$mean,
                     lo = pr$quantiles[, 1], hi = pr$quantiles[, 2])
    p <- ggplot(df, aes(x = x)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "gray20", alpha = 0.2) +
      geom_line(aes(y = mean), colour = "gray30", linewidth = 1) +
      geom_rug(data = data.frame(x = event_covariate(R, cv)), aes(x = x), sides = "b", inherit.aes = FALSE) +
      theme_minimal(base_family = "sans") +
      theme(panel.grid.major = element_line(colour = "gray80"), panel.grid.minor = element_line(colour = "gray90"),
            text = element_text(colour = "black"), axis.text = element_text(size = 11)) +
      labs(x = CANADA_LABELS[[cv]], y = expression(rho(z)))
    if (use_kernel) {
      df$kern <- kernel_1d(region, cv, u)
      p <- p + geom_line(data = df, aes(y = kern), colour = "gray20", linewidth = 1, linetype = "dashed")
    }
    plots[[cv]] <- p
  }
  if (length(plots)) save_plot(wrap_plots(plots, nrow = 1), sprintf("%s_1d%s.pdf", FILE_PREFIX[[region]],
                               if (region == "ontario") "_smoth" else ""), 4 * length(plots), 3.6)
}

# ------------------------------------------------------------------ spatial maps
map_plot <- function(R, values, year) {
  geo <- CANADA_GEO[[R$region]]
  xb <- seq(geo$lon_px[1], geo$lon_px[2], length.out = geo$n_xbreaks)
  yb <- seq(geo$lat_px[1], geo$lat_px[2], length.out = 5)
  lon <- function(x) geo$lon[1] + (x - geo$lon_px[1]) / diff(geo$lon_px) * diff(geo$lon)
  lat <- function(y) geo$lat[1] + (y - geo$lat_px[1]) / diff(geo$lat_px) * diff(geo$lat)
  df <- data.frame(x = R$xcol[R$valid_col], y = R$yrow[R$valid_row], value = values)
  ev <- R$events[[as.character(year)]]
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile(width = diff(R$xcol[1:2]), height = diff(R$yrow[1:2])) +
    scale_fill_gradientn(colours = PAL, na.value = "white", name = expression(lambda(x))) +
    coord_fixed() + theme_minimal() +
    annotate("rect", xmin = min(R$xcol), xmax = max(R$xcol), ymin = min(R$yrow), ymax = max(R$yrow),
             fill = NA, colour = "black", linewidth = 0.3) +
    scale_x_continuous(breaks = xb, labels = round(lon(xb), 2), name = "Longitude") +
    scale_y_continuous(breaks = yb, labels = round(lat(yb), 2), name = "Latitude") +
    geom_point(data = ev, aes(x = x, y = y), inherit.aes = FALSE, shape = 16, colour = "turquoise1", size = 0.9) +
    ggtitle(as.character(year))
}
plot_spatial <- function(region, model, years, name, ncol = 3) {
  pool <- get_pool(region, model); if (is.null(pool)) return(invisible(NULL))
  R <- get_region(region)
  plots <- lapply(years, function(y) map_plot(R, plugin_intensity(pool, R, y), y))
  nrow <- ceiling(length(plots) / ncol)
  save_plot(wrap_plots(plots, ncol = ncol), name, 4.6 * ncol, 4.2 * nrow)
}

# ------------------------------------------------------------------ 2D (temperature, precipitation)
plot_2d <- function(region = "ontario") {
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
  p2 <- ggplot(df, aes(x = x, y = y, fill = value)) + geom_raster() +
    scale_fill_gradientn(colours = PAL, na.value = "white", name = expression(rho(z))) + axes
  plots <- list(p1, p2)
  if (use_kernel) {
    df$kern <- ifelse(inside, kernel_2d(region, "temp", "prec", z), NA)
    plots[[3]] <- ggplot(df, aes(x = x, y = y, fill = kern)) + geom_raster() +
      scale_fill_gradientn(colours = PAL, na.value = "white", name = expression(rho(z))) + axes
  }
  save_plot(wrap_plots(plots, nrow = 1), sprintf("%s_2d_estim.pdf", FILE_PREFIX[[region]]), 4.4 * length(plots), 3.8)
}

# ------------------------------------------------------------------ 3D: (temp, prec) slices at wind quantiles
# Posterior mean on the whole (temp, prec) square at fixed wind quantiles.
plot_slices <- function(region = "ontario", probs = c(0.05, 0.5, 0.95)) {
  pool <- get_pool(region, "full"); if (is.null(pool)) return(invisible(NULL))
  R <- get_region(region); g <- seq(0, 1, length.out = 200)
  z12 <- as.matrix(expand.grid(z1 = g, z2 = g))
  u3 <- quantile(as.vector(R$U$wind), probs, names = FALSE)
  w_orig <- canada_backtransform(R, "wind", u3)
  vals <- lapply(u3, function(u) posterior_rho(pool, cbind(z12, u))$mean)
  lim <- range(unlist(vals))
  plots <- lapply(seq_along(probs), function(k)
    ggplot(data.frame(x = z12[, 1], y = z12[, 2], value = vals[[k]]), aes(x = x, y = y, fill = value)) +
      geom_raster() +
      scale_fill_gradientn(colours = PAL, limits = lim, name = expression(rho(z))) +
      coord_fixed() + theme_minimal() + labs(x = CANADA_LABELS[["temp"]], y = CANADA_LABELS[["prec"]],
      subtitle = sprintf("Wind speed: %s quantile (%.2f)", format(probs[k]), w_orig[k])) +
      scale_x_continuous(breaks = axis_breaks, labels = orig_labels(R, "temp")) +
      scale_y_continuous(breaks = axis_breaks, labels = orig_labels(R, "prec")))
  save_plot(wrap_plots(plots, nrow = 1, guides = "collect"), "2d_slicing.pdf", 13, 4)
  cat("Wind quantiles", paste(probs, collapse = "/"), "in km/h:", paste(round(w_orig, 2), collapse = ", "), "\n")
}

# ------------------------------------------------------------------ figures of the paper
run_fig(plot_1d("ontario"))                                                               # Fig. 5
run_fig(plot_spatial("ontario", "temp", c(2013, 2015, 2021), "ontario_temp.pdf"))          # Fig. 6
run_fig(plot_slices("ontario"))                                                           # Fig. 7
run_fig(plot_spatial("ontario", "full", c(2013, 2015, 2021), "ontario_3d.pdf"))            # Fig. 8
run_fig(plot_2d("ontario"))                                                               # Suppl. 2D
run_fig(plot_spatial("ontario", "temp_prec", c(2013, 2015, 2021), "ontario_2d_spat.pdf"))  # Suppl. 2D maps
run_fig(plot_spatial("ontario", "full", c(2006, 2008, 2010, 2011, 2018, 2019), "ontario_3d_add.pdf"))
run_fig(plot_spatial("ontario", "full", 2004:2009, "ontario_3d_1.pdf"))
run_fig(plot_spatial("ontario", "full", c(2010, 2011, 2012, 2014, 2016, 2017), "ontario_3d_2.pdf"))
run_fig(plot_spatial("ontario", "full", c(2018, 2019, 2020, 2022), "ontario_3d_3.pdf", ncol = 2))
nine <- c(2006, 2008, 2010, 2011, 2013, 2015, 2018, 2019, 2021)
for (region in c("saskatchewan", "british_columbia")) {
  run_fig(plot_1d(region))
  run_fig(plot_spatial(region, "full", nine, paste0(SPATIAL_3D_FILE[[region]], ".pdf")))
}
