###############################################################################
# canada_data.R
# Loading and pre-processing of the Canadian wildfire data (June, 2004-2022) for the
# revised analysis. Only base R and Matrix are needed: the spatstat objects stored in the
# .RData files are read as plain lists (im: $v, $xcol, $yrow, $xstep, $ystep; ppp: $x, $y).
#
# Pre-processing (as in the original CODE_Canada_analysis.R, plus the border-event rule):
#   * years selected by name (2004-2022, n = 19);
#   * each covariate mapped to [0, 1] by the empirical c.d.f. of its pixel values pooled over
#     the 19 years and all valid pixels (identical to the original `rescaling`);
#   * the likelihood integral uses every valid pixel of every year as a quadrature node, with
#     weight equal to the pixel area (xstep * ystep, i.e. the original area/pixels);
#   * events take the covariate value of the pixel containing them (nearest pixel centre, as
#     spatstat's `[.im`); events falling on an NA pixel (coast/border) are assigned to the
#     nearest valid pixel centre instead of being dropped.
###############################################################################

CANADA_DATA_DIR <- "Data"                  # Data Analysis/Application/Data (run from Application/)
CANADA_YEARS <- 2004:2022

# region name -> object prefix in the .RData files
CANADA_REGIONS <- c(ontario = "ontario", saskatchewan = "saska", british_columbia = "bt_columbia")
# covariate name -> object infix
CANADA_COVARIATES <- c(temp = "temp", prec = "precip", wind = "wind")
CANADA_LABELS <- c(temp = "Temperature", prec = "Precipitation", wind = "Wind speed")
CANADA_MODELS <- list(temp = "temp", prec = "prec", wind = "wind",
                      temp_prec = c("temp", "prec"), full = c("temp", "prec", "wind"))

# The 13 fits of the paper (main text + Supplement, Section 4)
CANADA_PAPER_FITS <- data.frame(
  region = c(rep("ontario", 5), rep("saskatchewan", 4), rep("british_columbia", 4)),
  model = c("temp", "prec", "wind", "temp_prec", "full",
            "temp", "prec", "wind", "full", "temp", "prec", "wind", "full"),
  stringsAsFactors = FALSE)
CANADA_PAPER_FITS$fit_id <- paste(CANADA_PAPER_FITS$region, CANADA_PAPER_FITS$model, sep = "__")

# Output folder of a fit inside <out>/Fits: capitalised region + model, e.g. ontario__temp ->
# Ontario_temp, british_columbia__temp_prec -> British_Columbia_temp_prec. The internal fit id
# (region__model: options, seeds, tables) is unchanged.
fit_folder <- function(fit_id) vapply(strsplit(fit_id, "__", fixed = TRUE), function(p) {
  w <- strsplit(p[1], "_", fixed = TRUE)[[1]]
  paste0(paste0(toupper(substring(w, 1, 1)), substring(w, 2), collapse = "_"), "_", p[2])
}, character(1))

# Pixel -> geographical coordinates for the axis labels (copied from the original script).
CANADA_GEO <- list(
  ontario = list(lon = c(-95.16, -74.34), lon_px = c(378.9401, 629.2343),
                 lat = c(41.68, 56.85), lat_px = c(29.52635, 263.85878), n_xbreaks = 5),
  saskatchewan = list(lon = c(-110, -101.3), lon_px = c(212.75, 319.87),
                      lat = c(49, 60), lat_px = c(132.68, 330.95), n_xbreaks = 3),
  british_columbia = list(lon = c(-139.03, -114.03), lon_px = c(32.538, 172.7),
                          lat = c(48.18, 60), lat_px = c(162.72, 438.09), n_xbreaks = 3))

canada_data_file <- function(region, data_dir = CANADA_DATA_DIR) {
  if (!region %in% names(CANADA_REGIONS)) stop("Unknown region: ", region)
  file.path(data_dir, paste0(CANADA_REGIONS[[region]], "_data_june_orig_scale.RData"))
}

# Returns a list with the pixel geometry, transformed and raw covariates at the valid pixels
# (one column per year), events (coordinates, pixel index, snapping record), and the pooled
# sorted raw values used for back-transformation to the original scale.
load_canada_region <- function(region, years = CANADA_YEARS, data_dir = CANADA_DATA_DIR) {
  file <- canada_data_file(region, data_dir)
  if (!file.exists(file)) stop("Data file not found: ", normalizePath(file, mustWork = FALSE))
  e <- new.env(parent = emptyenv()); load(file, envir = e)
  prefix <- CANADA_REGIONS[[region]]; yrs <- as.character(years)
  get_years <- function(infix) {
    L <- get(paste0(prefix, "_", infix, "_june"), envir = e)
    missing <- setdiff(yrs, names(L))
    if (length(missing)) stop(region, ": years missing for ", infix, ": ", paste(missing, collapse = ","))
    L[yrs]
  }
  ims <- lapply(CANADA_COVARIATES, get_years)
  pps <- get_years("ppp")
  ref <- ims$temp[[1]]
  mask <- !is.na(ref$v)
  for (cv in names(ims)) for (y in yrs) {
    im <- ims[[cv]][[y]]
    if (!identical(dim(im$v), dim(ref$v)) || !isTRUE(all.equal(im$xcol, ref$xcol)) ||
        !isTRUE(all.equal(im$yrow, ref$yrow)) || !identical(!is.na(im$v), mask))
      stop(region, ": pixel grid or mask differs for ", cv, " ", y)
  }
  valid <- which(mask)                       # column-major linear indices, as im$v[mask]
  nr <- nrow(mask)
  vr <- (valid - 1L) %% nr + 1L; vc <- (valid - 1L) %/% nr + 1L
  raw <- lapply(ims, function(L) vapply(L, function(im) im$v[valid], numeric(length(valid))))
  pooled <- lapply(raw, function(M) sort(as.vector(M)))
  U <- lapply(names(raw), function(cv) {
    Fz <- ecdf(pooled[[cv]])
    matrix(Fz(raw[[cv]]), nrow(raw[[cv]]), ncol(raw[[cv]]), dimnames = list(NULL, yrs))
  })
  names(U) <- names(raw)

  events <- lapply(yrs, function(y) {
    p <- pps[[y]]; x <- as.numeric(p$x); yy <- as.numeric(p$y)
    col <- round((x - ref$xcol[1]) / ref$xstep) + 1
    row <- round((yy - ref$yrow[1]) / ref$ystep) + 1
    col <- pmin(pmax(col, 1), ncol(mask)); row <- pmin(pmax(row, 1), nr)   # as spatstat: edge pixel
    lin <- (col - 1) * nr + row
    pix <- match(lin, valid)
    snapped <- is.na(pix); dist <- rep(0, length(x))
    for (i in which(snapped)) {
      d2 <- (ref$xcol[vc] - x[i])^2 + (ref$yrow[vr] - yy[i])^2
      pix[i] <- which.min(d2); dist[i] <- sqrt(min(d2))
    }
    data.frame(x = x, y = yy, pix = pix, snapped = snapped, snap_distance = dist)
  })
  names(events) <- yrs
  list(region = region, years = years, file = normalizePath(file), md5 = unname(tools::md5sum(file)),
       mask = mask, valid = valid, valid_row = vr, valid_col = vc,
       xcol = ref$xcol, yrow = ref$yrow, pixel_area = ref$xstep * ref$ystep,
       raw = raw, U = U, pooled = pooled, events = events,
       n_events = vapply(events, nrow, integer(1)),
       n_snapped = vapply(events, function(ev) sum(ev$snapped), integer(1)))
}

# Data in the format of the sampler (Zq, qw, Zev per year) for a set of covariates
canada_reps <- function(R, covs) {
  lapply(seq_along(R$years), function(i) {
    ev <- R$events[[i]]
    list(Zq = do.call(cbind, lapply(covs, function(cv) R$U[[cv]][, i])),
         qw = R$pixel_area,
         Zev = do.call(cbind, lapply(covs, function(cv) R$U[[cv]][ev$pix, i])))
  })
}

# Back-transform from [0, 1] to the original covariate scale (as in the original script)
canada_backtransform <- function(R, cv, u)
  unname(quantile(R$pooled[[cv]], probs = u, type = 1, names = FALSE))

canada_data_summary <- function(R) {
  data.frame(region = R$region, years = paste(range(R$years), collapse = "-"),
             n = length(R$years), valid_pixels = length(R$valid),
             pixel_area = R$pixel_area, region_area = R$pixel_area * length(R$valid),
             events = sum(R$n_events), min_events = min(R$n_events), max_events = max(R$n_events),
             snapped_events = sum(R$n_snapped),
             max_snap_distance = max(0, unlist(lapply(R$events, `[[`, "snap_distance"))),
             data_md5 = R$md5)
}

# Monitored evaluation points (probes) used for convergence diagnostics
canada_probes <- function(d) {
  if (d == 1) return(matrix(seq(0, 1, length.out = 21), ncol = 1))
  g <- if (d == 2) c(0.1, 0.3, 0.5, 0.7, 0.9) else c(0.2, 0.5, 0.8)
  as.matrix(expand.grid(rep(list(g), d)))[, d:1, drop = FALSE]   # z1 slowest
}
