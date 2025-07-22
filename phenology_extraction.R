# GIMMS NDVI3g+ phenology extraction using phenofit
# Compatible with NDVI 24-point half-monthly time series

library(terra)
library(phenofit)
library(imputeTS)
library(data.table)

# Define reference extent (e.g., Northern Hemisphere)
ref_rast <- rast(res = 1/12, xmin = -180, xmax = 180, ymin = 0, ymax = 90, crs = "EPSG:4326")

# Load and preprocess NDVI, QC, and FT data for one year
load_inputs <- function(year, vi_dir, qc_dir, ft_dir, replace_raster_path) {
  ndvi <- list.files(file.path(vi_dir, year), pattern = "\\.tif$", full.names = TRUE) %>%
    rast() %>% crop(ref_rast)

  qc <- list.files(file.path(qc_dir, year), pattern = "\\.tif$", full.names = TRUE) %>%
    rast() %>% crop(ref_rast)

  ft <- list.files(file.path(ft_dir, year), pattern = "\\.tif$", full.names = TRUE) %>%
    rast() %>% crop(ref_rast)

  winter_replace <- rast(replace_raster_path) %>% crop(ref_rast)

  # Apply masks: annual NDVI mean > 0.1 and max occurs between Mar–Oct
  mean_mask <- classify(app(ndvi, mean, na.rm = TRUE), matrix(c(-Inf, 0.1, NA), ncol = 3))
  max_pos_mask <- classify(app(ndvi, which.max, na.rm = TRUE), matrix(c(-Inf, 4.5, NA, 20.5, Inf, NA), ncol = 3, byrow = TRUE))
  ndvi <- mask(ndvi, mean_mask) %>% mask(max_pos_mask)
  qc <- mask(qc, mean_mask) %>% mask(max_pos_mask)
  ft <- mask(ft, mean_mask) %>% mask(max_pos_mask)
  winter_replace <- mask(winter_replace, mean_mask) %>% mask(max_pos_mask)

  return(c(ndvi, qc, ft, winter_replace))
}

# Main phenology extraction function (per pixel)
extract_phenology <- function(x) {
  nptperyear <- 24
  GIMMS_DOY <- c(8,24,39,53,67,83,98,113,128,144,159,174,189,205,220,236,251,266,281,297,312,327,342,358)
  dates <- as.Date(GIMMS_DOY, origin = "1982-12-31")

  x <- as.numeric(x)
  x[!is.finite(x)] <- NA
  ndvi_qc <- x[1:(2*nptperyear)]
  ft_code <- x[(2*nptperyear + 1):(3*nptperyear)]
  winter_replace <- x[3*nptperyear + 1]

  df <- data.frame(date = dates, y = ndvi_qc[seq(1, length(ndvi_qc), 2)], QC = ndvi_qc[seq(2, length(ndvi_qc), 2)])
  if (sum(complete.cases(df)) < 0.7 * nptperyear) return(rep(NA, 16))

  df$y <- na_interpolation(df$y)
  df$QC[is.na(df$QC)] <- 1
  df <- cbind(df, qc_NDVI3g(df$QC))
  df$t <- df$date
  df$y[ft_code == 0] <- pmax(df$y[ft_code == 0], winter_replace, na.rm = TRUE)

  d <- as.data.table(df[, c("date", "t", "y", "QC_flag", "w")])
  INPUT <- check_input(d$t, d$y, d$w, QC_flag = d$QC_flag, nptperyear = nptperyear, south = FALSE,
                       maxgap = 6, wmin = 0.2, wsnow = 0.8, ymin = 0.08, alpha = 0.02)

  fit_methods <- list("smooth_wWHIT", "smooth_wHANTS", "smooth_wSG")
  R2_list <- list()
  brks_list <- list()

  for (m in fit_methods) {
    res <- tryCatch(season_mov(INPUT, list(rFUN = m, wFUN = "wTSM", iters = 2, len_min = 60, maxExtendMonth = 12)), error = function(e) NULL)
    R2_list[[m]] <- if (!is.null(res)) res$GOF["R2"] else NA
    brks_list[[m]] <- if (!is.null(res)) res else NULL
  }

  best_fit <- which.max(unlist(R2_list))
  if (all(is.na(unlist(R2_list)))) return(rep(NA, 16))
  brks <- brks_list[[best_fit]]

  fit <- curvefits(INPUT, brks, options = list(methods = "Elmore", wFUN = "wTSM", iters = 3, wmin = 0.2, maxExtendMonth = 6))
  re_pheno <- tryCatch(get_pheno(fit, TRS = c(0.1, 0.15, 0.2, 0.5), IsPlot = FALSE), error = function(e) NULL)
  if (is.null(re_pheno)) return(rep(NA, 16))

  phenos <- re_pheno$doy$Elmore[, c("TRS1.sos", "TRS1.eos", "TRS2.sos", "TRS2.eos",
                                    "DER.sos", "DER.eos", "Greenup", "Senescence",
                                    "Maturity", "Dormancy")]
  return(unlist(phenos))
}
