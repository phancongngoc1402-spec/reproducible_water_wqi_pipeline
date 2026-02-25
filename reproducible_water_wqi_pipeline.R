#!/usr/bin/env Rscript
# ============================================================
# Reproducible pipeline for Water Quality Index (WQI),
# clustering (K-means), PCA, summary figures (boxplot & exceedance),
# and kriging maps.
#
# Outputs are written under --out_dir:
#   tables/  : normalized data, WQI tables, clustering table
#   figures/ : Fig2 boxplot, Fig4 exceedance, elbow, density, PCA, kriging maps
#   logs/    : sessionInfo.txt
#
# How to run (example):
#   Rscript reproducible_water_wqi_pipeline.R \
#     --data_water "data/data_water.xlsx" \
#     --coord_file "data/Water_coordinate_points.xlsx" \
#     --wqi_raw    "data/data_WQI.xlsx" \
#     --border_shp "data/border-nghean2.shp" \
#     --out_dir    "outputs" \
#     --k_opt      3 \
#     --grid_cellsize 250 \
#     --epsg_projected 9208 \
#     --seed 123
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c(
    "optparse", "readxl", "writexl", "dplyr", "janitor",
    "tidyr", "stringr", "forcats", "scales",
    "ggplot2", "viridis",
    "sf", "sp", "gstat", "purrr"
  )
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
  lapply(pkgs, library, character.only = TRUE)
})

# -----------------------
# 0) CLI arguments
# -----------------------
option_list <- list(
  optparse::make_option("--data_water", type = "character", default = "data/data_water.xlsx",
                        help = "Excel file containing raw water parameters (default: %default)"),
  optparse::make_option("--coord_file", type = "character", default = "data/Water_coordinate_points.xlsx",
                        help = "Excel file containing station coordinates (station, lon, lat) (default: %default)"),
  optparse::make_option("--wqi_raw", type = "character", default = "data/data_WQI.xlsx",
                        help = "Excel file containing WQI points for kriging (station, round, WQI100, lon, lat) (default: %default)"),
  optparse::make_option("--border_shp", type = "character", default = "data/border-nghean2.shp",
                        help = "Shapefile boundary of study area (default: %default)"),
  optparse::make_option("--out_dir", type = "character", default = "outputs",
                        help = "Output directory (default: %default)"),
  optparse::make_option("--k_opt", type = "integer", default = 3,
                        help = "Number of clusters for K-means (default: %default)"),
  optparse::make_option("--grid_cellsize", type = "double", default = 250,
                        help = "Grid cell size (meters) for kriging (default: %default)"),
  optparse::make_option("--epsg_projected", type = "integer", default = 9208,
                        help = "Projected CRS for kriging/maps (default: %default)"),
  optparse::make_option("--seed", type = "integer", default = 123,
                        help = "Random seed for reproducibility (default: %default)"),
  optparse::make_option("--use_round_labels", type = "logical", default = TRUE,
                        help = "Use manuscript-style labels for 4 rounds (default: %default)")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Output folders
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$out_dir, "logs"), showWarnings = FALSE, recursive = TRUE)

# Log session info
sink(file.path(opt$out_dir, "logs", "sessionInfo.txt"))
print(sessionInfo())
sink()

set.seed(opt$seed)

# -----------------------
# 1) Read and clean raw data
# -----------------------
message("1) Reading raw dataset: ", opt$data_water)
raw <- readxl::read_excel(opt$data_water)
dat <- raw %>% janitor::clean_names()

# -----------------------
# 2) Standards and variables
# -----------------------
std <- list(
  nitrite   = 0.05,
  ammonium  = 0.3,
  chloride  = 250,
  fluoride  = 1,
  arsenic   = 0.01,
  lead      = 0.02,
  manganese = 0.1,
  iron      = 0.5,
  e_coli    = 20,
  bod       = 4,
  cod       = 10,
  tss       = 25,
  coliform  = 1000
)
pollutant_vars <- names(std)

# -----------------------
# 3) Normalize indices I_x
# -----------------------
message("2) Normalizing indices I_x ...")
dat_norm <- dat %>%
  dplyr::mutate(
    dplyr::across(
      .cols  = dplyr::all_of(pollutant_vars),
      .fns   = ~ .x / std[[dplyr::cur_column()]],
      .names = "I_{.col}"
    )
  ) %>%
  dplyr::mutate(
    I_ph = dplyr::case_when(
      .data$ph < 6.5 ~ 6.5 / .data$ph,
      .data$ph > 8.5 ~ .data$ph / 8.5,
      TRUE           ~ 1
    ),
    I_do = dplyr::if_else(.data$do > 0, 6 / .data$do, NA_real_)
  )

# I_max
index_cols <- grep("^I_", names(dat_norm), value = TRUE)
dat_norm <- dat_norm %>%
  dplyr::rowwise() %>%
  dplyr::mutate(I_max = max(dplyr::c_across(dplyr::all_of(index_cols)), na.rm = TRUE)) %>%
  dplyr::ungroup()

writexl::write_xlsx(dat_norm, file.path(opt$out_dir, "tables", "data_water_normalized.xlsx"))

# -----------------------
# 4) Compute WQI (0–100)
# -----------------------
message("3) Computing WQI ...")
I_vars <- grep("^I_", names(dat_norm), value = TRUE)
I_vars <- setdiff(I_vars, "I_max")

w <- rep(1 / length(I_vars), length(I_vars))
names(w) <- I_vars

dat_wqi <- dat_norm %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    exceed_sum = sum(w * pmax(dplyr::c_across(dplyr::all_of(I_vars)) - 1, 0)),
    WQI   = 1 / (1 + exceed_sum),
    WQI100 = WQI * 100
  ) %>%
  dplyr::ungroup()

writexl::write_xlsx(dat_wqi, file.path(opt$out_dir, "tables", "data_water_wqi.xlsx"))

# -----------------------
# 5) K-means clustering using WQI100 + I_vars
# -----------------------
message("4) K-means clustering ...")
cluster_vars <- c("WQI100", I_vars)

X <- dat_wqi %>%
  dplyr::select(dplyr::all_of(cluster_vars)) %>%
  as.matrix()

X_scaled <- scale(X)

# Elbow (1..10)
wss <- numeric(10)
for (k in 1:10) {
  km <- stats::kmeans(X_scaled, centers = k, nstart = 20)
  wss[k] <- km$tot.withinss
}
elbow_df <- data.frame(k = 1:10, wss = wss)

p_elbow <- ggplot2::ggplot(elbow_df, ggplot2::aes(x = k, y = wss)) +
  ggplot2::geom_line(linewidth = 0.9, color = "black") +
  ggplot2::geom_point(size = 2.5, color = "black") +
  ggplot2::geom_vline(xintercept = opt$k_opt, linetype = "dashed", linewidth = 0.8, color = "black") +
  ggplot2::annotate(
    "text", x = opt$k_opt + 0.4,
    y = elbow_df$wss[elbow_df$k == opt$k_opt],
    label = paste("Selected k =", opt$k_opt),
    hjust = 0, size = 4, color = "black"
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::labs(x = "Number of clusters (k)", y = "Total within-cluster sum of squares", title = "")

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", paste0("elbow_kmeans_k", opt$k_opt, ".tiff")),
  plot = p_elbow,
  device = "tiff", dpi = 600,
  width = 16, height = 12, units = "cm",
  compression = "lzw"
)

km_model <- stats::kmeans(X_scaled, centers = opt$k_opt, nstart = 50)
dat_kmeans <- dat_wqi %>% dplyr::mutate(cluster = factor(km_model$cluster))

p_density <- ggplot2::ggplot(dat_kmeans, ggplot2::aes(x = WQI100, color = cluster)) +
  ggplot2::geom_density(linewidth = 1) +
  ggplot2::labs(x = "Water Quality Index (WQI)", y = "Density") +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(legend.position = "right")

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", paste0("cluster_density_k", opt$k_opt, ".tiff")),
  plot = p_density,
  device = "tiff", dpi = 600,
  width = 18, height = 16, units = "cm",
  compression = "lzw"
)

writexl::write_xlsx(dat_kmeans, file.path(opt$out_dir, "tables", "data_water_wqi_kmeans.xlsx"))

# -----------------------
# 6) PCA on I_vars
# -----------------------
message("5) PCA analysis ...")
pca_data <- dat_wqi %>%
  dplyr::select(dplyr::all_of(I_vars)) %>%
  tidyr::drop_na()

pca_res <- stats::prcomp(pca_data, scale. = TRUE)
loadings <- pca_res$rotation

loading_df <- as.data.frame(loadings) %>%
  dplyr::mutate(parameter = rownames(loadings)) %>%
  dplyr::select(parameter, PC1, PC2)

loading_long <- loading_df %>%
  tidyr::pivot_longer(cols = c("PC1", "PC2"), names_to = "PC", values_to = "loading") %>%
  dplyr::mutate(key = stringr::str_remove(parameter, "^I_"))

label_map <- c(
  cod       = "COD",
  bod       = "BOD",
  coliform  = "Total coliform",
  e_coli    = "E. coli",
  ammonium  = "NH\u2084\u207A\u2013N",
  nitrite   = "NO\u2082\u207B\u2013N",
  tss       = "TSS",
  do        = "DO",
  ph        = "pH",
  iron      = "Fe",
  manganese = "Mn",
  chloride  = "Cl\u207B",
  fluoride  = "F\u207B",
  arsenic   = "As",
  lead      = "Pb"
)

loading_long <- loading_long %>%
  dplyr::mutate(param_label = dplyr::if_else(key %in% names(label_map), label_map[key], key))

order_tbl <- loading_long %>%
  dplyr::group_by(param_label) %>%
  dplyr::summarise(max_abs = max(abs(loading), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(max_abs)

loading_long <- loading_long %>%
  dplyr::left_join(order_tbl, by = "param_label") %>%
  dplyr::mutate(param_label = factor(param_label, levels = order_tbl$param_label))

p_load_bar <- ggplot2::ggplot(loading_long, ggplot2::aes(x = param_label, y = loading, fill = PC)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.85), width = 0.78) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  ggplot2::coord_flip() +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::labs(x = "Water quality parameters", y = "Loading value", fill = "PC")

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", "pca_loadings_PC1_PC2.tiff"),
  plot = p_load_bar,
  device = "tiff", dpi = 600,
  width = 18, height = 14, units = "cm",
  compression = "lzw"
)

scores <- as.data.frame(pca_res$x)
dat_pca <- dat_wqi %>%
  dplyr::filter(!is.na(WQI100)) %>%
  dplyr::bind_cols(scores)

p_pca_scores <- ggplot2::ggplot(dat_pca, ggplot2::aes(x = PC1, y = PC2, color = WQI100)) +
  ggplot2::geom_point(alpha = 0.7) +
  viridis::scale_color_viridis(option = "D") +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::labs(color = "WQI", title = "")

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", "pca_scores_colored_by_wqi.tiff"),
  plot = p_pca_scores,
  device = "tiff", dpi = 600,
  width = 20, height = 14, units = "cm",
  compression = "lzw"
)

# -----------------------
# 7) Join coordinates to WQI table and save (this file is used for Fig.2 & Fig.4)
# -----------------------
message("6) Joining station coordinates ...")
coord <- readxl::read_excel(opt$coord_file) %>%
  janitor::clean_names()

coord <- coord %>%
  dplyr::mutate(station = toupper(trimws(as.character(station))))

dat_wqi2 <- dat_wqi %>%
  dplyr::mutate(station = toupper(trimws(as.character(station)))) %>%
  dplyr::left_join(coord, by = "station")

out_dat_wqi2 <- file.path(opt$out_dir, "tables", "dat_wqi_with_coordinates.xlsx")
writexl::write_xlsx(dat_wqi2, out_dat_wqi2)

missing_coords <- dat_wqi2 %>% dplyr::filter(is.na(lon) | is.na(lat))
if (nrow(missing_coords) > 0) {
  message("WARNING: Some stations have missing coordinates: ",
          paste(unique(missing_coords$station), collapse = ", "))
}

# -----------------------
# 8) Fig.2 Boxplot of WQI100 by monitoring round (monochrome EMA style)
# -----------------------
message("7) Creating Fig.2 (boxplot WQI100 by round) ...")

stopifnot("round" %in% names(dat_wqi2))
if (!("WQI100" %in% names(dat_wqi2))) stop("Column 'WQI100' not found after WQI computation.")

# If WQI100 accidentally stored as 0-1, rescale
if (max(dat_wqi2$WQI100, na.rm = TRUE) <= 1.2) {
  dat_wqi2 <- dat_wqi2 %>% dplyr::mutate(WQI100 = 100 * WQI100)
}

# Optional manuscript-style round labels (expects 4 campaigns coded as 1..4)
if (isTRUE(opt$use_round_labels)) {
  round_labels <- c(
    "Round 1 (Mar 2025)",
    "Round 2 (Jun 2025)",
    "Round 3 (Sep 2025)",
    "Round 4 (Dec 2025)"
  )
  # If round is numeric 1..4 => label; otherwise keep as-is
  if (all(sort(unique(dat_wqi2$round)) %in% 1:4)) {
    dat_wqi2 <- dat_wqi2 %>%
      dplyr::mutate(round = factor(round, levels = 1:4, labels = round_labels))
  } else {
    dat_wqi2 <- dat_wqi2 %>% dplyr::mutate(round = as.factor(round))
  }
} else {
  dat_wqi2 <- dat_wqi2 %>% dplyr::mutate(round = as.factor(round))
}

p_fig2 <- ggplot2::ggplot(dat_wqi2, ggplot2::aes(x = round, y = WQI100)) +
  ggplot2::geom_boxplot(
    width = 0.62,
    fill = "white",
    color = "black",
    outlier.shape = 1,
    outlier.size = 2
  ) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    shape = 17,  # triangle
    size = 3,
    color = "black"
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 100)) +
  ggplot2::labs(x = NULL, y = "WQI100") +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(size = 10),
    panel.grid.major.y = ggplot2::element_line(color = "gray90")
  )

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", "Fig2_boxplot_WQI100_monochrome.tiff"),
  plot     = p_fig2,
  device   = "tiff",
  dpi      = 600,
  width    = 20,
  height   = 14,
  units    = "cm",
  compression = "lzw"
)

# -----------------------
# 9) Fig.4 Exceedance frequency (%) by parameter (Ii > 1) (monochrome)
# -----------------------
message("8) Creating Fig.4 (exceedance frequency) ...")

i_cols <- names(dat_wqi2)[stringr::str_detect(names(dat_wqi2), "^I_")]
i_cols <- setdiff(i_cols, "I_max")
if (length(i_cols) == 0) stop("No columns starting with 'I_' found in dat_wqi2. Check normalization step.")

freq <- dat_wqi2 %>%
  dplyr::select(dplyr::all_of(i_cols)) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "param", values_to = "Ii") %>%
  dplyr::mutate(
    param = stringr::str_remove(param, "^I_"),
    exceed = .data$Ii > 1
  ) %>%
  dplyr::group_by(param) %>%
  dplyr::summarise(exceed_pct = mean(exceed, na.rm = TRUE) * 100, .groups = "drop") %>%
  dplyr::mutate(
    param_label = dplyr::if_else(param %in% names(label_map), label_map[param], param),
    param_label = forcats::fct_reorder(as.character(param_label), exceed_pct)
  )

p_fig4 <- ggplot2::ggplot(freq, ggplot2::aes(x = param_label, y = exceed_pct)) +
  ggplot2::geom_col(fill = "white", color = "black", width = 0.78) +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = ggplot2::expansion(mult = c(0, 0.06))
  ) +
  ggplot2::labs(x = NULL, y = "Exceedance frequency (%)") +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_line(color = "grey90")
  )

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", "Fig4_exceedance_frequency_monochrome.tiff"),
  plot     = p_fig4,
  device   = "tiff",
  dpi      = 600,
  width    = 20,
  height   = 14,
  units    = "cm",
  compression = "lzw"
)

# -----------------------
# 10) Kriging maps (by round)
# -----------------------
message("9) Kriging mapping ...")
wqi_raw <- readxl::read_excel(opt$wqi_raw) %>% janitor::clean_names()
border <- sf::st_read(opt$border_shp, quiet = TRUE)

# If lon/lat in wqi_raw are degrees -> EPSG:4326; if already projected meters -> set to opt$epsg_projected.
crs_wqi_input <- 4326

wqi_sf <- sf::st_as_sf(wqi_raw, coords = c("lon", "lat"), crs = crs_wqi_input)
wqi_sf <- sf::st_transform(wqi_sf, opt$epsg_projected)
border <- sf::st_transform(border, opt$epsg_projected)

bbox <- sf::st_bbox(border)

# Grid
cellsize <- opt$grid_cellsize
xs <- seq(bbox["xmin"], bbox["xmax"], by = cellsize)
ys <- seq(bbox["ymin"], bbox["ymax"], by = cellsize)
grid_df0 <- expand.grid(x = xs, y = ys)
grid_sf0 <- sf::st_as_sf(grid_df0, coords = c("x", "y"), crs = opt$epsg_projected)

inside <- sf::st_intersects(grid_sf0, border, sparse = FALSE)
grid_sf <- grid_sf0[inside, , drop = FALSE]
grid_sp <- as(grid_sf, "Spatial")

krige_one_round <- function(r, points_sf, grid_sp) {
  
  message("  - Processing round: ", r)
  pts_r <- points_sf %>% dplyr::filter(round == r)
  
  if (nrow(pts_r) < 3) {
    message("    Skipping: not enough points for variogram/kriging.")
    return(NULL)
  }
  
  pts_sp <- as(pts_r, "Spatial")
  v_emp <- gstat::variogram(WQI100 ~ 1, pts_sp)
  
  vgm_init <- gstat::vgm(
    psill  = stats::var(pts_sp$WQI100, na.rm = TRUE) * 0.5,
    model  = "Sph",
    range  = max(sp::spDists(pts_sp)) / 3,
    nugget = stats::var(pts_sp$WQI100, na.rm = TRUE) * 0.1
  )
  
  v_fit <- gstat::fit.variogram(v_emp, vgm_init, fit.method = 2)
  kr <- gstat::krige(WQI100 ~ 1, pts_sp, grid_sp, model = v_fit)
  
  coords <- sp::coordinates(kr)
  kr_df  <- as.data.frame(kr)
  kr_df$x <- coords[, 1]
  kr_df$y <- coords[, 2]
  kr_df$round <- r
  
  names(kr_df)[names(kr_df) == "var1.pred"] <- "WQI_pred"
  kr_df <- kr_df[, c("x", "y", "WQI_pred", "round")]
  kr_df
}

rounds <- sort(unique(wqi_sf$round))
pred_list <- purrr::map(rounds, ~ krige_one_round(.x, points_sf = wqi_sf, grid_sp = grid_sp))
pred_all_df <- dplyr::bind_rows(pred_list)

pts_df <- wqi_sf %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[, 1],
    y = sf::st_coordinates(.)[, 2]
  ) %>%
  sf::st_drop_geometry()

p_cont <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = pred_all_df, ggplot2::aes(x = x, y = y, fill = WQI_pred), alpha = 0.9) +
  ggplot2::geom_sf(data = border, fill = NA, color = "black", linewidth = 0.4) +
  ggplot2::geom_point(data = pts_df, ggplot2::aes(x = x, y = y), size = 1, color = "black") +
  viridis::scale_fill_viridis(name = "WQI", option = "C", direction = -1) +
  ggplot2::coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE,
    crs = sf::st_crs(opt$epsg_projected)
  ) +
  ggplot2::facet_wrap(~ round, ncol = 2) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::labs(title = "Kriging WQI surface (projected CRS)", x = "X (m)", y = "Y (m)")

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", "kriging_continuous.tiff"),
  plot = p_cont,
  device = "tiff", dpi = 600,
  width = 18, height = 16, units = "cm",
  compression = "lzw"
)

lvl <- c("Very Poor", "Poor", "Moderate", "Good", "Excellent")
cols_wqi <- c(
  "Very Poor" = "#440154",
  "Poor"      = "#3b528b",
  "Moderate"  = "#21918c",
  "Good"      = "#5ec962",
  "Excellent" = "#fde725"
)

pred_class <- pred_all_df %>%
  dplyr::mutate(
    WQI_class = cut(
      WQI_pred,
      breaks = c(-Inf, 25, 50, 75, 90, Inf),
      labels = lvl,
      include.lowest = TRUE
    ),
    WQI_class = factor(WQI_class, levels = lvl)
  )

dummy_legend <- data.frame(
  x = as.numeric(bbox["xmin"]) + 1,
  y = as.numeric(bbox["ymin"]) + 1,
  WQI_class = factor("Very Poor", levels = lvl),
  round = rounds[1]
)

p_class <- ggplot2::ggplot() +
  ggplot2::geom_point(
    data = dummy_legend,
    ggplot2::aes(x = x, y = y, fill = WQI_class),
    shape = 22, size = 6, alpha = 0, show.legend = TRUE
  ) +
  ggplot2::geom_raster(data = pred_class, ggplot2::aes(x = x, y = y, fill = WQI_class), alpha = 0.9) +
  ggplot2::geom_sf(data = border, fill = NA, color = "black", linewidth = 0.4) +
  ggplot2::geom_point(data = pts_df, ggplot2::aes(x = x, y = y), size = 1, color = "black") +
  ggplot2::scale_fill_manual(
    name = "WQI level",
    values = cols_wqi,
    breaks = lvl,
    limits = lvl,
    drop = FALSE
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
  ggplot2::coord_sf(
    xlim = c(as.numeric(bbox["xmin"]), as.numeric(bbox["xmax"])),
    ylim = c(as.numeric(bbox["ymin"]), as.numeric(bbox["ymax"])),
    expand = FALSE,
    crs = sf::st_crs(opt$epsg_projected)
  ) +
  ggplot2::facet_wrap(~ round, ncol = 2) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::labs(title = "", x = "X (m)", y = "Y (m)")

ggplot2::ggsave(
  filename = file.path(opt$out_dir, "figures", "kriging_classified.tiff"),
  plot = p_class,
  device = "tiff",
  dpi = 600,
  width = 16,
  height = 16,
  units = "cm",
  compression = "lzw"
)

message("DONE. Outputs saved to: ", opt$out_dir)