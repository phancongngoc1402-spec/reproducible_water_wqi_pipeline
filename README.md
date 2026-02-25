# reproducible_water_wqi_pipeline
# Reproducible Water Quality Index (WQI) Workflow (R)

This repository provides a fully reproducible R workflow to:
1) normalize measured water-quality parameters into dimensionless indices (I_x) relative to regulatory standards,  
2) compute a Water Quality Index (WQI) on a 0–100 scale,  
3) perform K-means clustering and PCA to explore patterns in water quality, and  
4) produce spatial WQI maps using ordinary kriging for each monitoring round.

The workflow is implemented as a **single R script** designed to be executed end-to-end and to generate all tables/figures in a standardized output folder.

---

## 1. Reproducibility Principles

- **Single-entry script**: the complete pipeline is in one file (`reproducible_water_wqi_pipeline.R`).
- **Command-line parameters**: file paths and key options can be passed via CLI arguments, ensuring portability across machines.
- **Deterministic outputs**: a fixed random seed is used for K-means and elbow computations (`--seed`, default `123`).
- **Session logging**: `sessionInfo()` is saved to `outputs/logs/sessionInfo.txt` for transparency and reproducibility.
- **No hard-coded absolute paths**: the script uses relative paths by default.

---

## 2. Requirements

### 2.1 R Version
- Recommended: **R ≥ 4.2.0**  
  (The script may work on older versions, but this is the recommended baseline for consistent package behavior.)

### 2.2 R Packages
The script will automatically install missing packages from CRAN, then load them:
- `optparse`, `readxl`, `writexl`, `dplyr`, `janitor`, `tidyr`, `stringr`
- `ggplot2`, `viridis`
- `sf`, `sp`, `gstat`
- `purrr`

> Note: Internet access is needed **only** the first time for package installation.

---

## 3. Repository Structure
project/
reproducible_water_wqi_pipeline.R
README.md
data/
data_water.xlsx
Water_coordinate_points.xlsx
data_WQI.xlsx
border-nghean2.shp
border-nghean2.dbf
border-nghean2.shx
border-nghean2.prj
outputs/ # generated automatically
tables/
figures/
logs/

---

## 4. Input Data Specifications

### 4.1 `data/data_water.xlsx` (raw water parameters)
Required columns (case-insensitive; script cleans names via `janitor::clean_names()`):
- `station` (station ID; string)
- `round` (sampling round; numeric or string)
- Water-quality parameters (examples):
  - `nitrite`, `ammonium`, `chloride`, `fluoride`, `arsenic`, `lead`,
  - `manganese`, `iron`, `e_coli`, `ph`, `do`, `bod`, `cod`, `tss`, `coliform`

**Units** must match the standards used in the script (see Section 5).

### 4.2 `data/Water_coordinate_points.xlsx` (station coordinates)
Required columns:
- `station` (station ID)
- `lon` and `lat`

> The script standardizes station IDs with `toupper(trimws())` before joining.

### 4.3 `data/data_WQI.xlsx` (WQI points for kriging)
Required columns:
- `station`
- `round`
- `WQI100` (WQI on 0–100 scale)
- `lon`, `lat`

**Coordinate assumption**: by default, the script assumes `lon/lat` are geographic coordinates in **EPSG:4326 (WGS84)** and then transforms them to the projected CRS used for kriging (default EPSG:9208).  
If your `lon/lat` are already projected (meters) in EPSG:9208, edit the script variable `crs_wqi_input` accordingly.

### 4.4 `data/border-nghean2.shp` (study area boundary)
A polygon shapefile of the study region. All shapefile sidecar files (`.shx`, `.dbf`, `.prj`, etc.) must be present.

---

## 5. Method Summary

### 5.1 Normalization to Dimensionless Indices (I_x)
For parameters with an upper permissible limit \( C_{std} \), an index is computed as:

\[
I_x = \frac{C_{measured}}{C_{std}}
\]

- If \( I_x > 1 \), the parameter exceeds the standard.
- A summary indicator \( I_{max} \) is computed per sample as the maximum of all \( I \) values.

**Special handling**:
- **pH**: normalized to the allowable range 6.5–8.5  
  - below 6.5: \( I_{pH} = 6.5/pH \)  
  - above 8.5: \( I_{pH} = pH/8.5 \)  
  - otherwise: \( I_{pH} = 1 \)
- **DO**: treated as a lower-bound standard (DO ≥ 6 mg/L):  
  \[
  I_{DO} = \frac{6}{DO}
  \]
  (with a safeguard against division by zero)

### 5.2 Water Quality Index (WQI)
This workflow uses an exceedance-based, equal-weight aggregation approach:

1) Define equal weights \( w_i = 1/n \) for \( n \) indices.  
2) Compute exceedance above 1:
\[
E = \sum_{i=1}^{n} w_i \cdot \max(I_i - 1, 0)
\]
3) Convert to bounded WQI in (0,1]:
\[
WQI = \frac{1}{1+E}
\]
4) Convert to 0–100 scale:
\[
WQI_{100} = 100 \cdot WQI
\]

> **Important**: This formulation is a transparent, reproducible exceedance-based WQI.  
> If your manuscript follows a different published WQI definition (e.g., NSF-WQI, weighted arithmetic WQI, CCME-WQI, etc.), you should update the formula section and/or script accordingly.

### 5.3 K-means Clustering
Clustering is performed on `WQI100` and all indices `I_*`:
- Standardization: z-score scaling
- Elbow method: within-cluster sum of squares for k=1..10
- Final k: user-defined via `--k_opt` (default 3)

### 5.4 PCA
PCA is computed on index variables `I_*` (excluding `I_max`) with scaling.
- Outputs include bar charts of loadings for PC1 and PC2, and a scores plot colored by WQI.

### 5.5 Spatial Interpolation (Ordinary Kriging)
For each sampling round:
- empirical variogram is computed and fitted (default initial model: spherical),
- ordinary kriging predicts WQI over a grid masked to the study boundary,
- both continuous and classified maps are produced.

---

## 6. How to Run

### 6.1 Run with default paths
From the project root:

```bash
Rscript reproducible_water_wqi_pipeline.R
Rscript reproducible_water_wqi_pipeline.R \
  --data_water "data/data_water.xlsx" \
  --coord_file "data/Water_coordinate_points.xlsx" \
  --wqi_raw    "data/data_WQI.xlsx" \
  --border_shp "data/border-nghean2.shp" \
  --out_dir    "outputs" \
  --k_opt      3 \
  --grid_cellsize 250 \
  --epsg_projected 9208 \
  --seed 123
7. Outputs

All outputs are written to outputs/ (or the folder specified by --out_dir).

7.1 Tables (outputs/tables/)

data_water_normalized.xlsx
Raw data + normalized indices I_* + I_max

data_water_wqi.xlsx
Adds WQI (0–1) and WQI100 (0–100)

data_water_wqi_kmeans.xlsx
Adds cluster assignment

data_wqi_with_coordinates.xlsx
WQI100 joined with station coordinates

7.2 Figures (outputs/figures/)

elbow_kmeans_k<k>.tiff
Elbow curve for choosing k

cluster_density_k<k>.tiff
Density plot of WQI by cluster

pca_loadings_PC1_PC2.tiff
PCA loadings bar chart (PC1/PC2)

pca_scores_colored_by_wqi.tiff
PCA scores colored by WQI

kriging_continuous.tiff
Continuous kriged WQI surfaces by round

kriging_classified.tiff
Classified WQI surfaces by round

7.3 Logs (outputs/logs/)

sessionInfo.txt
R version + package versions
