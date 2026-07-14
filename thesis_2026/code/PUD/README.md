# code/PUD/

**PUD = Photo-User-Days**, a standard recreation-ecology metric for landscape-use intensity
(number of distinct user-days a location was photographed).

| Notebook | What it does |
|---|---|
| `01072026_pud_calculation_landscape_1km_grid.ipynb.ipynb` | Takes the 3,422 "landscape"-classified photos (from `Final_classification/`), joins them with the original Flickr metadata (owner, upload date), filters to 2004–2024 (3,349 images), projects coordinates to UTM 37S, assigns each photo to a 1 km grid cell, and counts unique (user, date) pairs per cell. Result: 228 grid cells with ≥1 PUD, max PUD = 61, mean annual PUD = 0.220. |

**Inputs:** `data/Classification_Individual_CoEco_Landscape/.../landscape/`,
`data/final_flickr_dataset_with_metadata.csv`
**Outputs:** `data/PUD_Analysis/` — PUD statistics CSV, GeoPackage grid polygons (UTM + WGS84), and
map visualizations.
