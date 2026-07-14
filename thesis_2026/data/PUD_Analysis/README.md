# data/PUD_Analysis/

Photo-User-Days (PUD) recreation-intensity analysis on the 3,422 `landscape/`-classified photos,
produced by `code/PUD/01072026_pud_calculation_landscape_1km_grid.ipynb.ipynb`. 228 grid cells
have ≥1 PUD; max PUD = 61; mean annual PUD = 0.220.

| File | Description |
|---|---|
| `landscape_puds_1km_grid_polygons.csv` | Per-cell statistics: `cell_id`, `cell_x`/`cell_y`, `PUD_total`, `n_images`, `n_unique_users`, `lat_mean`/`lon_mean`, `year_min`/`year_max`, `PUD_annual_mean`, `PUD_per_user_year`. |
| `landscape_puds_1km_grid_polygons_utm.gpkg` | Same grid as GeoPackage polygons, UTM 37S projection (used for the area-based PUD calculation). |
| `landscape_puds_1km_grid_polygons_wgs84.gpkg` | Same grid reprojected to WGS84 (for web-map display). |
| `landscape_puds_map_final.png`, `landscape_puds_map_osm.png` | PUD-intensity choropleth maps (with/without OSM basemap). |
| `landscape_puds_per_user_map_osm.png` | Map normalized by unique users per cell. |
| `landscape_unique_users_map_final.png` | Map of unique-user counts per cell. |
