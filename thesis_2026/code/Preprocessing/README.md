# code/Preprocessing/

Prepares the raw Flickr dataset before any classification happens.

| Notebook | What it does |
|---|---|
| `Flickr_Photo_Map.ipynb` | Visualizes the spatial distribution of all photos around Kilimanjaro on a 1 km×1 km UTM grid, using `geopandas`/`folium`/`contextily` with an OSM basemap (log-scaled photo density). Output feeds `figures/maps/kilimanjaro_photo_map.png`. |
| `Duplicate_Detection.ipynb` | Detects near-duplicate images via perceptual hashing (pHash, threshold=10, `imagehash` library) on a manually sampled subset; writes a duplicate-pairs list and a cleaned image list. |
| `Random_selection_Outdoor_Images.ipynb` | Draws a stratified random sample of 300 images (seed=44) from the outdoor subset for manual annotation; writes the sample + blank annotation-template CSV/Excel. |

**Inputs:** `data/final_flickr_dataset_with_metadata.csv`, `data/All_images/`
**Outputs:** `data/Duplicates/duplicates_all_images.csv`, random-sample folders under `data/Annotation_*/`
