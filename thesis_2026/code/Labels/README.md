# code/Labels/

| Notebook | What it does |
|---|---|
| `creating_photo_labels_17042026.ipynb` | Scans an `indoor/`/`outdoor/` folder structure and generates a label CSV from the folder names; skips augmented copies (`augmented_*`). Produced `photo_labels_indoor_outdoor_20042026.csv` — 1,889 unique photos (450 indoor, 1,439 outdoor), 3,111 augmented duplicates excluded. |

This is the ground-truth label generation step used before the indoor/outdoor gold-standard
annotation and evaluation work in `code/Annotation/`.
