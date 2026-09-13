# Master's Thesis 2026 — Mapping Landscape Aesthetics from Social Media Images

*Which landscape variables are significant predictors of perceived landscape aesthetics in the
Mount Kilimanjaro region, modelled from georeferenced Flickr images and deep-learning methods?*

This repository documents the full pipeline of a master's thesis that classifies geotagged Flickr
photos from the Kilimanjaro region (Tanzania) using vision-language models (CLIP variants,
MobileCLIP, EVA-CLIP) and a ResNet-50 baseline, then uses Photo-User-Days (PUD) over the resulting
"landscape" photo subset as a **proxy for perceived landscape aesthetic value** and models it
against a set of landscape / scenery indicators.

This README is the entry point of a **documentation cascade**: every major folder (and most
subfolders) contains its own `README.md` explaining exactly what is inside. A machine-readable
catalog of every folder/file, with description, creation date and file counts, is kept in
[`mother_table.csv`](mother_table.csv) at the repository root — open it in Excel/Sheets for a
sortable/filterable overview.

---

## How to navigate this repo

| Folder | What it is | Local README |
|---|---|---|
| [`code/`](code/README.md) | Jupyter notebooks (classification, PUD) + R scripts (landscape indicators, modelling), grouped by pipeline stage | [code/README.md](code/README.md) |
| [`data/`](data/README.md) | All raw data, intermediate results, images and final outputs | [data/README.md](data/README.md) |
| [`figures/`](figures/README.md) | Standalone figures used outside the `data/` pipeline | [figures/README.md](figures/README.md) |
| [`literature/`](literature/README.md) | Bibliography (.bib) | [literature/README.md](literature/README.md) |
| [`texts/`](texts/README.md) | LaTeX thesis manuscript & literature-review document | [texts/README.md](texts/README.md) |
| [`mother_table.csv`](mother_table.csv) | Master catalog of every folder/file in this repo | — |

---

## The classification pipeline, in one picture

The core task is a **cascading 3-level classification** of every photo in the deduplicated
dataset, each level narrowing down the previous level's output. A parallel, independent branch
scores aesthetic quality; a final branch computes Photo-User-Days (PUD) on the resulting landscape
photos as a **proxy for perceived landscape aesthetic value** and **models it against ~18 landscape /
scenery indicators** per 1 km grid cell.

```
data/final_flickr_dataset_with_metadata.csv   (12,332 Flickr photos, metadata + GPS)
                │
                ▼
        data/All_images/                       (raw JPGs downloaded)
                │  Preprocessing: duplicate detection (perceptual hashing)
                ▼
        data/All_images_deduplicated/           (8,819 images — working dataset)
                │
                ▼
   ══════ LEVEL 1 — Indoor vs. Outdoor ══════
        data/Classification_Indoor_Outdoor/
                │  best model: MobileCLIP-S1, simple prompts
                ▼
        …/final_class_indoor_outdoor_sorted/{indoor, outdoor}
                │  (outdoor branch continues)
                ▼
   ══════ LEVEL 2 — Human vs. Non-human ══════
        data/Classification_Human_non_Human/
                │  4 prediction variants: MobileCLIP / ViT-L14 × zero-shot / linear-probe
                ▼
        …/{…}_sorted/{human, non_human}
                │  (non_human branch continues)
                ▼
   ══════ LEVEL 3 — Individual vs. Community/Ecosystem vs. Landscape ══════
        data/Classification_Individual_CoEco_Landscape/
                │  2 prediction variants: MobileCLIP mixed-prompt / MobileCLIP linear-probe
                ▼
        …/final_class_individual_coeco_landscape_sorted/{individual, community_ecosystem, landscape}
                │  (landscape branch continues)
                ▼
   ══════ Landscape-aesthetics proxy (PUD) ══════
        code/generate_PUD/         →  data/PUD_Analysis/  (Photo-User-Days per 1 km cell)
        code/generate_number_images/ →  data/PUD_Analysis/  (raw image count per cell)
                │
                ▼
   ══════ Landscape-indicator modelling (R) ══════
        code/Landscape_Indicators/  →  data/Landscape_indicators/  (~18 predictors per 1 km cell:
                                       viewshed, landscape metrics, Sentinel indices, distances,
                                       structural & colour diversity)   [git-ignored, local only]
                │  code/C_Screening/   (Spearman + VIF collinearity screening)
                │  code/Variable_data_distrubution/  (distribution diagnostics)
                ▼
        code/Modeling_PUD/  →  data/Model_training_testing/
        LM / GLM / GAM / Random Forest / XGBoost  for  avg_annual_PUD  and  n_images
        (backward selection, spatial + random CV)
```

Each classification level was first **validated on a gold-standard subset** before being applied
to the full corpus:

```
Random sample of images  →  manual annotation (2–3 raters)  →  inter-annotator agreement
        →  zero-shot model evaluation vs. gold standard (simple/detailed/mixed prompts)
        →  linear-probe training on model embeddings (all 4 models)
        →  dev/test split evaluation, best model+prompt selected per task
        →  applied at full-dataset scale ("Final_classification")
```

This validation work lives in `data/Annotation_Indoor_Outdoor/`, `data/Annotation_Human_Non_Human/`,
`data/Annotation_Individual_CoEco_Landscape/`, and `data/Model_Evaluation_*_28062026/` /
`_30062026/`. See [data/README.md](data/README.md) for the full breakdown.

### Models used

| Model | Type | Role |
|---|---|---|
| CLIP ViT-B/32 (`openai/clip-vit-base-patch32`) | Vision-language, zero-shot | Compared across all 3 tasks |
| CLIP ViT-L/14 (`openai/clip-vit-large-patch14`) | Vision-language, zero-shot | Compared across all 3 tasks; used for some final classifications |
| EVA01-CLIP-g/14 | Vision-language, zero-shot | Compared across all 3 tasks |
| MobileCLIP-S1 | Vision-language, zero-shot, lightweight | Best-performing model on most tasks; used for final large-scale classification |
| Linear probe (Logistic Regression on frozen embeddings) | Supervised, trained per model | Trained per model per task on the gold-standard sample; often outperforms pure zero-shot |
| ResNet-50 (fine-tuned) | CNN, supervised baseline | Indoor/outdoor only, for comparison against CLIP-based zero-shot approaches |

### Aesthetic scoring (independent branch)

`code/Aesthetic_score_mod/` scores outdoor/landscape images on perceived aesthetic quality using
three independent methods (NIMA, CLIP-IQA+, Aesthetic Predictor V2.5). This is not part of the
indoor/outdoor/human/landscape cascade — it is a separate quality dimension computed on the same
images.

### Landscape indicators & aesthetics modelling (R)

Once the `landscape` photo subset exists, perceived landscape aesthetic value is proxied per 1 km
grid cell by **Photo-User-Days** (`code/generate_PUD/`) and raw **image count**
(`code/generate_number_images/`).
`code/Landscape_Indicators/` then computes ~18 scenery / landscape predictors per cell
(viewshed openness, landscape-metric diversity, visible land-cover composition, Sentinel-2
indices, distances to summit / water / trails, vertical structural heterogeneity, colour
diversity, …) from DEMs, Sentinel-2 rasters and vector layers held in `data/Landscape_indicators/`
(**git-ignored**, several GB). After collinearity screening (`code/C_Screening/`) and
distribution checks (`code/Variable_data_distrubution/`), `code/Modeling_PUD/` fits and compares
LM, GLM, GAM, Random Forest and XGBoost for both response variables, with backward selection and
spatial + random cross-validation. Results: `data/Model_training_testing/`.

---

## Reproducing the pipeline

The classification/evaluation **notebooks** were developed on **Google Colab** and contain
`drive.mount()` cells and Google-Drive-relative paths — adjust to local paths before running
outside Colab. The landscape-indicator and modelling **R scripts** were developed locally and
use absolute Windows paths (`C:\Users\Lukas\…`) — adjust before running. See
[code/README.md](code/README.md) for the stage-by-stage run order.

## Requirements

### Python (classification, aesthetic scoring, PUD grids)

```
torch torchvision
transformers open_clip_torch
tensorflow keras          # NIMA
pandas numpy scikit-learn joblib
Pillow imagehash tqdm
matplotlib seaborn
geopandas folium contextily shapely   # maps, PUD grid
```

```bash
pip install torch torchvision transformers open_clip_torch tensorflow \
            pandas numpy scikit-learn joblib pillow imagehash tqdm matplotlib seaborn \
            geopandas folium contextily shapely
```

### R (landscape indicators & aesthetics modelling)

```
sf terra raster exactextractr           # geodata, zonal stats
landscapemetrics                        # SHDI, contagion, shape/edge metrics
mgcv                                    # GAM
MASS                                    # GLM (Negative Binomial)
randomForest xgboost                    # ML models
blockCV                                 # spatial cross-validation
car                                     # VIF
moments                                 # skewness
dplyr
```

Rasters are produced by a Google Earth Engine script
(`code/Landscape_Indicators/files/files/22_2026-07-17_GEE_Datadownload_FINAL_All7Exports.js`)
run in the Earth Engine Code Editor.

---

## Note on repository documentation

The **repository documentation** (the per-folder `README.md` files and
[`mother_table.csv`](mother_table.csv)), together with the `.gitignore`, was drafted in August
2026 with the help of an AI assistant (Anthropic Claude, via Claude Code) and then reviewed by
the author.
