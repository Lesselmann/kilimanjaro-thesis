# Master's Thesis 2026 — Landscape Photo Classification & Recreation-Use Analysis

This repository documents the full pipeline of a master's thesis that classifies geotagged Flickr
photos from the Kilimanjaro region (Tanzania) using vision-language models (CLIP variants,
MobileCLIP, EVA-CLIP) and a ResNet-50 baseline, then derives a landscape recreation-intensity
metric (Photo-User-Days, PUD) from the resulting "landscape" photo subset.

This README is the entry point of a **documentation cascade**: every major folder (and most
subfolders) contains its own `README.md` explaining exactly what is inside. A machine-readable
catalog of every folder/file, with description, creation date and file counts, is kept in
[`mother_table.csv`](mother_table.csv) at the repository root — open it in Excel/Sheets for a
sortable/filterable overview.

---

## How to navigate this repo

| Folder | What it is | Local README |
|---|---|---|
| [`code/`](code/README.md) | All Jupyter notebooks, grouped by pipeline stage | [code/README.md](code/README.md) |
| [`data/`](data/README.md) | All raw data, intermediate results, images and final outputs | [data/README.md](data/README.md) |
| [`figures/`](figures/README.md) | Standalone figures used outside the `data/` pipeline | [figures/README.md](figures/README.md) |
| [`literature/`](literature/README.md) | Bibliography (.bib) | [literature/README.md](literature/README.md) |
| [`texts/`](texts/README.md) | LaTeX thesis manuscript & literature-review document | [texts/README.md](texts/README.md) |
| [`mother_table.csv`](mother_table.csv) | Master catalog of every folder/file in this repo | — |

---

## The classification pipeline, in one picture

The core task is a **cascading 3-level classification** of every photo in the deduplicated
dataset, each level narrowing down the previous level's output. A parallel, independent branch
scores aesthetic quality; a final branch computes recreation-use intensity (PUD) on the resulting
landscape photos.

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
   ══════ Recreation-use analysis ══════
        code/PUD/  →  data/PUD_Analysis/
        Photo-User-Days per 1 km grid cell (landscape photos only)
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

---

## Reproducing the pipeline

All classification/evaluation notebooks were developed on **Google Colab** and contain
`drive.mount()` cells and Google-Drive-relative paths — adjust these to local paths before running
outside Colab. See [code/README.md](code/README.md) for the notebook-by-notebook run order.

## Requirements

Core packages used across notebooks:

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
