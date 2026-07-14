# data/ — Overview

~84,000 files. See the [root README](../README.md) for how these folders chain together into the
3-level classification cascade. Each subfolder below has its own `README.md` with details; for
folders full of thousands of category images, only a folder-level description is given (see
[mother_table.csv](../mother_table.csv) at the repo root for per-folder file counts).

## Root-level files

| File | Description |
|---|---|
| `final_flickr_dataset_with_metadata.csv` | The full raw dataset — 12,332 geotagged Flickr photos with metadata (photo ID, owner, GPS coordinates, upload/taken date, view/fave/comment counts, tags, direct image URL). Source of truth for all downstream joins (e.g. PUD analysis). |
| `README_FinalDATASET.xlsx` | Documentation of the dataset fields and the Flickr collection/query process. |
| `Model_Eva_Results_indoor_outdoor.xlsx` | Aggregated side-by-side comparison of all models' indoor/outdoor evaluation results. |

## Subfolders — pipeline order

| # | Folder | Stage | README |
|---|---|---|---|
| 1 | [`All_images/`](All_images/README.md) | Raw dataset | [README](All_images/README.md) |
| 2 | [`All_images_deduplicated/`](All_images_deduplicated/README.md) | Deduplicated working dataset | [README](All_images_deduplicated/README.md) |
| 3 | [`Duplicates/`](Duplicates/README.md) | Duplicate-detection output | [README](Duplicates/README.md) |
| 4 | [`manual_deleted_images_outdoor_26062026/`](manual_deleted_images_outdoor_26062026/README.md) | Manually removed outdoor images | [README](manual_deleted_images_outdoor_26062026/README.md) |
| 5 | [`Annotation_Indoor_Outdoor/`](Annotation_Indoor_Outdoor/README.md) | Level 1 gold standard: sampling, annotation, agreement, evaluation | [README](Annotation_Indoor_Outdoor/README.md) |
| 6 | [`Annotation_Human_Non_Human/`](Annotation_Human_Non_Human/README.md) | Level 2 gold standard | [README](Annotation_Human_Non_Human/README.md) |
| 7 | [`Annotation_Individual_CoEco_Landscape/`](Annotation_Individual_CoEco_Landscape/README.md) | Level 3 gold standard | [README](Annotation_Individual_CoEco_Landscape/README.md) |
| 8 | [`Model_Evaluation_Indoor_Outdoor_28062026/`](Model_Evaluation_Indoor_Outdoor_28062026/README.md) | Level 1 dev/test model benchmark | [README](Model_Evaluation_Indoor_Outdoor_28062026/README.md) |
| 9 | [`Model_Evaluation_Human_Non_Human_28062026/`](Model_Evaluation_Human_Non_Human_28062026/README.md) | Level 2 dev/test model benchmark | [README](Model_Evaluation_Human_Non_Human_28062026/README.md) |
| 10 | [`Model_Evaluation_Individual_CoEco_Landscape_30062026/`](Model_Evaluation_Individual_CoEco_Landscape_30062026/README.md) | Level 3 dev/test model benchmark | [README](Model_Evaluation_Individual_CoEco_Landscape_30062026/README.md) |
| 11 | [`LinearProbe_AllModels_Comparison/`](LinearProbe_AllModels_Comparison/README.md) | Linear-probe artifacts, Level 2, all 4 models | [README](LinearProbe_AllModels_Comparison/README.md) |
| 12 | [`LinearProbe_AllModels_Comparison_IndividualCoEcoLandscape/`](LinearProbe_AllModels_Comparison_IndividualCoEcoLandscape/README.md) | Linear-probe artifacts, Level 3, all 4 models | [README](LinearProbe_AllModels_Comparison_IndividualCoEcoLandscape/README.md) |
| 13 | [`Classification_Indoor_Outdoor/`](Classification_Indoor_Outdoor/README.md) | **Final** Level 1 classification, full dataset | [README](Classification_Indoor_Outdoor/README.md) |
| 14 | [`Classification_Human_non_Human/`](Classification_Human_non_Human/README.md) | **Final** Level 2 classification, full dataset | [README](Classification_Human_non_Human/README.md) |
| 15 | [`Classification_Individual_CoEco_Landscape/`](Classification_Individual_CoEco_Landscape/README.md) | **Final** Level 3 classification, full dataset | [README](Classification_Individual_CoEco_Landscape/README.md) |
| 16 | [`PUD_Analysis/`](PUD_Analysis/README.md) | Photo-User-Days recreation-intensity grid (landscape photos) | [README](PUD_Analysis/README.md) |
| 17 | [`weights/`](weights/README.md) | Downloaded model weights | [README](weights/README.md) |

## Common file patterns you'll see repeated across folders

- **`*_sorted/` folders** contain the *same images*, physically copied into per-class
  subfolders (e.g. `human/`, `non-human/`, `indoor/`, `outdoor/`, `individual/`,
  `community_ecosystem/`, `landscape/`) — this is the classification *output*, not new data.
- **`file_name, file_path, pred_label, confidence`** is the standard prediction-CSV schema used
  by (almost) every classification notebook.
- **Model folder names** (`CLIP-ViT-B32`, `CLIP-ViT-L14`, `EVA01-CLIP-g14`, `MobileCLIP-S1`)
  recur throughout `Model_Evaluation_*` and `LinearProbe_*` folders — same 4 models, compared per
  task.
