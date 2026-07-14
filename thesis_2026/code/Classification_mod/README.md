# code/Classification_mod/

Single model/prompt classification runs on the gold-standard samples (~300 images each) — these
are exploratory experiments, distinct from the full-dataset runs in `Final_classification/`.

| Notebook | Task | Model | Prompt style | Result |
|---|---|---|---|---|
| `01072026_class_individual_coeco_landscape_4models_annotation_predictions.ipynb.ipynb` | Individual / CoEco / Landscape | All 4 models compared | mixed | 297 images; best: MobileCLIP-S1 + mixed prompts (Macro-F1 = 0.787) |
| `28062026_class_human_nonhuman_eva01clip_detailed_prompts.ipynb.ipynb` | Human / Non-human | EVA01-CLIP-g/14 | detailed | 300 images → 140 human / 160 non-human |
| `28062026_class_human_nonhuman_vitl14clip_detailed_prompts.ipynb.ipynb` | Human / Non-human | CLIP ViT-L/14 | detailed | 300 images → 124 human / 176 non-human |
| `28062026_class_indoor_outdoor_mobileclips1_simple_prompts.ipynb.ipynb` | Indoor / Outdoor | MobileCLIP-S1 | simple | 300 images → 14 indoor / 286 outdoor |
| `28062026_linear_probe_human_nonhuman_training.ipynb` | Human / Non-human | MobileCLIP embeddings + Logistic Regression | — | Trains the first linear probe on 298 gold-standard embeddings (5-fold CV AUC = 0.981); writes `code/Linear_probe/` artifacts |

Each classification notebook writes a prediction CSV (`file_name`, `pred_label`, `confidence`) and
copies the source images into per-class subfolders (e.g. `human/`, `non-human/`) under the
matching `data/Annotation_*/` folder.
