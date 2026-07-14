# code/Evalu_Class_mod/

Systematic evaluation across **all 4 models × all 3 prompt styles**, with proper dev/test splits,
plus the full linear-probe benchmark. This is where the "best model + method per task" decisions
(used later in `Final_classification/`) were made.

| Notebook | Task | What it does |
|---|---|---|
| `20062026_human_nonhuman_4models_prompttuning_evaluation.ipynb.ipynb` | Human / Non-human | 3 prompt variants (simple/detailed/mixed) × 4 models, stratified 35% dev / 65% test split. Per-model prompt-comparison CSV + test-set confusion matrix. Best: CLIP ViT-L/14 + mixed prompts. |
| `24062026_indoor_outdoor_4models_prompttuning_evaluation.ipynb.ipynb` | Indoor / Outdoor | Same dev/test framework for indoor/outdoor. |
| `30062026_individual_coeco_landscape_4models_prompttuning_evaluation.ipynb.ipynb` | Individual / CoEco / Landscape | Same dev/test framework (35% dev split). |
| `29062026_Evaluation_CLIP_with_linear_probe.ipynb` | Human / Non-human | Trains linear probes for all 4 models on the 298 gold-standard samples (5-fold CV, no dev/test split). Best: CLIP ViT-L/14 (AUC = 0.9914, Accuracy = 94.3%). Writes `data/LinearProbe_AllModels_Comparison/`. |
| `01072026_linear_probe_allmodels_individual_coeco_landscape.ipynb.ipynb` | Individual / CoEco / Landscape | Same linear-probe benchmark for all 4 models (297 gold-standard images). Ranking (Macro-F1): MobileCLIP-S1 (0.819) > ViT-B/32 (0.811) > EVA01 (0.795) > ViT-L/14 (0.788). Writes `data/LinearProbe_AllModels_Comparison_IndividualCoEcoLandscape/`. |
| `Eva_Resnet50.ipynb` | Indoor / Outdoor | Evaluates the fine-tuned ResNet-50 (fold 2) on 700 held-out images (223 indoor, 477 outdoor). Accuracy = 74.57% — used as the supervised-baseline comparison point against the CLIP-based zero-shot/probe approaches. |

**Inputs:** gold-standard annotation CSVs from `data/Annotation_*/`.
**Outputs:** per-model metrics/prompt-comparison CSVs and linear-probe artifacts, feeding into
`data/Model_Evaluation_*/` and `data/LinearProbe_AllModels_Comparison*/`.
