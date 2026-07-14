# code/ — Notebooks by Pipeline Stage

All code is written as Jupyter notebooks, developed on Google Colab (they contain `drive.mount()`
cells — adjust Drive paths to local paths to run outside Colab). See the [root README](../README.md)
for the overall pipeline picture.

## Suggested run order

| Stage | Folder | Purpose | README |
|---|---|---|---|
| 1 | [`Preprocessing/`](Preprocessing/README.md) | Map visualization, duplicate detection, random sampling | [Preprocessing/README.md](Preprocessing/README.md) |
| 2 | [`Labels/`](Labels/README.md) | Build label CSVs from folder structure | [Labels/README.md](Labels/README.md) |
| 3 | [`Annotation/`](Annotation/README.md) | Inter-annotator agreement & gold-standard vs. model evaluation | [Annotation/README.md](Annotation/README.md) |
| 4 | [`Classification_mod/`](Classification_mod/README.md) | Single model/prompt classification experiments + linear-probe training | [Classification_mod/README.md](Classification_mod/README.md) |
| 5 | [`Linear_probe/`](Linear_probe/README.md) | Saved linear-probe artifacts (embeddings, labels, trained model) | [Linear_probe/README.md](Linear_probe/README.md) |
| 6 | [`Evalu_Class_mod/`](Evalu_Class_mod/README.md) | Multi-model / multi-prompt evaluation, dev/test splits, linear-probe benchmarks | [Evalu_Class_mod/README.md](Evalu_Class_mod/README.md) |
| 7 | [`Resnet50_Training/`](Resnet50_Training/README.md) | Fine-tune ResNet-50 supervised baseline (indoor/outdoor) | [Resnet50_Training/README.md](Resnet50_Training/README.md) |
| 8 | [`Final_classification/`](Final_classification/README.md) | Apply the best model/method per task to the full dataset (thousands of images) | [Final_classification/README.md](Final_classification/README.md) |
| 9 | [`PUD/`](PUD/README.md) | Photo-User-Days recreation-intensity calculation on a 1 km grid | [PUD/README.md](PUD/README.md) |
| — | [`Aesthetic_score_mod/`](Aesthetic_score_mod/README.md) | Independent aesthetic-quality scoring (not part of the classification cascade) | [Aesthetic_score_mod/README.md](Aesthetic_score_mod/README.md) |

## Naming convention

Most notebook filenames start with a `DDMMYYYY_` date prefix indicating when that experiment was
run, followed by the task (`indoor_outdoor`, `human_nonhuman`, `individual_coeco_landscape`), the
model, and the prompt style (`simple_prompts` / `detailed_prompts` / `mixed_prompts`).
