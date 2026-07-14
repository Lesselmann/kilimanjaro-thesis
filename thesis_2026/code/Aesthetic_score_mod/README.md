# code/Aesthetic_score_mod/

Independent branch: scores images on perceived **aesthetic quality**, separate from the
indoor/outdoor/human/landscape classification cascade. All three notebooks were run on
outdoor/landscape image subsets loaded from Google Drive.

| Notebook | Method | What it does |
|---|---|---|
| `Aesthetic_Predictor_V2_5.ipynb` | CLIP ViT-L/14 + MLP head (trained on SAC+LOGOS+AVA1) | Scores 1,009 photos on a ~4–10 aesthetic scale. Writes `aesthetic_v2_scores.csv`, distribution plot, top/bottom-5 visualizations. |
| `CLIP_IQA+.ipynb` | CLIP-IQA+ (learned prompts) | Scores 897 outdoor photos on a 0–1 quality scale (mean = 0.567). Writes `clipiqa_scores.csv`. |
| `NIMA_class.ipynb` | NIMA (MobileNet, AVA-trained) | Scores 897 outdoor photos on a 1–10 scale with uncertainty (mean = 5.27). Writes `nima_scores.csv`. |

These scores are not currently tied to a `data/` output folder in this repo (outputs were saved to
Google Drive during the Colab runs) — copy the resulting CSVs into `data/` if they should be
version-controlled alongside the rest of the pipeline.
