# code/Linear_probe/

Saved artifacts from the **first** linear-probe experiment (Human/Non-human task, MobileCLIP only),
produced by `code/Classification_mod/28062026_linear_probe_human_nonhuman_training.ipynb`. This is
not a notebook folder — it only holds the trained artifacts.

| File | Description |
|---|---|
| `embeddings_298_mobileclip.npy` | MobileCLIP-S1 image embeddings for the 298 gold-standard human/non-human samples |
| `labels_298_gold.npy` | Corresponding gold-standard labels (human / non_human) |
| `linear_probe_mobileclip_human_nonhuman.joblib` | Trained scikit-learn `LogisticRegression` probe (5-fold CV AUC = 0.981); applied to new MobileCLIP embeddings to classify unseen images |

The later, more complete linear-probe comparison across **all 4 models** and **both remaining
tasks** lives in `data/LinearProbe_AllModels_Comparison/` and
`data/LinearProbe_AllModels_Comparison_IndividualCoEcoLandscape/`, produced by notebooks in
`code/Evalu_Class_mod/`.
