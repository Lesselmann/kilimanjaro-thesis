# code/Annotation/

Builds and validates the **gold standard** for the Indoor/Outdoor task: human annotation agreement,
then zero-shot model performance against that gold standard.

| Notebook | What it does |
|---|---|
| `28062026_interannotator_agreement_indoor_outdoor_2raters.ipynb` | Measures agreement between 2 human annotators (NB, LE) on 300 indoor/outdoor images. Cohen's κ = 0.813, raw agreement = 98%. Writes disagreement list for reconciliation. |
| `28062026_interannotator_agreement_indoor_outdoor_3raters.ipynb` | Extends to a 3rd rater (MED); computes pairwise Cohen's κ and Fleiss' κ. Note: this run hit a `KeyError` because the input CSV was missing the `annotator_MED` column — check `data/Annotation_Indoor_Outdoor/annotations_indoor_outdoor.csv` has all 3 annotator columns before re-running. |
| `28062026_goldstandard_evaluation_indoor_outdoor_clip.ipynb` | Evaluates zero-shot CLIP predictions against the reconciled gold standard (300 images). Accuracy = 96.33%, F1 = 95.91%, Cohen's κ = 60.26%. Writes confusion matrix + agreement/disagreement CSVs to `data/Annotation_Indoor_Outdoor/evaluation_indoor_outdoor/`. |

The equivalent annotation/agreement work for the other two tasks (Human/Non-human,
Individual/CoEco/Landscape) is done directly inside `code/Classification_mod/` and
`code/Evalu_Class_mod/` notebooks rather than as separate notebooks here — see those folders'
READMEs.
