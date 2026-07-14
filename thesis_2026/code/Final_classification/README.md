# code/Final_classification/

Applies the **best-performing model + method per task** (as determined in `Evalu_Class_mod/`) to
the full, large-scale image corpus. This is the production run of the classification cascade
described in the [root README](../../README.md).

| Notebook | Task / stage | Model + method | Input images | Output distribution |
|---|---|---|---|---|
| `28062026_class_indoor_outdoor_mobileclips1_simple_prompts.ipynb.ipynb` | Level 1: Indoor/Outdoor | MobileCLIP-S1, simple prompts | 300-image check run | 14 indoor / 286 outdoor |
| `28062026_Final_class_mobile_clip_human_non_human_linear_probe.ipynb` | Level 2: Human/Non-human | MobileCLIP-S1 + linear probe | 8,219 outdoor photos | 3,125 human / 5,094 non-human |
| `30062026_Final_class_vitl14_human_non_human_linear_probe.ipynb` | Level 2: Human/Non-human (alt.) | CLIP ViT-L/14 + linear probe | 8,219 outdoor photos | 2,933 human / 5,286 non-human |
| `01072026_final_class_individual_coeco_landscape_mobileclip_mixed.ipynb.ipynb` | Level 3: Individual/CoEco/Landscape | MobileCLIP-S1, mixed prompts (zero-shot) | 5,240 non-human photos | 992 community_ecosystem / 562 individual / 3,686 landscape |
| `01072026_final_class_individual_coeco_landscape_mobileclip_linearprobe.ipynb.ipynb` | Level 3: Individual/CoEco/Landscape (alt.) | MobileCLIP-S1 + linear probe | 5,240 non-human photos | 1,301 community_ecosystem / 517 individual / 3,422 landscape |

Two variants were run in parallel for Levels 2 and 3 (zero-shot/mixed-prompt vs. linear-probe) so
the results can be compared — see `data/Classification_Human_non_Human/` and
`data/Classification_Individual_CoEco_Landscape/` for the resulting sorted image folders and CSVs.
The `landscape` output feeds directly into `code/PUD/`.
