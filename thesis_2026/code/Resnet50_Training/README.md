# code/Resnet50_Training/

| Notebook | What it does |
|---|---|
| `Resnet50_training_indoor_outdoor.ipynb` | Fine-tunes ResNet-50 (ImageNet-pretrained, frozen backbone + dense head) on 5,000 balanced indoor/outdoor images (2,500 each) using 5-fold cross-validation, data augmentation (rotation/flip/shift/zoom) and mixed-precision GPU training. Writes 5 trained models (`.keras`), per-fold confusion matrices, loss curves and a metrics CSV. Average across folds: Accuracy ≈ 75–76%, F1 ≈ 0.73–0.75. |

Serves as the **supervised CNN baseline** for the Indoor/Outdoor task, evaluated in
`code/Evalu_Class_mod/Eva_Resnet50.ipynb` against the CLIP-based zero-shot/linear-probe approaches.
