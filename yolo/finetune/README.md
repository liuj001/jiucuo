### README

### INTRODUCTION

This is a model using YOLOv8 to detect adapter-like features in the field of biology. The model takes self-generated images as input, which contain SV and adapter features. It detects all targets to prepare for the subsequent screening of adapters.

### Environment

python>=3.8

` pip install -e .`

### DATASET

In the project, we have provided a small portion of the dataset that has been annotated. You need to modify the `data.yaml` file to configure it.

```yaml
train: /root/autodl-tmp/ada/train/data/adapter/images/train
test:  /root/autodl-tmp/ada/train/data/adapter/images/test
val:   /root/autodl-tmp/ada/train/data/adapter/images/val


nc: 1 # only 1 class
names: ['object'] 

```

### run

`python train.py`

```
from ultralytics import YOLO

# Load a COCO-pretrained YOLOv8n model
model = YOLO("yolov8n.pt")

# Display model information (optional)
model.info()

# Train the model on the COCO8 example dataset for 100 epochs
results = model.train(data="data.yaml", epochs=400, imgsz=1471)
```

you can use yolov8x.pt to train a better model.
