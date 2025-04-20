import argparse
from ultralytics import YOLO

# 解析命令行参数
parser = argparse.ArgumentParser(description='Train YOLOv8 model with custom parameters.')
parser.add_argument('-lr', type=float, default=0.001, help='Initial learning rate')
parser.add_argument('-b', type=int, default=1, help='Batch size for training')
parser.add_argument('-epochs', type=int, default=400, help='Number of epochs to train')
parser.add_argument('-imgsz', type=int, default=1471, help='Image size for training')
parser.add_argument('-data', type=str, default="data.yaml", help='Path to data configuration file')
args = parser.parse_args()

# Load a COCO-pretrained YOLOv8n model
model = YOLO("yolov8n.pt")

# Display model information (optional)
model.info()

# Train the model with parameters from command line
results = model.train(
    data=args.data,
    epochs=args.epochs,
    imgsz=args.imgsz,
    lr0=args.lr,
    batch=args.b
)