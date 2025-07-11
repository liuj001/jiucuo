import argparse
from ultralytics import YOLO
import os
import shutil

# 获取当前脚本所在的目录
script_dir = os.path.dirname(os.path.abspath(__file__))

# 解析命令行参数
parser = argparse.ArgumentParser(description='Train YOLOv8 model with custom parameters.')
parser.add_argument('-lr', type=float, default=0.001, help='Initial learning rate')
parser.add_argument('-b', type=int, default=1, help='Batch size for training')
parser.add_argument('-epochs', type=int, default=4, help='Number of epochs to train')
parser.add_argument('-imgsz', type=int, default=1471, help='Image size for training')
parser.add_argument('-data', type=str, default=os.path.join(script_dir, "data.yaml"), help='Path to data configuration file')  # 关键修改
parser.add_argument('-model', type=str, default=os.path.join(script_dir, "../process/best.pt"), help='Path of pretrained model')
args = parser.parse_args()

# 加载模型
model = YOLO(args.model)

# 训练模型
results = model.train(
    data=args.data,  # 现在 data.yaml 的路径是绝对路径
    epochs=args.epochs,
    imgsz=args.imgsz,
    lr0=args.lr,
    batch=args.b,
    amp=False,
    project="runs",
    name="weights",
    exist_ok=True
)

# 移动权重文件
src_dir = results.save_dir
dst_dir = os.path.join(script_dir, "../process")
os.makedirs(dst_dir, exist_ok=True)

weight_file = "best.pt"
src_path = os.path.join(src_dir, "weights", weight_file)
if os.path.exists(src_path):
    shutil.copy(src_path, dst_dir)