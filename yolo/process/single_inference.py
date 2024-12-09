from ultralytics import YOLO
import os
import pandas as pd
import cv2

# 加载预训练的 YOLO 模型
model = YOLO("best.pt")  # 预训练 YOLOv8n 模型

# 输入图像目录
image_dir = "../images"
image_paths = [os.path.join(image_dir, fname) for fname in os.listdir(image_dir) if fname.endswith(('.jpg', '.png', '.jpeg'))]

# 输出结果保存目录
output_image_dir = "../output_images"
os.makedirs(output_image_dir, exist_ok=True)  # 如果目录不存在，则创建

# 存储检测结果
all_results = []

# 遍历图像路径逐张处理
for image_path in image_paths:
    result = model(image_path)  # 对单张图像运行推断
    image_name = os.path.basename(image_path)  # 获取图片文件名
    boxes = result[0].boxes  # 获取边界框对象（取第一个结果，因为每次只处理一张图像）
    if boxes is not None:
        for box in boxes:
            # 将 Tensor 数据转换为标准 Python 类型
            x1, y1, x2, y2 = map(int, box.xyxy[0].tolist())  # 边界框的左上角和右下角坐标
            confidence = float(box.conf[0])  # 检测置信度
            class_id = int(box.cls[0])  # 类别ID
            class_name = model.names[class_id]  # 类别名称
            all_results.append([image_name, x1, y1, x2, y2, confidence, class_id, class_name])
    
    # 保存生成的检测图片
    output_image_path = os.path.join(output_image_dir, image_name)
    annotated_image = result[0].plot()  # 返回绘制结果的图像 (numpy array)
    cv2.imwrite(output_image_path, annotated_image)  # 使用 OpenCV 保存图像

# 将结果保存为 DataFrame
df = pd.DataFrame(all_results, columns=["image_id", "x1", "y1", "x2", "y2", "confidence", "hard_score", "Class Name"])
new_df = df.iloc[:, :7]
# 保存为 CSV 文件
output_csv_path = "../data/detection_results.csv"
new_df.to_csv(output_csv_path, index=False, encoding="utf-8-sig")

print(f"Detection results saved to {output_csv_path}")
print(f"Annotated images saved in {output_image_dir}")
