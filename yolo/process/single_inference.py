from ultralytics import YOLO
import os
import pandas as pd
import cv2

# Load the pre-trained YOLO model
model = YOLO("best.pt")  

# Input image directory
image_dir = "../images"
image_paths = [os.path.join(image_dir, fname) for fname in os.listdir(image_dir) if fname.endswith(('.jpg', '.png', '.jpeg'))]

# Output directory to save results
output_image_dir = "../output_images"
os.makedirs(output_image_dir, exist_ok=True)  # 如果目录不存在，则创建

# Store detection results
all_results = []

# Process each image in the image directory
for image_path in image_paths:
    result = model(image_path)   # Run inference on a single image
    image_name = os.path.basename(image_path)   
    boxes = result[0].boxes  
    if boxes is not None:
        for box in boxes:
            # Convert Tensor data to standard Python types
            x1, y1, x2, y2 = map(int, box.xyxy[0].tolist()) 
            confidence = float(box.conf[0])  
            class_id = int(box.cls[0])  
            class_name = model.names[class_id]  
            all_results.append([image_name, x1, y1, x2, y2, confidence, class_id, class_name])
    
    # save images
    output_image_path = os.path.join(output_image_dir, image_name)
    annotated_image = result[0].plot()  
    cv2.imwrite(output_image_path, annotated_image)  

# Save the generated annotated image
df = pd.DataFrame(all_results, columns=["image_id", "x1", "y1", "x2", "y2", "confidence", "hard_score", "Class Name"])
new_df = df.iloc[:, :7]
output_csv_path = "../data/detection_results.csv"
new_df.to_csv(output_csv_path, index=False, encoding="utf-8-sig")

print(f"Detection results saved to {output_csv_path}")
print(f"Annotated images saved in {output_image_dir}")
