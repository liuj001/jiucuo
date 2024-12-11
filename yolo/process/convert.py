import pandas as pd


df = pd.read_csv('../result/adapter.csv')


selected_columns = ['image_id', 'x1', 'y1', 'x2', 'y2']
new_df = df[selected_columns].copy()


new_df['confidence'] = 0
new_df['hard_score'] = 0


new_df.to_csv('../result/adapter.csv', index=False)