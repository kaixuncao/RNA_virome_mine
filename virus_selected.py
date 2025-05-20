import os
import pandas as pd

# 步骤 1: 读取当前文件夹下的RdRp_DB_Read_length.txt
db_read_length_file = "RdRp_DB_Read_length.txt"
df_db_read_length = pd.read_csv(db_read_length_file, sep="\t", header=None)

# 获取当前文件夹下的所有子文件夹
subfolders = [f for f in os.listdir('.') if os.path.isdir(f)]

# 步骤 2: 检索每个文件夹
for folder in subfolders:
    rdrp_list_file = os.path.join(folder, folder + "_rdrp_list.txt")
    # 检查文件是否存在
    if os.path.isfile(rdrp_list_file):
        # 步骤 3: 读取文件，根据第二列与RdRp_DB_Read_length.txt中第一列匹配的内容
        df_rdrp_list = pd.read_csv(rdrp_list_file, sep="\t", header=None)
        # 合并两个DataFrame，条件是RdRp_DB_Read_length.txt的第一列与_rdrp_list.txt的第二列相匹配
        merged_df = pd.merge(df_db_read_length, df_rdrp_list, left_on=0, right_on=1)

        # 计算 y/x 的比值，保留比值大于0.9的行
        merged_df['y/x'] = merged_df.iloc[:, 6] / merged_df.iloc[:, 2]
        result_df = merged_df[merged_df['y/x'] > 0.9]

        # 选取所需的列
        result_df = result_df.iloc[:, [3, 4, 9, 10]]
        result_df['y/x'] = merged_df['y/x']

        # 保存为当前文件夹名加 "_covmatch.txt" 的文件
        result_df.to_csv(os.path.join(folder, folder + "_covmatch.txt"), sep="\t", index=False, header=False)

print("Process completed.")
