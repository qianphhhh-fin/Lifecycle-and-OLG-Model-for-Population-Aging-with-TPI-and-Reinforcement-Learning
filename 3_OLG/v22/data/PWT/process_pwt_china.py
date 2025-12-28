import pandas as pd
import os

# 获取脚本所在目录并构建文件路径
script_dir = os.path.dirname(os.path.abspath(__file__))
excel_file = os.path.join(script_dir, 'pwt1001.xlsx')

print(f"正在读取文件: {excel_file}")

# 读取Excel文件，第一行作为列名
df = pd.read_excel(excel_file, header=0, index_col=None,sheet_name='Data')

print(f"原始数据形状: {df.shape}")
print(f"列名: {list(df.columns)}")

# 筛选country=China的数据
china_data = df[df['country'] == 'China'].copy()

print(f"中国数据行数: {len(china_data)}")

if len(china_data) == 0:
    print("警告：未找到country='China'的数据，尝试查看country列的唯一值...")
    print("Country列的唯一值（前20个）:")
    print(df['countryy'].unique()[:20])
    
    # 尝试查找包含China的行
    china_like = df[df['country'].str.contains('China', case=False, na=False)]
    if len(china_like) > 0:
        print(f"找到包含'China'的数据: {len(china_like)}行")
        print("相关country值:")
        print(china_like['country'].unique())
        china_data = china_like.copy()

# 检查所需列是否存在
required_columns = ['year','rtfpna', 'cn', 'cgdpo'] # 年份，全要素生产率，资本存量，国内生产总值
missing_columns = [col for col in required_columns if col not in df.columns]

if missing_columns:
    print(f"警告：以下列不存在：{missing_columns}")
    print("可用列名:")
    print([col for col in df.columns if any(keyword in col.lower() for keyword in ['year','rtfpna', 'capital', 'gdp', 'output'])])

# 提取所需列（如果存在的话）
available_columns = ['country'] + [col for col in required_columns if col in df.columns]



if len(china_data) > 0:
    # 提取中国数据的指定列
    china_selected = china_data[available_columns].copy()

    china_selected['K-Y']=china_selected['cn']/china_selected['cgdpo'] # 资本产出比
    
    print(f"提取的列: {available_columns}")
    print(f"数据预览:")
    print(china_selected.head())
    
    # 保存为新的xlsx文件
    output_file = os.path.join(script_dir, 'china_pwt_data.xlsx')
    china_selected.to_excel(output_file, index=False)
    
    print(f"数据已保存到: {output_file}")
    print(f"保存的数据形状: {china_selected.shape}")
else:
    print("错误：未找到中国的数据") 