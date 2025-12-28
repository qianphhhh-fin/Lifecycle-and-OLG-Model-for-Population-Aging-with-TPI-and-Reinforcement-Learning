# 读取data/china_population_by_age.xlsx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 读取人口数据
population_data = pd.read_excel('china_population_by_age.xlsx')
population_data.rename(columns={'100+': 100}, inplace=True)
# 显示数据前几行
print("人口数据概览：")
print(population_data.head())

# 获取所有unique的年份
unique_years = population_data['Year'].unique()
print(f"数据包含年份：{unique_years}")

# 从22岁开始，每5岁一个年龄组
age_groups = list(range(22, 101, 5))

# 创建用于存储所有年份结果的DataFrame
all_population_data = pd.DataFrame()
all_normalized_data = pd.DataFrame()

# 对每个年份进行处理
for year in unique_years:
    print(f"处理年份：{year}")
    
    # 筛选当前年份的数据
    year_data = population_data[population_data['Year'] == year]
    
    # 按年龄组求和
    population_by_age_group = []
    for i in age_groups:
        try:
            population_by_age_group.append(year_data.loc[:, i:i+4].sum(axis=1).values[0])
        except:
            population_by_age_group.append(year_data.loc[:, i:].sum(axis=1).values[0])
            break
    
    # 创建当前年份的DataFrame
    year_column_name = f'y{year}'
    year_population_df = pd.DataFrame(population_by_age_group, index=age_groups, columns=[year_column_name])
    
    # 归一化人口数量使总和为1000
    year_normalized_df = pd.DataFrame(
        year_population_df[year_column_name] / year_population_df[year_column_name].sum() * 1000, 
        index=age_groups, 
        columns=[year_column_name]
    )
    
    # 合并到总的DataFrame中
    if all_population_data.empty:
        all_population_data = year_population_df
        all_normalized_data = year_normalized_df
    else:
        all_population_data = pd.concat([all_population_data, year_population_df], axis=1)
        all_normalized_data = pd.concat([all_normalized_data, year_normalized_df], axis=1)

# 保存到Excel文件的两个工作表
with pd.ExcelWriter('population_by_age_group_all_years.xlsx') as writer:
    # 重置索引并将年龄组列命名为AgeGroup
    all_population_data_reset = all_population_data.reset_index()
    all_population_data_reset.rename(columns={'index': 'AgeGroup'}, inplace=True)
    
    all_normalized_data_reset = all_normalized_data.reset_index()
    all_normalized_data_reset.rename(columns={'index': 'AgeGroup'}, inplace=True)
    
    all_population_data_reset.to_excel(writer, sheet_name='population', index=False)
    all_normalized_data_reset.to_excel(writer, sheet_name='pop_normalized', index=False)

print("数据已保存到 population_by_age_group_all_years.xlsx")
print("population sheet：原始人口数据")
print("pop_normalized sheet：归一化人口数据")

# 处理population_data，将数字列名前加上'age_'前缀
population_data_fixed = population_data.copy()
rename_dict = {}
for col in population_data.columns:
    if isinstance(col, (int, float)) or (isinstance(col, str) and col.replace('.', '', 1).isdigit()):
        rename_dict[col] = f'age_{col}'

population_data_fixed.rename(columns=rename_dict, inplace=True)
population_data_fixed.to_excel('china_population_by_age_headerFix.xlsx', index=False)
print("已将数字列名添加'age_'前缀并保存为'china_population_by_age_headerFix.xlsx'")
# 处理UN_PPP2024_Output_CBR.xlsx文件
print("\n开始处理UN_PPP2024_Output_CBR.xlsx文件...")

try:
    # 读取出生率数据文件，列名在第17行
    birth_rate_data = pd.read_excel('UN_PPP2024_Output_CBR.xlsx', header=16)
    
    print("出生率数据原始列名：")
    print(list(birth_rate_data.columns))
    
    # 查找包含China的行
    china_column = None
    for col in birth_rate_data.columns:
        if 'Region' in str(col) or 'country' in str(col) or 'area' in str(col):
            china_column = col
            break
    
    if china_column is None:
        print("未找到包含地区/国家信息的列")
        raise ValueError("未找到包含地区/国家信息的列")
    
    print(f"找到地区列：{china_column}")
    
    # 筛选中国数据
    china_data = birth_rate_data[birth_rate_data[china_column] == 'China']
    
    if china_data.empty:
        print("未找到中国数据")
        raise ValueError("未找到中国数据")
    
    print(f"找到中国数据行数：{len(china_data)}")
    
    # 识别年份列并创建新的列名映射
    new_column_names = {}
    year_columns = []
    
    for col in china_data.columns:
        # 检查列名是否包含4位数字的年份
        if re.search(r'\b(19|20|21)\d{2}\b', str(col)):
            # 提取年份数字
            year_match = re.search(r'\b(19|20|21)\d{2}\b', str(col))
            if year_match:
                year = year_match.group()
                new_column_names[col] = f'y{year}'
                year_columns.append(col)
    
    print(f"找到年份列数：{len(year_columns)}")
    print("年份列映射：")
    for old_name, new_name in new_column_names.items():
        print(f"  {old_name} -> {new_name}")
    
    # 只保留年份列
    china_year_data = china_data[year_columns].copy()
    
    # 重命名列
    china_year_data.rename(columns=new_column_names, inplace=True)
    
    print("重命名后的列名：")
    print(list(china_year_data.columns))
    
    # 保存为UN_PPP2024_CBR_birthper1000_China.xlsx
    china_year_data.to_excel('UN_PPP2024_CBR_birthper1000_China.xlsx', index=False)
    print("已保存处理后的UN_PPP2024_CBR_birthper1000_China.xlsx文件")
    
except Exception as e:
    print(f"处理UN_PPP2024_Output_CBR.xlsx文件时出错：{e}")
