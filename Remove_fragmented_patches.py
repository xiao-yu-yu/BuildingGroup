import geopandas as gpd
import pandas as pd

def filter_elements(element_a, element_b, min_required):
    # 创建一个新的 GeoDataFrame 以存储满足条件的 B 元素
    filtered_elements_b = gpd.GeoDataFrame(columns=element_b.columns)

    # 遍历 B 元素
    for index_b, geom_b in element_b.iterrows():
        # 统计 B 元素中包含在 A 元素中的数量
        count = 0
        for index_a, geom_a in element_a.iterrows():
            if geom_b['geometry'].contains(geom_a['geometry']):
                count += 1

        # 如果数量大于等于 min_required，则将 B 元素添加到新的 GeoDataFrame 中
        if count >= min_required:
            filtered_elements_b = pd.concat([filtered_elements_b, element_b.loc[[index_b]]], ignore_index=True)

    return filtered_elements_b

# 读取 A 和 B 的 Shapefile
file_path_a = 'E:/建筑物合并/不同比例尺矢量数据集/geball_EliminatePolygonPart1_CreateBuffers1.shp'
file_path_b = 'E:/建筑物合并/不同比例尺矢量数据集/geball25_AggregatePolygons.shp'

element_a = gpd.read_file(file_path_a)
element_b = gpd.read_file(file_path_b)

# 指定每个 B 元素需要包含的最小 A 元素数量
min_required_count = 10

# 过滤 B 元素
filtered_b = filter_elements(element_a, element_b, min_required_count)

# 打印结果
print("Filtered Element B:")
print(filtered_b)

# 将结果保存到新的 Shapefile 文件
filtered_b.to_file('E:/建筑物合并/不同比例尺矢量数据集/Filtered_B.shp')
