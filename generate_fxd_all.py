import sys
import json
import re
import rasterio
from rasterio.mask import mask
from shapely.geometry import shape, mapping, MultiPolygon, Polygon, Point
from shapely.ops import transform, nearest_points, unary_union
from pyproj import Proj, Transformer
from rasterstats import zonal_stats
import math

def filter_features(geojson_data):
    filtered_features = []
    for feature in geojson_data['features']:
        elev = feature['properties'].get('elev', None)
        if elev is None or elev < -0.5:
            filtered_features.append(feature)
    return filtered_features

def merge_polygons(geojson_data):
    polygons = []
    for feature in geojson_data['features']:
        geom = shape(feature['geometry'])
        if isinstance(geom, Polygon):
            polygons.append(geom)
        elif isinstance(geom, MultiPolygon):
            polygons.extend(geom.geoms)
    if polygons:
        merged = unary_union(polygons)
        if isinstance(merged, Polygon):
            merged = [merged]
        elif isinstance(merged, MultiPolygon):
            merged = merged.geoms
        return [mapping(p) for p in merged]
    return []

def get_properties(geojson_data):
    mmin_values = []
    mmax_values = []
    elev_values = []
    for feature in geojson_data['features']:
        properties = feature['properties']
        if 'mmin' in properties:
            mmin_values.append(properties['mmin'])
        if 'mmax' in properties:
            mmax_values.append(properties['mmax'])
        if 'elev' in properties:
            elev_values.append(properties['elev'])
    properties = {
        'mmin': min(mmin_values) if mmin_values else None,
        'mmax': max(mmax_values) if mmax_values else None,
        'elev': max(elev_values) if elev_values else None
    }
    return properties

def split_multipolygons(geojson_data):
    new_features = []
    for feature in geojson_data['features']:
        geom = shape(feature['geometry'])
        if isinstance(geom, MultiPolygon):
            for poly in geom.geoms:
                new_feature = feature.copy()
                new_feature['geometry'] = poly.__geo_interface__
                new_features.append(new_feature)
        elif isinstance(geom, Polygon):
            new_features.append(feature)
    return new_features

def generate_centroids(geojson_data):
    centroid_features = []
    for feature in geojson_data['features']:
        geom = shape(feature['geometry'])
        centroid = geom.centroid
        centroid_feature = {
            'type': 'Feature',
            'geometry': centroid.__geo_interface__,
            'properties': feature['properties']
        }
        centroid_features.append(centroid_feature)
    return centroid_features

def transform_geom(geom, transformer):
    return transform(transformer.transform, geom)

def transform_geojson(geojson_data, src_proj, dst_proj):
    transformer = Transformer.from_proj(src_proj, dst_proj, always_xy=True)
    for feature in geojson_data['features']:
        geom = shape(feature['geometry'])
        transformed_geom = transform_geom(geom, transformer)
        feature['geometry'] = transformed_geom.__geo_interface__
    return geojson_data

def find_nearest_line(centroid, lines):
    min_distance = float('inf')
    nearest_line = None
    for line in lines:
        line_geom = shape(line['geometry'])
        p1, p2 = nearest_points(centroid, line_geom)
        distance = p1.distance(p2)
        if distance < min_distance:
            min_distance = distance
            nearest_line = line
    return nearest_line

def extract_time_period_from_filename(filename):
    match = re.search(r'(\d{6}_\d{6})', filename)
    if match:
        return match.group(1)
    else:
        return ''

def add_nearest_line_info(centroids_data, lines_data, transformer_to_4326, time_period):
    lines = lines_data['features']
    for feature in centroids_data['features']:
        centroid = shape(feature['geometry'])
        nearest_line = find_nearest_line(centroid, lines)
        if nearest_line:
            feature['properties']['毗邻断面'] = nearest_line['properties'].get('section_name', '')

        # 将中心点转换回4326坐标系以获取经纬度
        centroid_4326 = transform(transformer_to_4326.transform, centroid)
        lon, lat = centroid_4326.x, centroid_4326.y
        feature['properties']['风险点中心坐标'] = [lon, lat]

        # 添加比较时段
        feature['properties']['比较时段'] = time_period

    return centroids_data

def calculate_stats(polygons_geojson, dem_file, dem_crs):
    src_proj = Proj('epsg:4326')
    dst_proj = Proj(dem_crs)
    transformer = Transformer.from_proj(src_proj, dst_proj, always_xy=True).transform

    stats_list = []
    for feature in polygons_geojson['features']:
        geom = shape(feature['geometry'])
        geom_transformed = transform(transformer, geom)

        with rasterio.open(dem_file) as src:
            out_image, out_transform = mask(src, [mapping(geom_transformed)], crop=True, nodata=-999)
            out_image = out_image[0]

            # 使用rasterstats计算统计数据
            stats = zonal_stats([geom_transformed], out_image, affine=out_transform, nodata=-999, stats=["mean", "min"])

            area = geom_transformed.area  # 面积，单位为平方米

            stats_list.append({
                "最大冲深(m)": stats[0]["min"],
                "平均冲深(m)": stats[0]["mean"],
                "风险面积(m^2)": area
            })

    return stats_list

def update_centroids_with_stats(input_centroids_file, stats, output_centroids_file):
    with open(input_centroids_file, 'r', encoding='utf-8') as f:
        centroids_geojson = json.load(f)

    for centroid_feature, stat in zip(centroids_geojson['features'], stats):
        centroid_feature['properties'].update(stat)

    with open(output_centroids_file, 'w', encoding='utf-8') as f:
        json.dump(centroids_geojson, f, ensure_ascii=False, indent=4)

def save_geojson(data, filename):
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

def main(input_polygons_file, input_lines_file, input_dem_file, output_polygons_file, output_centroids_file):
    # 步骤一：处理面要素
    with open(input_polygons_file, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)

    filtered_features = filter_features(geojson_data)
    geojson_data['features'] = filtered_features

    merged_polygons = merge_polygons(geojson_data)

    merged_geojson = {
        'type': 'FeatureCollection',
        'features': [{'type': 'Feature', 'geometry': geom, 'properties': get_properties(geojson_data)} for geom in merged_polygons]
    }

    splitted_features = split_multipolygons(merged_geojson)
    geojson_data['features'] = splitted_features

    save_geojson(geojson_data, output_polygons_file)

    # 步骤二：处理中心点和毗邻断面、风险点中心坐标、比较时段
    centroids_data = {
        'type': 'FeatureCollection',
        'features': generate_centroids(geojson_data)
    }

    src_proj = Proj('epsg:4326')
    dst_proj = Proj('epsg:3857')

    time_period = extract_time_period_from_filename(input_polygons_file)

    with open(input_lines_file, 'r', encoding='utf-8') as f:
        lines_data = json.load(f)

    lines_data_transformed = transform_geojson(lines_data, src_proj, dst_proj)
    centroids_data_transformed = transform_geojson(centroids_data, src_proj, dst_proj)

    transformer_to_4326 = Transformer.from_proj(dst_proj, src_proj, always_xy=True)

    centroids_data_updated = add_nearest_line_info(centroids_data_transformed, lines_data_transformed,
                                                   transformer_to_4326, time_period)

    centroids_data_final = transform_geojson(centroids_data_updated, dst_proj, src_proj)

    save_geojson(centroids_data_final, output_centroids_file)

    # 步骤三：处理统计数据并更新中心点GeoJSON
    with rasterio.open(input_dem_file) as dem_src:
        dem_crs = dem_src.crs.to_string()

    stats = calculate_stats(geojson_data, input_dem_file, dem_crs)
    update_centroids_with_stats(output_centroids_file, stats, output_centroids_file)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python generate_fxd_all.py <input_polygons_geojson> <input_lines_geojson> <input_dem_tif> <output_polygons_geojson> <output_centroids_geojson>")
        sys.exit(1)

    input_polygons_file = sys.argv[1]
    input_lines_file = sys.argv[2]
    input_dem_file = sys.argv[3]
    output_polygons_file = sys.argv[4]
    output_centroids_file = sys.argv[5]
    main(input_polygons_file, input_lines_file, input_dem_file, output_polygons_file, output_centroids_file)
