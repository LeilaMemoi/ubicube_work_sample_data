# %% [code] {"jupyter":{"outputs_hidden":false},"execution":{"iopub.status.busy":"2024-11-13T23:31:11.311882Z","iopub.execute_input":"2024-11-13T23:31:11.312343Z","iopub.status.idle":"2024-11-13T23:31:34.049797Z","shell.execute_reply.started":"2024-11-13T23:31:11.312287Z","shell.execute_reply":"2024-11-13T23:31:34.048559Z"}}

import tifffile as tiff
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from skimage.segmentation import slic
from skimage.segmentation import mark_boundaries
from skimage.measure import label, regionprops, regionprops_table
from skimage.measure import regionprops
from skimage.util import map_array
import rasterio
from rasterio.features import geometry_mask
import rasterio.mask
import rasterstats
from rasterio.mask import mask

# Read the data
def read_data(rgb_path, nir_path, parcel_path, building_path):
    # Read the raster data
    with rasterio.open(rgb_path) as src:
        rgb = src.read()  # Read all bands (rows, columns)
        rgb = np.transpose(rgb, (1, 2, 0))  # Change the order of the bands to (rows, columns, bands)

    with rasterio.open(nir_path) as src:
        nir = src.read(1) 

    # Read the vector data
    parcels = gpd.read_file(parcel_path)
    buildings = gpd.read_file(building_path)
    
    # Ensure the data have the same coordinate system
    buildings = buildings.to_crs(src.crs)
    parcels = parcels.to_crs(src.crs)

    return rgb, nir, parcels, buildings, src

# Function to normalize image bands
def normalize_bands(red, green, blue, nir):
    red_norm = (red - np.min(red)) / (np.max(red) - np.min(red))
    green_norm = (green - np.min(green)) / (np.max(green) - np.min(green))
    blue_norm = (blue - np.min(blue)) / (np.max(blue) - np.min(blue))
    nir_norm = (nir - np.min(nir)) / (np.max(nir) - np.min(nir))
    return red_norm, green_norm, blue_norm, nir_norm

# Function to compute NDVI
def compute_ndvi(nir_norm, red_norm):
    denominator_zero_mask = (nir_norm + red_norm) == 0
    ndvi = np.where(denominator_zero_mask, 0, (nir_norm - red_norm) / (nir_norm + red_norm))
    return ndvi

# Object-based image segmentation, using the SLIC algorithm
def image_segmentation(image, algorithm='slic', algorithm_params=None, channel_axis=-1):
    if algorithm_params is None:
        algorithm_params = {}

    # Perform segmentation
    segments_result = slic(image, **algorithm_params, channel_axis=channel_axis)
    
    # Mark boundaries on the original image using the segmented result
    boundaries = mark_boundaries(image[:, :, [0, 1, 2]], segments_result, (1, 0, 0), mode="thick")

    return segments_result, boundaries

# Function to plot the segments
def plot_segments(image, segments, boundaries):
    # Display the original image and the segmentation result
    fig, ax = plt.subplots(1, 3, figsize=(12, 6))

    ax[0].imshow(image[:, :, [0, 1, 2]])
    ax[0].set_title(f'Tile')

    ax[1].imshow(segments)
    ax[1].set_title(f'Segments_Tile')

    ax[2].imshow(boundaries)
    ax[2].set_title(f'Segments with boundaries')

    plt.show()

# Function to map NDVI values back to the image segments
def map_segment_ndvi(segments, seg_ndvi):
    mapped_ndvi = map_array(segments, np.array(seg_ndvi['label']), np.array(seg_ndvi['intensity_mean']))
    return mapped_ndvi

def process_data(rgb_path, nir_path, parcel_path, building_path, output_path):
    rgb, nir, parcels, buildings, src = read_data(rgb_path, nir_path, parcel_path, building_path)
    
    # Extract bands from RGB image
    red = rgb[:, :, 0]
    green = rgb[:, :, 1]
    blue = rgb[:, :, 2]

    # Normalize the bands
    red_norm, green_norm, blue_norm, nir_norm = normalize_bands(red, green, blue, nir)

    # Compute NDVI
    ndvi = compute_ndvi(nir_norm, red_norm)

    # Perform slic image segmentation
    slic_params = {'n_segments': 1000, 'compactness': 0.5}
    segments, boundaries = image_segmentation(rgb, algorithm='slic', algorithm_params=slic_params)

    # Calculate NDVI values for each segment
    seg_ndvi = calc_segment_ndvi(segments, ndvi)

    # Map NDVI values back to segments
    mapped_ndvi = map_segment_ndvi(segments, seg_ndvi)

    # Define NDVI range for pervious surfaces (vegetated)
    ndvi_min = -1
    ndvi_max = 0.1
    # Eliminate these pervious surfaces
    impervious_surfaces = np.where((mapped_ndvi >= ndvi_min) & (mapped_ndvi < ndvi_max), mapped_ndvi, np.nan)
    
    # Create a building mask where pixels within the building are 0 and the rest 1
    building_mask = geometry_mask([geom for geom in buildings.geometry], transform=src.transform, out_shape=image.shape[1:])

    # Mask the impervious surfaces array with the building mask
    impervious_surfaces_no_buildings = np.where(building_mask, impervious_surfaces, np.nan)
    
    impervious_mask = impervious_surfaces > 0  # Binary mask for impervious areas

    # Calculate zonal statistics - count of impervious pixels within each parcel polygon
    stats = rasterstats.zonal_stats(parcels, impervious_surfaces_no_buildings, affine=src.transform, stats="count", nodata=src.nodata)

    # Add the impervious area (count of impervious pixels) to the GeoDataFrame
    parcels["imperv_pixels"] = [s['count'] if s['count'] is not None else 0 for s in stats]
    parcels["imperv_area"] = parcels["imperv_pixels"] * 0.4  # To compute area in real world units(m2) by multiplying by a pixel area (0.2m resolution)

    # Calculate impervious area percentage
    parcels["unsealing_potential"] = (parcels["imperv_area"] / parcels.area) * 100
    
    # Plot the parcels showing unsealing potential
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    parcels.plot(column="imperv_area", cmap="OrRd", legend=True, ax=ax)
    plt.title("Unsealing Potential of Parcels")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.show()

    # Save the updated GeoDataFrame with unsealing potential to a new shapefile
    parcels.to_file(output_path)

    print(f"Shapefile with unsealing potential saved to {output_path}")

# Input file paths
rgb_path = '/kaggle/input/ubicube-work-sample/aerial/AOI_ID_1_ortho_cog_RGB.tif'
nir_path = '/kaggle/input/ubicube-work-sample/aerial/AOI_ID_1_ortho_cog_NIR.tif'
bld_path = '/kaggle/input/ubicube-work-sample/building_footprints/building_footprints.geojson'
par_path = '/kaggle/input/ubicube-work-sample/parcels/parcels.shp'

# Output path (for saving updated parcels unsealing potential shapefile )
output_path = "/kaggle/working/parcels_unsealing_potential.shp"

# Run the processing function
process_data(rgb_path, nir_path, par_path, bld_path, output_path)