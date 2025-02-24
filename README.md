# Automated Workflow for Assessing Urban Unsealing Potential
A Python script that calculates the "unsealing potential" (sealed soil without buildings, measured in m²) for individual real estate properties (land parcels). The script enables ranking properties within a given geographical area based on their unsealing potential, allowing for the identification of the most promising properties for soil unsealing. This ranking can support the design of targeted measures to reduce soil sealing effectively.

### Data

Aerial imagery: 20cm spatial resolution, RGB + NIR.

Building Footprints: OpenStreetMap data.

Parcel/property boundaries

## Workflow overview
### 1. Data Ingestion 

Reads RGB and NIR aerial imagery using rasterio.

Loads vector data (parcel boundaries and building footprints) using geopandas.

Ensures coordinate reference systems (CRS) are consistent across datasets.

### 2. Preprocessing

Extracts individual bands (Red, Green, Blue, and NIR) from the imagery.

Normalizes each band for consistent analysis.

Computes the Normalized Difference Vegetation Index (NDVI) to differentiate vegetated and impervious surfaces.

### 3. Image Segmentation

Applies the SLIC (Simple Linear Iterative Clustering) algorithm to segment the image into meaningful regions.

Generates boundary-marked visualizations of the segmented image.

### 4. Feature Extraction

Calculates NDVI statistics for each segment.

Maps segment-level NDVI values back to the image.

Identifies impervious surfaces by setting NDVI thresholds.

### 5. Building Masking

Generates a building mask to exclude built-up areas from impervious surface analysis.

Applies the mask to remove buildings from the impervious surface map.

### 6. Zonal Statistics & Unsealing Potential Calculation

Uses rasterstats to compute the count of impervious pixels within each parcel.

Converts pixel counts into real-world area measurements.

Calculates the percentage of unsealing potential (impervious area as a proportion of parcel size).

### 7.  Visualization & Output

Plots the parcels with unsealing potential highlighted using a color-coded map.

Saves the final GeoDataFrame with calculated unsealing potential as a shapefile.

## Usage

### Required Libraries:

Ensure the following Python libraries are installed:

```bash
pip install rasterio rasterstats tifffile numpy pandas geopandas matplotlib scikit-image
```

### Running the Script
The script takes the paths to the image data, building footprints, and parcel boundaries as input and writes the results (unsealing potential in m²) to a shapefile.
```bash
process_data(rgb_path, nir_path, par_path, bld_path, output_path)
```
Modify file paths accordingly and execute the script:
