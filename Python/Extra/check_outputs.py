""" Test the integrity of processed data
    May also generate metadata; quicklook images etc
    old name: test_outputs.py
"""
import os
import datetime
# import xarray
import rioxarray
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import logging
"""
# ##############################################################################
from sentinel_hydrogen import check_geotiff_integrity
log_name = 'file_integrity_check.log'
logging.basicConfig(filename=f'{log_name}', level=logging.INFO, format='%(asctime)s:%(message)s')
logging.info('  ')
logging.info('-------------------------START LOGGING--------------------------------------')

# List of tiles that cover North America (our AOI)
modis_tile_list = [
    'h07v03', 'h07v05', 'h07v06', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h10v02', 'h10v03', 
    'h10v04', 'h10v05', 'h10v06', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h13v01', 
    'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v01', 'h16v02']
logging.info(f"Number modis tiles required: {len(modis_tile_list)}")

# 0. Check the integrity of downloaded files
root_dir = "C:"  # "/mnt/c"
base_folder = f"{root_dir}/Github/coressd"
base_folder = "/discover/nobackup/projects/coressd"
modis_download_folder = f"{base_folder}/OSU/MOD10A1F.061/MODIS_Proc/download_snow"  # /2016001/001
# # Generate the date range for WY2016 and convert to DOY format of MODIS naming convention
# start_date = datetime.datetime.strptime("2015-10-01", "%Y-%m-%d")
# end_date = datetime.datetime.strptime("2016-09-30", "%Y-%m-%d")
# date_generated = pd.date_range(start_date, end_date)
# year_doy_list = list(date_generated.strftime("%Y%j"))
year_doy_list = os.listdir(modis_download_folder)
for year_doy in year_doy_list:
    logging.info(year_doy)
    year = year_doy[:4]  # "2016"
    doy = year_doy[4:]  # "215"
    download_folder = f"{modis_download_folder}/{year}{doy}/{doy}"
    files = [f for f in os.listdir(download_folder) if f.endswith(".hdf") and f.startswith("MOD10A1F")]
    files = [f for f in files if f.split(".")[2] in modis_tile_list]
    # for tif in files:
    #     check_geotiff_integrity(f'{download_folder}/{tif}')
    _ = Parallel(n_jobs=-1)(delayed
        (check_geotiff_integrity)
        (f'{download_folder}/{tif}')
        for tif in files)
logging.info('Finished checking integrity')
# ##############################################################################
"""


# Check if sinusoidal projection from all regions are same: else we cannot merge
# 1. Check the MODIS_CGF clipped to SEUP_NorthAmerica Resolution
# 1.1. Resolution
# def region_metadata_csv(aoi, band='B02', save_csv=False):
def region_metadata_csv(save_csv=True):
    """ Create metdata csv for each region
        Get bounding coordinates, extent, and resolution of each raster
        We will use this csv later to troubleshoot.
        NB:
        ---
            Changed region to aoi [dec 31, 2021]
    """
    #for aoi in aoi_list[:1]:
    # Get Clipped TIFFs for a region
    #tifs = os.listdir(f'/{base_folder}/clipped/{region}') 
    clipped_folder = '/discover/nobackup/projects/coressd/Blender/Modis/CGF_NDSI_Snow_Cover/NA_mosaic'
    # clipped_folder = '/discover/nobackup/projects/coressd/Blender/Modis/MOD44B/Percent_Tree_Cover/NA_mosaic'
    tifs = os.listdir(clipped_folder) 
    tifs.sort()
    # print(aoi, len(tifs))
    # Get some shape attributes from each raster to analyze what is going on
    coord_list = []
    tif_list = []
    utm_list = []
    shape_list = []
    res_list = [] #resolution; if all 10m or not
    kb_list = [] #raster size in kilobytes
    for tif in tifs:
        # try:
        ds = rioxarray.open_rasterio(f'{clipped_folder}/{tif}').squeeze()
        # ds = ds.sel(band=1)
        # ds = ds["CGF_NDSI_Snow_Cover"]
        kb = ds.nbytes/1e3
        kb_list.append(kb)
        coord_list.append(ds.rio.bounds())
        shape_list.append(ds.shape)
        res_list.append(ds.rio.resolution())
        tif_name = tif.split('.nc')[0]
        tif_list.append(tif_name)
        # utm_list.append(tif_name.split('_')[-1])
        # except:
        #     print(f'Error/Missing: {clipped_folder}/{tif}')
    coord_list = np.array(coord_list) #convert tuple to numpy array so we can stack them horizontally
    shape_list = np.array(shape_list)
    res_list = np.array(res_list)
    # utm_list = np.array(utm_list)
    coord_arr = np.hstack((coord_list, shape_list, res_list))
    # Convert to Dataframe
    # df = pd.DataFrame([shape_list, coord_list], index=tif_list)
    df = pd.DataFrame(coord_arr, index=tif_list)
    df.columns=['x0', 'y0', 'x1', 'y1', 'sizey', 'sizex', 'resx', 'resy']
    df['kb'] = kb_list
    # df['utm'] = utm_list
    df.index.name = 'Tile'
    if save_csv:
        # df.to_csv('/discover/nobackup/projects/coressd/Blender/Modis/tree_mosaic_meta.csv')  # creation of clipped folder defined below; not necessary for data integrity check
        df.to_csv('/discover/nobackup/projects/coressd/Blender/Modis/cgf_snow_meta.csv')  # creation of clipped folder defined below; not necessary for data integrity check
    return df


region_metadata_csv()

# 1.2 x and y coordinate values 

# 1.3 projection/crs

# 1. 

""" Extra helpful codes
dem.rio.resolution()

# Double check if size matches
for idx in idx_list:
    tif_name = os.listdir(f'{s2clip_folder}/{idx}')[0]
    s2 = rioxarray.open_rasterio(f'{s2clip_folder}/{idx}/{tif_name}/{tif_name}_B08.tif').squeeze()
    dem = rioxarray.open_rasterio(f'{clip_folder}/{idx}.tif').squeeze()
    if s2.x.size!=dem.x.size:
        print(idx, s2.x.size, dem.x.size, s2.y.size, s2.y.size)
    if s2.y.size!=dem.y.size:
        print(idx, s2.x.size, dem.x.size, s2.y.size, s2.y.size)


"""