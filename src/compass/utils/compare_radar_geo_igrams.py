import argparse
import json
import glob
import isce3
import os
import numpy as np
from osgeo import gdal
from s1reader.s1_reader import load_bursts
from s1reader.s1_orbit import get_orbit_file_from_list
from compass.utils.geo_grid import generate_geogrids

# Paths to directories
data_path = '/mnt/aurora-r0/vbrancat/data/S1/data/Rosamond'
stack_path = '/mnt/aurora-r0/vbrancat/data/S1/stack_processor/Rosamond'
orbit_dir = f'{data_path}/orbits'
dem_path = f'{data_path}/dem_4326.tiff'
radar_path = f'{stack_path}/radar'
geo_path = f'{stack_path}/geocoded/geo_bursts'


def command_line_parser():
    parser = argparse.ArgumentParser(description="""
                                      Compare geocoded and radar igrams""",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r', '--ref-date', type=str, dest='ref_date',
                        help='Reference date')
    parser.add_argument('-s', '--sec-date', type=str, dest='sec_date',
                        help='Secondary date')
    parser.add_argument('-b', '--burst-id', type=str, dest='burst_id',
                        help='burst id')
    parser.add_argument('-p', '--pol', type=str, dest='pol', default='VV',
                        help='burst id')
    parser.add_argument('-o', '--outdir', type=str, dest='outdir',
                        help='Output directory')
    return parser.parse_args()


def create_igram(ref_path, sec_path, ref_date,
                 sec_date, burst_id, outdir,
                 is_geocoded=True):

    ref_ds = gdal.Open(ref_path, gdal.GA_ReadOnly)
    ref = ref_ds.GetRasterBand(1).ReadAsArray()

    sec_ds = gdal.Open(sec_path, gdal.GA_ReadOnly)
    sec = sec_ds.GetRasterBand(1).ReadAsArray()
    igram = ref * np.conj(sec)
    length, width = igram.shape
    if is_geocoded:
        out_path = f'{outdir}/geo_igram_{burst_id}_{ref_date}_{sec_date}.tiff'
    else:
        out_path = f'{outdir}/radar_igram_{burst_id}_{ref_date}_{sec_date}.tiff'

    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(out_path, width, length, 1, gdal.GDT_CFloat32)
    out_ds.SetGeoTransform(ref_ds.GetGeoTransform())
    out_ds.SetProjection(ref_ds.GetProjection())
    out_ds.GetRasterBand(1).WriteArray(igram)
    out_ds.FlushCache()

    return out_path


def geocode_radar_igram(burst, metadata, radar_igram_path,
                        dem_path, out_dir):
    geo2rdr_cfg = metadata['runconfig']['processing']['geo2rdr']
    geo_cfg = metadata['runconfig']['processing']['geocoding']

    # Instantiate geocode object and initialize it
    dem_raster = isce3.io.Raster(dem_path)
    epsg = dem_raster.get_epsg()
    proj = isce3.core.make_projection(epsg)
    ellipsoid = proj.ellipsoid
    
    # Generate geogrid
    geo_grid = generate_geogrids([burst], geo_cfg, dem_path)
    geo_obj = isce3.geocode.GeocodeCFloat32()
    geo_obj.orbit = burst.orbit
    geo_obj.ellipsoid = ellipsoid
    geo_obj.doppler = isce3.core.LUT2d()
    geo_obj.threshold_geo2rdr = geo2rdr_cfg['threshold']
    geo_obj.numiter_geo2rdr = geo2rdr_cfg['numiter']
    geo_obj.lines_per_block = geo2rdr_cfg['lines_per_block']
    geo_obj.data_interpolator = "SINC"
    geo_grid = geo_grid[burst.burst_id]
    geo_obj.geogrid(geo_grid.start_x,
                    geo_grid.start_y,
                    geo_grid.spacing_x,
                    geo_grid.spacing_y,
                    geo_grid.width,
                    geo_grid.length,
                    geo_grid.epsg)

    # Allocate input/output rasters
    rdr_igram_raster = isce3.io.Raster(radar_igram_path)
    out_path = f'{out_dir}/geo_{os.path.basename(radar_igram_path)}'
    geo_igram_raster = isce3.io.Raster(out_path, geo_grid.width,
                                       geo_grid.length, 1, gdal.GDT_CFloat32,
                                       'GTiff')
    geo_obj.geocode(radar_grid=burst.as_isce3_radargrid(), input_raster=rdr_igram_raster,
                    output_raster=geo_igram_raster, dem_raster=dem_raster,
                    output_mode=isce3.geocode.GeocodeOutputMode.INTERP)

    geotransform = [geo_grid.start_x, geo_grid.spacing_x, 0,
                    geo_grid.start_y, 0, geo_grid.spacing_y]
    geo_igram_raster.set_geotransform(geotransform)
    geo_igram_raster.set_epsg(geo_cfg['output_epsg'])
    del geo_igram_raster


def geocode_layover_shadow(burst, metadata, mask_path,
                           dem_path, out_dir):
    geo2rdr_cfg = metadata['runconfig']['processing']['geo2rdr']
    geo_cfg = metadata['runconfig']['processing']['geocoding']

    # Instantiate geocode object and initialize it
    dem_raster = isce3.io.Raster(dem_path)
    epsg = dem_raster.get_epsg()
    proj = isce3.core.make_projection(epsg)
    ellipsoid = proj.ellipsoid
    
    geo_grid = generate_geogrids([burst], geo_cfg, dem_path) 
    geo_grid = geo_grid[burst.burst_id]   
    geo_obj = isce3.geocode.GeocodeFloat32()
    geo_obj.orbit = burst.orbit
    geo_obj.ellipsoid = ellipsoid
    geo_obj.doppler = isce3.core.LUT2d()
    geo_obj.threshold_geo2rdr = geo2rdr_cfg['threshold']
    geo_obj.numiter_geo2rdr = geo2rdr_cfg['numiter']
    geo_obj.lines_per_block = geo2rdr_cfg['lines_per_block']
    geo_obj.data_interpolator = "NEAREST"
    geo_obj.geogrid(geo_grid.start_x,
                    geo_grid.start_y,
                    geo_grid.spacing_x,
                    geo_grid.spacing_y,
                    geo_grid.width,
                    geo_grid.length,
                    geo_grid.epsg)

    mask_raster = isce3.io.Raster(mask_path)
    path = os.path.basename(mask_path).split('.rdr')[0]
    out_path = f'{out_dir}/geo_{path}_{burst.burst_id}.tiff'
    geo_mask_raster = isce3.io.Raster(out_path, geo_grid.width,
                                      geo_grid.length, 1, gdal.GDT_Byte,
                                      'GTiff')
    geo_obj.geocode(radar_grid=burst.as_isce3_radargrid(),
                    input_raster=mask_raster,
                    output_raster=geo_mask_raster, dem_raster=dem_raster,
                    output_mode=isce3.geocode.GeocodeOutputMode.INTERP)

    geotransform = [geo_grid.start_x, geo_grid.spacing_x, 0,
                    geo_grid.start_y, 0, geo_grid.spacing_y]
    geo_mask_raster.set_geotransform(geotransform)
    geo_mask_raster.set_epsg(geo_cfg['output_epsg'])
    del geo_mask_raster


if __name__ == '__main__':
    # Command line parser
    cmd = command_line_parser()

    # Get paths for reference and secondary geocoded/radar
    geo_ref_path = f'{geo_path}/{cmd.ref_date}/{cmd.burst_id}_{cmd.ref_date}_{cmd.pol}.slc'
    geo_sec_path = f'{geo_path}/{cmd.sec_date}/{cmd.burst_id}_{cmd.sec_date}_{cmd.pol}.slc'
    rdr_ref_path = f'{radar_path}/{cmd.burst_id}/{cmd.ref_date}/{cmd.burst_id}_{cmd.ref_date}_{cmd.pol}.slc'
    rdr_sec_path = f'{radar_path}/{cmd.burst_id}/{cmd.sec_date}/{cmd.burst_id}_{cmd.sec_date}_{cmd.pol}.slc'
    mask_path = f'{radar_path}/{cmd.burst_id}/{cmd.ref_date}/layoverShadowMask.rdr'

    # Get reference SAFE file, polarization, and subswath
    ref_safe_file = glob.glob(f'{data_path}/S*SLC*{cmd.ref_date}*zip')[0]
    pol = cmd.pol
    i_subswath = int(cmd.burst_id.split('iw')[1].split('_')[0])

    # Create and save geocoded interferogram
    geo_igram_path = create_igram(geo_ref_path, geo_sec_path,
                                  cmd.ref_date, cmd.sec_date,
                                  cmd.burst_id, cmd.outdir)
    radar_igram_path = create_igram(rdr_ref_path, rdr_sec_path,
                                    cmd.ref_date, cmd.sec_date,
                                    cmd.burst_id, cmd.outdir,
                                    is_geocoded=False)

    # Open metadata for geocoded reference file
    with open(geo_ref_path.replace('slc', 'json'), 'r') as f:
        ref_metadata = json.load(f)

    # Extract reference burst
    orbit_files = sorted(glob.glob(f'{orbit_dir}/*EOF'))
    orbit_path = get_orbit_file_from_list(ref_safe_file, orbit_files)
    bursts = load_bursts(ref_safe_file, orbit_path, i_subswath, pol)
    burst = [dummy for dummy in bursts if dummy.burst_id == cmd.burst_id][0]
    print('burst ID', burst.burst_id, cmd.burst_id)

    # Geocode radar igram using geocode CSLC metadata
    dem_raster = isce3.io.Raster(dem_path)
    ds = gdal.Open(geo_igram_path, gdal.GA_ReadOnly)
    width = ds.RasterXSize
    length = ds.RasterYSize
    ds = None
    geocode_radar_igram(burst, ref_metadata, radar_igram_path,
                        dem_path, cmd.outdir)

    # Geocode layover Shadow mask
    geocode_layover_shadow(burst, ref_metadata, mask_path,
                           dem_path, cmd.outdir)
