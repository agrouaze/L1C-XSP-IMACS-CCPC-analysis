__author__ = 'Lisa Maillard'
__email__ = 'lisa.maillar@ifremer.fr'
__date__ = '2023-09'

import xarray as xr
import numpy  as np
import re 
from datetime import datetime
from shapely  import geometry
#For Figures
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature
import matplotlib as mpl 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
import numpy.ma as ma
import geopandas as gpd

#Lisa's tools
from l1canalysis.map_SAR_variables_toolkit.plot_tools.maps import plot_map_with_data

def preprocess_sar(ds0,varn=['sigma0_filt','normalized_variance_filt','macs_Re','macs_Im'],dims_sel={'channels_used':'VV,VH'}):
    """
    preprocess of SARWAVE L1b NetCDF
    
    	Parameters: 
		ds0 (xarray.DataSet): initial dataset of SAR
    	Return: 
		ds  (xarray.DataSet): preprocessed dataset of SAR

    """
    if ('freq_sample' in ds0.dims) and ('freq_line' in ds0.dims): ds = ds0.drop_dims(['freq_sample','freq_line'])
    else: ds = ds0.copy()

    # adding tile line and tile sample as coordinates

    if (not 'tile_line'   in ds.coords): ds = ds.assign_coords({'tile_line'  :ds.tile_line.values})
    if (not 'tile_sample' in ds.coords): ds = ds.assign_coords({'tile_sample':ds.tile_sample.values})

    # selection of wanted dimensions
    for dims in dims_sel.keys():
        if dims in ds.dims  : ds = ds.sel({dims:dims_sel[dims]})
        if dims in ds.coords: ds = ds.drop(dims)
    # keeps only needed data variables
    kept_vars = []
    for var in varn:
        if var in ds.data_vars: kept_vars+=[var]
    if 'corner_latitude'  in ds.data_vars: kept_vars = kept_vars+['corner_latitude']
    if 'corner_longitude' in ds.data_vars: kept_vars = kept_vars+['corner_longitude']
    if 'land_flag'        in ds.data_vars: kept_vars = kept_vars+['land_flag']
    if len(kept_vars): ds = ds[kept_vars]
    else             : ds = ds[[]]
    # get footrpint information from global attributes (extrait tous les nombres (positifs ou négatifs, entiers ou décimaux) de la chaîne ds0.footprint et les stocke sous forme de valeurs flottantes dans la liste fprt)
    fprt = [float(s) for s in re.findall(r'-?\d+\.?\d*', ds0.footprint)]
    # create a DataArray of footprint information
    da = xr.DataArray( fprt, dims=("fprt"),coords={"fprt": range(len(fprt))})
    # add footprint into dataset 
    ds = ds.assign(footprt =  da)
    # get start and stop dates of each tiles (second precision) and add into dataset as datetime
    dt_start = xr.DataArray( datetime.strptime(ds.start_date[:19], '%Y-%m-%d %H:%M:%S'))
    dt_stop  = xr.DataArray( datetime.strptime(ds.stop_date [:19], '%Y-%m-%d %H:%M:%S'))
    ds = ds.assign(start =  dt_start)
    ds = ds.assign(stop  =  dt_stop)
    return ds

def create_polygon(ds_sar,idx_sw,idx_l,idx_s):
    """
    Creates shapely polygons for a specific SAR tile 
    
    	Parameters: 
		ds_sar (xarray.DataSet): dataset containing SARWAVE L2 informations 
		idx_sw (int): swath index of tile
		idx_l  (int): line index of tile
		idx_s  (int): sample index of tile

    	Return: 
		geometry.Polygon: shapely polygon of the SAR tile (if no values of positions of the tile, returns NaN)

    """
    # if swath dimension exists
    if idx_sw is not None: 
        if not np.any(np.isnan(ds_sar.corner_longitude.isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).values)):
            return geometry.Polygon ([(ds_sar.corner_longitude.isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=0),
              ds_sar.corner_latitude .isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=0)),
             (ds_sar.corner_longitude.isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=0),
              ds_sar.corner_latitude .isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=0)),
             (ds_sar.corner_longitude.isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=1),
              ds_sar.corner_latitude .isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=1)),
             (ds_sar.corner_longitude.isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=1),
              ds_sar.corner_latitude .isel(swath=idx_sw,tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=1))])
        else: 
           return geometry.Polygon ([])
    else: 
        if not np.any(np.isnan(ds_sar.corner_longitude.isel(tile_line=idx_l,tile_sample=idx_s).values)):
            return geometry.Polygon ([(ds_sar.corner_longitude.isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=0),
              ds_sar.corner_latitude .isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=0)),
             (ds_sar.corner_longitude.isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=0),
              ds_sar.corner_latitude .isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=0)),
             (ds_sar.corner_longitude.isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=1),
              ds_sar.corner_latitude .isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=1,c_line=1)),
             (ds_sar.corner_longitude.isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=1),
              ds_sar.corner_latitude .isel(tile_line=idx_l,tile_sample=idx_s).isel(c_sample=0,c_line=1))])
        else: 
           return geometry.Polygon ([])

def get_swath_tiles_polygons_from_l2(ds_sar):
    """
    Creates geopandas geoseries of polygons of all SAR tiles.
    
        Parameters: 
                ds_sar (xarray.DataSet): dataset containing SARWAVE L2 informations
        Return: 
                polygons (geopandas.Geoseries): geoseries of tiles as shapely polygons

    """

    polygons_tmp = []
    for idx_sw in range(len(ds_sar.swath)):
        #print('processing swath number %d' %idx_sw)

        for idx_l  in range(len(ds_sar.tile_line  )):
            for idx_s  in range(len(ds_sar.tile_sample)):
                polygons_tmp.append(create_polygon(ds_sar,idx_sw,idx_l,idx_s))
    
    polygons = gpd.GeoSeries(polygons_tmp,crs=4326)
    return polygons



def get_footprint_polygons_from_l2(ds_sar):
    """
    Creates geopandas geoseries of polygons of all SAR footprints
    
        Parameters: 
                ds_sar (xarray.DataSet): dataset containing SARWAVE L2 informations
        Return: 
                polygons_fp (geopandas.Geoseries): geoseries of footprints as shapely polygons

    """

    polygons_fp = gpd.GeoSeries([geometry.Polygon
                   ([(ds_sar.footprt.isel(swath=idx_sw)[0],ds_sar.footprt.isel(swath=idx_sw)[1]),
                     (ds_sar.footprt.isel(swath=idx_sw)[2],ds_sar.footprt.isel(swath=idx_sw)[3]),
                     (ds_sar.footprt.isel(swath=idx_sw)[4],ds_sar.footprt.isel(swath=idx_sw)[5]),
                     (ds_sar.footprt.isel(swath=idx_sw)[6],ds_sar.footprt.isel(swath=idx_sw)[7]),
                     (ds_sar.footprt.isel(swath=idx_sw)[8],ds_sar.footprt.isel(swath=idx_sw)[9])])
                    for idx_sw in range(len(ds_sar.swath))],crs=4326)
    return polygons_fp


def extract_wanted_sar_values(ds,name,idx_sw=None,idx_l=None,idx_s=None):
    """
    Extracts a data variable from the SAR Dataset for a specific tile.   
    
        Parameters: 
                ds     (xarray.DataSet): dataset containing SARWAVE L2 informations
                name   (str): name of the wanted data variable
                idx_sw (int): swath index of tile
                idx_l  (int): line index of tile
                idx_s  (int): sample index of tile

        Return: 
                (xarray.Dataset): data variable at the specific SAR tile

    """
    # if swath dimension exists
    if len(ds.swath)>1:
        return np.squeeze([ds[name].isel(swath       = idx_sw[idx],
                                         tile_line   = idx_l [idx],
                                         tile_sample = idx_s [idx]).values for idx in range(len(idx_sw))])
    else:
        return np.squeeze(ds[name].isel(tile_line   = idx_l ,
                                         tile_sample = idx_s).values)

def create_gdfs_for_sar(ds_sar,valn,polygons,**kwargs):
    """
    create geodataframes of SAR polygons and associated values and properties. 
    if polygons of footprint (polygons_fp) are in kwargs, add them in the geodataframe
    if values of land (land_mask) and/or quality flag (qflag_mask) are given in kwargs, geodataframe of masked polygons are added, respectively in white and black.

    e.g., 
          ds = xr.open_dataset(filename_sar)
          polygons = get_swath_tiles_polygons_from_l2_local(ds)
          polygons_fp = get_footprint_polygons_from_l2_local(ds)
          gdfs,gdfs_ppties = create_gdfs_for_sar(ds,'sar_hs',polygons,polygons_fp=polygons_fp,land_mask=ds.land_flag,qflag_mask=ds.sar_radar_parameters_flag)
    """    
    gdfs,gdfs_ppties = [],[]

    if 'land_mask' in kwargs:
        polygons_plot = polygons[kwargs['land_mask'].values.flatten()==0]
        val_plot = ds_sar[valn].values.flatten()[kwargs['land_mask'].values.flatten()==0]
    else:
        polygons_plot = polygons
        val_plot = ds_sar[valn].values.flatten()

    # main polygons
    gdf = gpd.GeoDataFrame(geometry=polygons_plot)
    gdf = gdf.assign(color=val_plot)
    gdf_ppties = {'column':'color','zorder':4,'linewidth':0.4}
    gdfs.append(gdf)
    gdfs_ppties.append(gdf_ppties)
    
    # footprint
    if 'polygons_fp' in kwargs:
        gdf = gpd.GeoDataFrame(geometry=kwargs['polygons_fp'])
        gdf_ppties = {'color':'k','alpha':0.2,'edgecolor':'w','linewidth':0.4,'zorder':3}
        gdfs.append(gdf)
        gdfs_ppties.append(gdf_ppties)

    # flagmask
    if 'qflag_mask' in kwargs:
        gdf = gpd.GeoDataFrame(geometry=polygons[kwargs['qflag_mask'].values.flatten()==1])
        if not gdf.empty:
            gdf_ppties = {'color':'k','linewidth':0.4,'zorder':4}
            gdfs.append(gdf)
            gdfs_ppties.append(gdf_ppties)

    # landmask
    if 'land_mask' in kwargs:  
        gdf = gpd.GeoDataFrame(geometry=polygons[kwargs['land_mask'].values.flatten()==1])
        if not gdf.empty:
            gdf_ppties = {'color':'w','linewidth':0.4,'zorder':4,'alpha':0.5}
            gdfs.append(gdf)
            gdfs_ppties.append(gdf_ppties)

    return gdfs,gdfs_ppties



def map_sar(ds,varn,
               display_sar = True,
               pts=None,
               trajs=None,txt_along_traj=None,
               bkgd = None,
               arrow_field=None,
               streamline=None,
               cmap = plt.cm.jet, cmin=None,cmax=None, cval=None,
               percentile_color=None,diverging = False,
               clbl=None,
               lonextent =None,
               latextent = None,
               dlon=3, dlat=3,
               qflag_mask=None,
               land_mask =None,
               title=None,
               ax  = None,
               figsize=(8,8),
               return_fig=False):

    """
    Plots a map (cartopy land) of the SAR tiles and footprints, with other plotting options (add a background, trajectories, points)

        Parameters: 
                ds (xarray.DataSet): dataset containing SARWAVE L2 informations
                varn (str): name of the data variable in ds_sar to be displayed
        Optional Parameters: 
            display_sar (bool)                           : if True, display SAR tiles and footprints
            pts              (list of arrays)            : points (lon,lat,color) or (lon,lat,value) if color must change with value
            trajs            (list of arrays)            : trajectories (i.e., sequence of points) (lon,lat,color) or  (lon,lat,value) if color must change with value
            txt_along_traj   (list of arrays)            : texts to be plotted along trajectories 
            bkgd             (list of arrays)            : background field (lon,lat,values)
            arrow_field      (list of arrays)            : field of arrows (lon,lat,U,V)
            streamline       (list of arrays)            : field of streamlines (lon,lat,U,V,color)
            cmap             (plt.cm)                    : colormap
            cmin             (float)                     : minimum value of colorbar
            cmax             (float)                     : maximum value of colorbar
            cval             (array)                     : values to use for colorbar range calibration
            percentile_color (float)                     : percentile of cval used to set colorbar range (51-100)
            diverging        (boolean)                   : if True, cmin is set to  -cmax
            clbl             (string)                    : colorbar label
            lonextent        (array)                     : longitudinal extent  (degree East)
            latextent        (array)                     : latitudinal extent
            title            (string)                    : title of map
            ax               (axes)                      : axe on which the map is plotted
            figsize          (array)                     : figure size
            return_fig       (bool)                      : if True, return the axis
        Return: 
            ax               (axes)                      : axe on which the map is plotted
            

    """
    polygons    = get_swath_tiles_polygons_from_l2(ds)
    if 'footprt' in ds.data_vars:
        polygons_fp = get_footprint_polygons_from_l2(ds)
        geom_args = {'polygons_fp':polygons_fp}
    else: 
        geom_args = {}
    if qflag_mask is not None: geom_args['qflag_mask'] = qflag_mask
    if land_mask  is not None: geom_args['land_mask' ] = land_mask
    gdfs,gdfs_ppties = create_gdfs_for_sar(ds,varn,polygons,**geom_args)

    #Map coordinates extent
    if latextent is None: latextent = [np.nanmin(ds.latitude)  - dlat,np.nanmax(ds.latitude)  + dlat]
    if lonextent is None: lonextent = [np.nanmin(ds.longitude) - dlon,np.nanmax(ds.longitude) + dlon]

    plotting_args = {'trajs'           :trajs,
                     'pts'             :pts,
                     'txt_along_traj'  :txt_along_traj,
                     'bkgd'            :bkgd,
                     'arrow_field'     :arrow_field,
                     'streamline'      :streamline,
                     'cmap'            :cmap,
                     'cmin'            :cmin,
                     'cmax'            :cmax,
                     'percentile_color':percentile_color,
                     'diverging'       :diverging,
                     'clbl'            :clbl,
                     'lonextent'       :lonextent,
                     'latextent'       :latextent,
                     'ax'              :ax,
                     'title'           :title,
                     'figsize'         :figsize,
                     'return_fig'      :True}
    
    if display_sar: ax = plot_map_with_data(gdfs=gdfs,
                                            gdfs_ppties=gdfs_ppties,cval=ds[varn], **plotting_args)
    else          : ax = plot_map_with_data(**plotting_args)

    # end    
    if return_fig:
        return ax
    else:
        return


