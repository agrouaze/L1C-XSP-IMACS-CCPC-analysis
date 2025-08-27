from matplotlib import pyplot as plt
import os
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader

# Define local data path
# local_shapefile_dir = '/home1/datahome/ljessel/Cartopy/'  # Remplacez par le chemin vers vos fichiers décompressés
local_shapefile_dir = '/home/datawork-cersat-public/cache/project/sarwave/tools/landmask/cartopy/shapefiles/natural_earth/physical'
shapefile_path = os.path.join(local_shapefile_dir, 'ne_110m_coastline.shp')

# Read the local data by using Reader
coastline = Reader(shapefile_path).geometries()
coastline_feature = cfeature.ShapelyFeature(coastline, ccrs.PlateCarree())

def map_macs_scatter(df,satellite='S1A+B',part='Re',burstkind='intraburst',polarization='vv',lambda_val='50'):
    """

    :param df:
    :param satellite:
    :param part:
    :param burstkind:
    :param polarization:
    :param lambda_val:
    :return:
    """
    sup18ms_mask = df['Wspeed'] > 18
    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))
    imacs_sup18ms = df[varname][sup18ms_mask]
    lat_sup18ms = df['latitude'][sup18ms_mask];
    lon_sup18ms = df['longitude'][sup18ms_mask]
    # display(lat_sup18ms,lon_sup18ms)

    # Create the map
    fig = plt.figure(figsize=(36, 20))
    ax = plt.axes(projection=ccrs.PlateCarree())
    scat = plt.scatter(lon_sup18ms, lat_sup18ms, s=1, c=imacs_sup18ms.values, cmap='seismic',
                       vmin=-0.02, vmax=0.02)
    cbar = plt.colorbar(scat, label='IMACS %s %sm'%(polarization,lambda_val), extend='both', shrink=0.5, aspect=35)  # colorbar
    ax.add_feature(coastline_feature, facecolor='none', edgecolor='black')
    gl = ax.gridlines(
        draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    plt.title(satellite+' | IMACS where wind speed > 18 m/s | %s'%burstkind, fontsize=25)
    plt.show()