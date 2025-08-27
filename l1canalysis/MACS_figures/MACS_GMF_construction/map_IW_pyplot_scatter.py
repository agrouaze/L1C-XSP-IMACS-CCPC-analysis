"""
A Grouazel
March 2025
working for IW and WV (new version 1 nc = 1 SAFE)
"""
from matplotlib import pyplot as plt
import cartopy
import os

os.environ["CARTOPY_DATA_DIR"] = "/home1/datahome/agrouaze/.local/share/cartopy/"
# cartopy.config["pre_existing_data_dir"] = "/home1/datahome/agrouaze/.local/share/cartopy/shapefiles/natural_earth/physical/"
cartopy.config["pre_existing_data_dir"] = "/home1/datahome/agrouaze/.local/share/cartopy/shapefiles/natural_earth/physical/"
# cartopy.config["downloaders"] = None
print(cartopy.config["data_dir"])
import cartopy.crs as ccrs
variables_real_names = {'phase_ma_crosssp':"Phase of MACS"}
import numpy as np


from matplotlib.colors import LinearSegmentedColormap

# Define the colors: green for the low end, white in the middle, and purple for the high end
colors = [
    (0.0, "green"),  # Start color (low value)
    (0.5, "cyan"),  # Middle color (neutral value)
    (1.0, "yellow")  # End color (high value)
]

# Create the custom colormap
custom_cmap = LinearSegmentedColormap.from_list("CustomCmap", colors)


# Define the colors: green for the low end, white in the middle, and purple for the high end
colors = [
    (0.0, "blue"),  # Start color (low value)
    (0.5, "grey"),  # Middle color (neutral value)
    (1.0, "orange")  # End color (high value)
]

# Create the custom colormap
custom_cmap_im = LinearSegmentedColormap.from_list("CustomCmap", colors)
PuOr = LinearSegmentedColormap.from_list("", ["darkgoldenrod","white","purple"])
variables_real_names = {'phase_ma_crosssp':"Phase of MACS"}
clims_intra = {'phase_ma_crosssp':None,# (-3,3) # (-30,30)
             'doppler_centroid':(-3,3),
             'macs_Re':(-0.8,0.8),
             'macs_Im':(-0.3,0.3), #(-0.5,0.5)
               'IMACS_GMF':(-0.3,0.3), #(-0.5,0.5)
               'IMACS_anomaly':(-0.1,0.1),
             'abs_macs':(0,0.8),
            'sigma0_filt': (0,0.3),
            'sigma0': (0,0.3),
               'windspeed': (0,30),
            'incidence': (22,47),
               'ground_heading':(0,360),
            'azimuth_cutoff':(100,400),
               'CCCP_Re':(-0.08,0.08),
            'CCCP_Im':(-0.1,0.1),
            }
clims_inter = {'phase_ma_crosssp':None, # (-3,3)
             'doppler_centroid':(-3,3),
             'macs_Re':(-0.02,0.02),
             'macs_Im':(-0.01,0.01),
            'IMACS_GMF':(-0.01,0.01), #(-0.5,0.5)
               'IMACS_anomaly':(-0.01,0.01),
             'abs_macs':(0,0.05),
            'sigma0_filt': (0,0.1),
            'sigma0': (0,0.1),
            'incidence': (22,47),
            'ground_heading':(0,360),
               'azimuth_cutoff': (100, 400),
               'CCCP_Re': (-0.08, 0.08),
               'CCCP_Im': (-0.1, 0.1),
            }
colormaps = {'phase_ma_crosssp':'jet',
             'doppler_centroid':'bwr',
             'macs_Re':custom_cmap,
             'macs_Im':'bwr',
             'IMACS_GMF':'bwr',
             'IMACS_anomaly':'bwr',
             'abs_macs':'viridis',
             'sigma0_filt':'Greys_r',
            'sigma0':'Greys_r',
             'windspeed':'jet',
             'incidence':'magma',
            'ground_heading':'plasma',
             'azimuth_cutoff': 'jet',
             'CCCP_Re': PuOr,
             'CCCP_Im': PuOr,
            }

def plot_maps_gulf(ax,variables_vect,variable,l1c_lon,l1c_lat,burst_bounds,grp,lon_coarse=None,lat_coarse=None,coarse_param=None):
    if variable in variables_real_names:
        vava = variables_real_names[variable]
    else:
        vava = variable
    if grp=='intraburst':
        clims = clims_intra
    else:
        clims = clims_inter
    plt.title(vava)
    # ax.set_extent([-179.99999, 180, -80, 80])
    # ax.set_extent([-67,-59,33,46])
    ax.coastlines(antialiased=True)
    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
    # if clims[variable] is not None:
    #     clicli = clims[variable]
    # else:
    #     clicli
    im = plt.scatter(l1c_lon,l1c_lat,c=variables_vect[variable],s=1,clim=clims[variable],cmap=colormaps[variable]) # clim=(-0.01,0.01) clim=(-0.1,0.1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabel_style = {'color': 'gray', 'weight': 'bold','rotation':45}
    # Modify the properties of the contour lines
    levels = np.arange(-50, 50, 10)
    if lon_coarse is not None:
        cc = plt.contour(lon_coarse,lat_coarse,coarse_param,levels=levels,colors='m')
        linewidths = [0.1+abs(l/40) for l in levels]  # Linewidth depends on level
        alphas = [0.001+abs(l)/55 for l in levels]  # Alpha depends on level
        for i, collection in enumerate(cc.collections):
            # print(dir(collection))
            collection.set_linewidth(linewidths[i])
            collection.set_alpha(alphas[i])
    cb = plt.colorbar(im,fraction=0.02, pad=0.04)
    # cb.set_label(variable)
    for rr in burst_bounds:
        # plt.plot(rr[:,0],rr[:,1],'-')
        plt.plot(*rr.exterior.xy,'-',color='black',lw=0.7,alpha=0.7)



def plot_map_tiles(gdf,ax,burst_bounds,variable,grp):
    if grp=='intraburst':
        clims = clims_intra
    else:
        clims = clims_inter
    vmin,vmax = clim=clims[variable]
    gdf.plot(column='value', cmap=colormaps[variable], legend=True, edgecolor='black', ax=ax,vmin=vmin,vmax=vmax) #,
    ax.coastlines(antialiased=True)
    ax.add_feature(cartopy.feature.LAND, zorder=100, edgecolor='k')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.ylabels_right = False
    gl.ylabels_left = False
    gl.xlabel_style = {'color': 'gray', 'weight': 'bold','rotation':45}
    # Modify the properties of the contour lines
    # levels = np.arange(-50, 50, 10)
    # cb = plt.colorbar(im,fraction=0.02, pad=0.04)
    # cb.set_label(variable)
    for rr in burst_bounds:
        # plt.plot(rr[:,0],rr[:,1],'-')
        plt.plot(*rr.exterior.xy,'-',color='black',lw=0.7,alpha=0.7)


