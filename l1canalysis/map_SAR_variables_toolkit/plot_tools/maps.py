__author__ = 'Lisa Maillard'
__email__ = 'lisa.maillar@ifremer.fr'
__date__ = '2023-09'

import numpy as np 
import matplotlib.pyplot as plt
import cartopy.crs       as ccrs
import cartopy.feature   as cfeature
import matplotlib as mpl
from matplotlib import patches
from matplotlib import text as mtext
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes


def plot_map(ax=None,
             figsize=(10,10),
             proj=ccrs.PlateCarree(),
             latextent=[-40,40],
             lonextent=[-180,179.9]
             ):
    """
    Plots a map with land from cartopy. 
 
            Parameters:
                    ax        (axes or None)      : axe on which to plot the figure. If None, create new axe and figure.
                    figsize   (list)              : figure size
                    proj      (cartopy projection): map projection
                    lonextent (array)             : array with minimum and maximum longitude extent (from -180 to 180°E)
                    latextent (array)             : array with minimum and maximum latitude  extent

            Returns:
                    ax        (axes or None)      : axe on which to plot the figure. If None, create new axe and figure.
    """
    ax=init_map(ax,figsize,proj)
    ax=set_map_extent(ax,lonextent,latextent)
    ax=end_map(ax)
    return ax
    

def plot_background_map(ax,bkgd,**kwargs):
    """
    Plot a background to an existing map
            Parameters:
                    ax        (axis)       : axe on which to plot the figure.
                    bkgd      (3d list of arrays) : background information (lon,lat,values). If values has shape 2, use pcolormesh. If values has shape 1, use tripcolor

            Returns:

    """
    if   len(np.shape(bkgd[2])) == 1: ax.tripcolor(bkgd[0],bkgd[1],bkgd[2],**kwargs)
    elif len(np.shape(bkgd[2])) == 2: ax.pcolormesh(bkgd[0],bkgd[1],bkgd[2],**kwargs)
    return

def plot_geodataframe(ax,gdf,**kwargs):
    """
    Plot a geodataframe to an existing axis
            Parameters:
                    ax        (axis)       : axe on which to plot the figure.
                    gdf       (pandas.geodataframe) : Geodataframe containing geometries to plot. To specify colors, add it in kwargs as {'column': X}
                                                                                                    example : 
                                                                                                        gdf = gpd.GeoDataFrame(geometry=polygons)
                                                                                                        gdf = gdf.assign(color=myvals)
                                                                                                        gdf_ppties = {'column':'color'}
                                                                                                        plot_geodataframe(ax,gdf,**gdf_ppties)

            Returns:
    """
    gdf.to_crs(ccrs.PlateCarree()).plot(ax=ax,**kwargs)
    return

def plot_arrowfield(ax,arr,**kwargs):
    """
    Plot a field of arrows to an existing axis
            Parameters:
                    ax        (axis)       : axe on which to plot the figure.
                    arr      (4d list of arrays)  : arrows information (lon,lat,U,V).

            Returns:

    """

    ax.quiver(arr[0],arr[1],arr[2],arr[3],**kwargs)
    return
    
def plot_streamlines(ax,stml,**kwargs):
    """
    Plot streamlines to an existing axis
            Parameters: 
                    ax        (axis)       : axe on which to plot the figure.
                    stml      (5d list of arrays)  : streamlines information (lon,lat,U,V,color).

            Returns:

    """
    ax.streamplot(stml[0],stml[1],stml[2],stml[3],color=stml[4],**kwargs)
    return

def plot_point(ax,pt,**kwargs):
    """
    Plot a unique point to an existing axis
            Parameters: 
                    ax        (axis)       : axe on which to plot the figure.
                    pt       (3d list of arrays)  : point information (lon,lat,color).

            Returns:
    """

    ax.plot   (pt[0],pt[1],c=pt[2],linestyle='',**kwargs)
    return

def plot_trajectory(ax,traj,argscat={'s':8},**kwargs):
    """	
            Parameters: 
                    ax        (axis)       : axe on which to plot the figure.
                    traj      (3d list of arrays)  : trajectory information (lon,lat,colors).
                    argscat   (dictionnary) : properties of the scatter 

            Returns:
    """
    ax.scatter(traj[0],traj[1],c='k',s=17,edgecolors='k',linewidths = 0.05,**kwargs)
    ax.scatter(traj[0],traj[1],c=traj[2],**{**argscat,**kwargs})
    return

def init_map(ax,figsize,proj):
    """
    initialise a cartopy map
            Parameters: 
                    ax        (axis)       : axe on which to plot the figure.
                    figsize   (list)              : figure size
                    proj      (cartopy projection): map projection
            Returns:
                    ax        (axis)       : axe on which to plot the figure.
    """
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1,1,1, projection=proj)
    ax.add_feature(cfeature.LAND,color='lightgrey',alpha=1,zorder=2)
    #ax.add_feature(cfeature.COASTLINE,alpha=1,zorder=2)
    return ax

def set_map_extent(ax,lonextent,latextent):
    """
    set coordinate extent of map
            Parameters: 
                    ax        (axis)       : axe on which to plot the figure.
                    lonextent (array)             : array with minimum and maximum longitude extent (from -180 to 180°E)
                    latextent (array)             : array with minimum and maximum latitude  extent

            Returns:
                    ax        (axis)       : axe on which to plot the figure.
    """
    ax.set_extent([lonextent[0], lonextent[1],
                          latextent[0], latextent[1]], crs=ccrs.PlateCarree())
    return ax 

def end_map(ax):
    """
    finish drawing of gridlines and label styles of a cartopy map
            Parameters: 
                    ax        (axis)              : axe on which to plot the figure.
            Returns:
                    ax        (axis)       : axe on which to plot the figure.
    """

    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1, color='gray', alpha=0.8)
    gl.xlabel_style = {'size': 12,'color': 'k'}
    gl.ylabel_style = {'size': 12,'color': 'k'}
    gl.top_labels = gl.right_labels = False
    return ax

def add_colorbar(ax,cmap,norm,clbl):
    """
    add a colorbar onto a new axis                               
            Parameters: 
                    ax        (axis)              : axe on which to plot the figure.
                    cmap (matplotlib.pyplot.cm) : colormap
                    norm  (mpl.colors.Normalize) : normalisation used (e.g., mpl.colors.Normalize(vmin=0, vmax=10))
                    clbl (str): label of colorbar
            Returns:
                    ax        (axis)       : axe on which to plot the figure.
    """

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05,axes_class=maxes.Axes)
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='proportional')
    # colorbar label
    if clbl is not None: cb.set_label(clbl)
    return ax,cb

def plot_map_hexbin(lonval,latval,val,
                    cmap,norm,crs,
                    lonextent,latextent,figsize=(8,8),gridsize=[400,100],
                    reduce_C_function=np.mean,clabel='my val',fnsave=None):
    """
    to generate hexbin maps with all data falling in hexbin being reduced using the reduce_C_function (e.g., np.mean, np.std,...)
    """
    ax = plot_map(proj=crs,
                  latextent=latextent,
                  lonextent=lonextent,
                  figsize=figsize
                 )
    
    ax.hexbin(lonval,latval,C=val,reduce_C_function=reduce_C_function, transform=crs,cmap=cmap,norm=norm,gridsize=gridsize)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='1%', pad=0.05,axes_class=maxes.Axes)
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='proportional')
    cb.set_label(clabel)
    if fnsave is not None:
        plt.savefig(fnsave,dpi=300)
    return ax

def define_colorbar_ppties(cval=None,cmin=None,cmax=None,cmap=plt.cm.jet,percentile_color=None,diverging=False):
    """
    define properties of colorbar                                
            Parameters: 
                    cval  (array) : values to use for colorbar range calibration
                    cmin  (float) : minimum value of colorbar range
                    cmax  (float) : maximum value of colorbar range
                    cmap (matplotlib.pyplot.cm): colormap
                    percentile_color (float): percentile of cval used to set colorbar range (0-100)
                    diverging (bool) : if True, set cmin = - cmax
            Returns:
                    (dictionnary) : dictionnary containing cmap and norm information
    """

    if cval is not None: 
        if cmax is None:
            if percentile_color is None: cmax = np.nanmax(cval)
            else                       : cmax = float(np.nanpercentile(cval,percentile_color))
        if cmin is None:
            if percentile_color is None: cmin = np.nanmin(cval)
            else                       : cmin = float(np.nanpercentile(cval,100.-percentile_color))
    else: 
        if (cmax is None) or (cmin is None): 
            print('User should either give colorbar ranges (cmin,cmax) or values on which the colorbar will be based (cval)')
            return 

    if (diverging) and (cmin != -cmax):
        if np.abs(cmin)<np.abs(cmax): cmin = -cmax
        else                        : cmax = -cmin
            
    norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

    return {'cmap':cmap,'norm':norm}

def plot_map_with_data(gdfs=None, gdfs_ppties=None,
                  pts=None,
                  trajs=None,txt_along_traj=None,
                  bkgd = None,
                  arrow_field=None,
                  streamline=None,
                  cmap = plt.cm.jet, cmin=None,cmax=None, cval=None,
                  percentile_color=None,diverging = False,
                  clbl=None, 
                  lonextent = [-180,179.9],
                  latextent = [-60,60],
                  title=None,
                  ax  = None,
                  figsize=(8,8),
                  return_fig=False):

    """
    Plot a cartopy map with wanted background,trajectories, geopandas geometries, points, arrow fields, streamlines
        Optional Parameters: 
            gdfs             (list of gpd Geodataframes) : geodataframes containing geometries 
            gdfs_ppties      (list of dict)              : plotting properties of geodataframes
                                                        example : 
                                                            gdf = gpd.GeoDataFrame(geometry=polygons)
                                                            gdf = gdf.assign(color=myvals)
                                                            gdf_ppties = {'column':'color'}
                                                            plot_map_with_data([gdf],[gdf_ppties])
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
            percentile_color (float)                     : percentile of cval used to set colorbar range (0-100)
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
    proj = ccrs.PlateCarree()
    proj_dict = {'transform':proj}
    
    # Start cartopy figure      
    ax = init_map(ax,figsize,proj)
    
    # define colorbar min and max values
    if (cval is not None) or (cmin is not None and cmax is not None): cmap_dict = define_colorbar_ppties(cval,cmin,cmax,cmap,percentile_color,diverging)
    
    # Map in background
    if bkgd is not None: plot_background_map(ax,bkgd,**{**cmap_dict,**proj_dict,**{'zorder':1}})

    # geometries stored in geodataframes
    if gdfs is not None :
        for i,gdf in enumerate(gdfs):
            if 'color' in gdfs_ppties[i]: plot_geodataframe(ax,gdf,**{**gdfs_ppties[i],**proj_dict})
            else                        : plot_geodataframe(ax,gdf,**{**gdfs_ppties[i],**proj_dict,**cmap_dict})

    # arrow field
    if arrow_field is not None:
        plot_arrowfield(ax,arrow_field,**proj_dict)

    # streamlines
    if streamline is not None:
        plot_streamlines(ax,streamline,**{**proj_dict,**cmap_dict})

    # Buoy (or points)
    if pts is not None:
        argplot = {'marker':'x'}
        argscat = {'s':10,'marker':'o','edgecolor':'k'}
        for i,pt in enumerate(pts): 
                if type(pt[2]) == str: plot_point(ax,pt,**{**proj_dict,**{'zorder':6},**argplot})
                else                 : plot_trajectory(ax,pt,{**cmap_dict,**argscat},**{**proj_dict,**{'zorder':6}})

    # Track of trajectory
    if trajs is not None:
        trajs_ppties={'s':8}
        for i,traj in enumerate(trajs):
            if type(traj[2])==str: plot_trajectory(ax,traj,trajs_ppties,**{**proj_dict,**{'zorder':5}})
            else                 : plot_trajectory(ax,traj,{**trajs_ppties,**cmap_dict},**{**proj_dict,**{'zorder':5}},)
            if txt_along_traj is not None: 
                text = CurvedText(x = traj[0],y = traj[1],text=txt_along_traj[i],va = 'bottom',axes = ax)

    # title
    if title is not None: ax.set_title(title)

    # map extent
    ax.set_extent([lonextent[0], lonextent[1],
                   latextent[0], latextent[1]], crs=ccrs.PlateCarree())

    # finish plotting map 
    ax = end_map(ax)

    # Colorbar
    if (cval is not None) or (cmin is not None and cmax is not None): ax,cb = add_colorbar(ax,cmap_dict['cmap'],cmap_dict['norm'],clbl)

    # end    
    if return_fig:
        return ax
    else:
        return


class CurvedText(mtext.Text):
    """
    A text object that follows an arbitrary curve.
    """
    def __init__(self, x, y, text, axes, **kwargs):
        super(CurvedText, self).__init__(x[0],y[0],' ', **kwargs)

        axes.add_artist(self)

        ##saving the curve:
        self.__x = x
        self.__y = y
        self.__zorder = self.get_zorder()

        ##creating the text objects
        self.__Characters = []
        for c in text:
            if c == ' ':
                ##make this an invisible 'a':
                t = mtext.Text(0,0,'a')
                t.set_alpha(0.0)
            else:
                t = mtext.Text(0,0,c, **kwargs)

            #resetting unnecessary arguments
            t.set_ha('center')
            t.set_rotation(0)
            t.set_zorder(self.__zorder +3)

            self.__Characters.append((c,t))
            axes.add_artist(t)


    ##overloading some member functions, to assure correct functionality
    ##on update
    def set_zorder(self, zorder):
        super(CurvedText, self).set_zorder(zorder)
        self.__zorder = self.get_zorder()
        for c,t in self.__Characters:
            t.set_zorder(self.__zorder+3)

    def draw(self, renderer, *args, **kwargs):
        """
        Overload of the Text.draw() function. Do not do
        do any drawing, but update the positions and rotation
        angles of self.__Characters.
        """
        self.update_positions(renderer)

    def update_positions(self,renderer):
        """
        Update positions and rotations of the individual text elements.
        """

        #preparations

        ##determining the aspect ratio:
        ##from https://stackoverflow.com/a/42014041/2454357

        ##data limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        ## Axis size on figure
        figW, figH = self.axes.get_figure().get_size_inches()
        ## Ratio of display units
        _, _, w, h = self.axes.get_position().bounds
        ##final aspect ratio
        aspect = ((figW * w)/(figH * h))*(ylim[1]-ylim[0])/(xlim[1]-xlim[0])

        #points of the curve in figure coordinates:
        x_fig,y_fig = (
            np.array(l) for l in zip(*self.axes.transData.transform([
            (i,j) for i,j in zip(self.__x,self.__y)
            ]))
        )

        #point distances in figure coordinates
        x_fig_dist = (x_fig[1:]-x_fig[:-1])
        y_fig_dist = (y_fig[1:]-y_fig[:-1])
        r_fig_dist = np.sqrt(x_fig_dist**2+y_fig_dist**2)

        #arc length in figure coordinates
        l_fig = np.insert(np.cumsum(r_fig_dist),0,0)

        #angles in figure coordinates
        rads = np.arctan2((y_fig[1:] - y_fig[:-1]),(x_fig[1:] - x_fig[:-1]))
        degs = np.rad2deg(rads)


        rel_pos = 10
        for c,t in self.__Characters:
            #finding the width of c:
            t.set_rotation(0)
            t.set_va('center')
            bbox1  = t.get_window_extent(renderer=renderer)
            w = bbox1.width
            h = bbox1.height

            #ignore all letters that don't fit:
            if rel_pos+w/2 > l_fig[-1]:
                t.set_alpha(0.0)
                rel_pos += w
                continue

            elif c != ' ':
                t.set_alpha(1.0)

            #finding the two data points between which the horizontal
            #center point of the character will be situated
            #left and right indices:
            il = np.where(rel_pos+w/2 >= l_fig)[0][-1]
            ir = np.where(rel_pos+w/2 <= l_fig)[0][0]

            #if we exactly hit a data point:
            if ir == il:
                ir += 1

            #how much of the letter width was needed to find il:
            used = l_fig[il]-rel_pos
            rel_pos = l_fig[il]

            #relative distance between il and ir where the center
            #of the character will be
            fraction = (w/2-used)/r_fig_dist[il]

            ##setting the character position in data coordinates:
            ##interpolate between the two points:
            x = self.__x[il]+fraction*(self.__x[ir]-self.__x[il])
            y = self.__y[il]+fraction*(self.__y[ir]-self.__y[il])

            #getting the offset when setting correct vertical alignment
            #in data coordinates
            t.set_va(self.get_va())
            bbox2  = t.get_window_extent(renderer=renderer)

            bbox1d = self.axes.transData.inverted().transform(bbox1)
            bbox2d = self.axes.transData.inverted().transform(bbox2)
            dr = np.array(bbox2d[0]-bbox1d[0])

            #the rotation/stretch matrix
            rad = rads[il]
            rot_mat = np.array([
                [math.cos(rad), math.sin(rad)*aspect],
                [-math.sin(rad)/aspect, math.cos(rad)]
            ])

            ##computing the offset vector of the rotated character
            drp = np.dot(dr,rot_mat)

            #setting final position and rotation:
            t.set_position(np.array([x,y])+drp)
            t.set_rotation(degs[il])

            t.set_va('center')
            t.set_ha('center')

            #updating rel_pos to right edge of character
            rel_pos += w-used
