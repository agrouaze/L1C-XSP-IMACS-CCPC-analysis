import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os


def plot_one_xspec(l1b_path, tile_line, tile_sample, group='intraburst', output_dir='.', png_name='xspectra.png'):
    """
    Do two plots af one Cross-spectra (imaginary parts and real parts) in the provided L1B product tile.
    Return the two plots describe above.
    
    Args:
        l1b_path (str) : absolute path to L1B product
        group (str) : 'intraburst' or 'interburst'
        tile_line : position of the tile in the azimuth direction
        tile_sample : position of the tile in the range direction
        the first tile is located on the upper-right corner of the sub-swath & the last tile is located on the lower-left corner of the sub-swath   
    """
    import datatree
    from xsarslc.processing.xspectra import symmetrize_xspectrum
    from xsarslc.tools import xndindex
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

    l1b = datatree.open_datatree(l1b_path)  # open the l1b product as a datatree
    if group == 'intraburst':
        xss = l1b[group]['xspectra_2tau_Re']+1j*l1b[group]['xspectra_2tau_Im']
    elif group == 'interburst':
        xss = l1b[group]['xspectra_Re']+1j*l1b[group]['xspectra_Im']
    else:
        raise ValueError('Unknown group: {}'.format(group))
    xss = xss.assign_coords({'k_rg':xss['k_rg'].mean(dim=set(xss['k_rg'].dims)-set(['freq_sample']), keep_attrs=True)}).swap_dims({'freq_sample':'k_rg', 'freq_line':'k_az'})   # simplify the notation on the data set and set a new dimension, k_rg, wich is the mean of k_rg over all the tiles.
    if '2tau' in xss.dims:
        xss = xss.squeeze(dim='2tau')
    xss = symmetrize_xspectrum(xss)
    xs = xss.sel(tile_line=tile_line, tile_sample=tile_sample)   # select the tile desired
    
    heading = np.radians(l1b[group]['ground_heading'].sel(tile_line=tile_line, tile_sample=tile_sample).data)
    incidence = np.round(l1b[group]['incidence'].sel(tile_line=tile_line, tile_sample=tile_sample).item(),2)
    tau = np.round(l1b[group]['tau'].sel(tile_line=tile_line, tile_sample=tile_sample).item(),3)
    cutoff = np.rint(l1b[group].ds.sel(tile_line=tile_line, tile_sample=tile_sample)['azimuth_cutoff'].data)
    cutoff = int(cutoff) if np.isfinite(cutoff) else cutoff # handle nan cutoff
    lon = np.round(l1b[group].ds.sel(tile_line=tile_line, tile_sample=tile_sample)['longitude'].item(), 2)
    lat = np.round(l1b[group].ds.sel(tile_line=tile_line, tile_sample=tile_sample)['latitude'].item(), 2)
    heading = np.degrees(np.arctan2(np.sin(heading), np.cos(heading)))

    # here are where the plots are made :
    figreal, figimag = xs_figures(xs, heading = heading, incidence = incidence, tau = tau, cutoff = cutoff, lon=lon, lat=lat)  

    
    filename, ext = os.path.splitext(png_name)   # split png_name as follow : ('wallpaper_xspectra', '.png')
    grp = os.path.splitext(group)[0][:-5]        # 'intraburst' became 'intra'
    real_path = os.path.join(output_dir, filename+'_real'+'_'+grp+'_tline'+str(tile_line)+'_tsampl'+str(tile_sample)
                                         +'_lon'+str(lon)+'_lat'+str(lat)+ext)
    imag_path = os.path.join(output_dir, filename+'_imag'+'_'+grp+'_tline'+str(tile_line)+'_tsampl'+str(tile_sample)
                                         +'_lon'+str(lon)+'_lat'+str(lat)+ext)
    figreal.savefig(real_path, dpi='figure',bbox_inches='tight')
    figimag.savefig(imag_path, dpi='figure',bbox_inches='tight')
    

def xs_figures(xs, heading = 0, incidence = None, tau = None, cutoff = None, lon=None, lat=None, kmax = 0.07, kmin = 2*np.pi/1000., **kwargs):
    """
    Return two figures of respectively real and imaginary part of Cross-spectrum. Use as many information as provided

    Args:
        xs (xarray.DataArray): xarray of one cross specrtum with k_rg and k_az coordinate
        heading (float) : heading of the satellite
        incidence (float) : incidence angle [deg] at cross-spectrum location
        tau (float) : reference time delay
        cutoff (float) : SAR cut-off value [m]
        lon (float) : longitude [deg] at cross-spectrum location
        lat (float) : latitude [deg] at cross-spectrum location
        kmax (float) : maximum wavenumber to be plotted
        kmin (float) : wavelength larger than kmin are not plotted
    Return:
        (figure, figure) : figures of real and imaginary cross-spectrum

    """

    from matplotlib import colors as mcolors
    cmap = mcolors.LinearSegmentedColormap.from_list("", ["white","violet","mediumpurple","cyan","springgreen","yellow","red"])
    PuOr = mcolors.LinearSegmentedColormap.from_list("", ["darkgoldenrod","white","purple"])
    
    xs = xs.where((xs['k_rg']**2+xs['k_az']**2)>kmin**2, 0j)
    heading = np.radians(heading)
    keast = xs['k_rg']*np.cos(heading)+xs['k_az']*np.sin(heading)
    knorth = xs['k_az']*np.cos(heading)-xs['k_rg']*np.sin(heading)
    keast.attrs.update({'long_name':'wavenumber in East direction', 'units':'rad/m'})
    knorth.attrs.update({'long_name':'wavenumber in North direction', 'units':'rad/m'})
    xs = xs.assign_coords({'k_east':keast,'k_north':knorth})
    heading = np.arctan2(np.sin(heading), np.cos(heading))
    range_rotation = -np.degrees(heading) if np.abs(np.degrees(heading))<=90 else -np.degrees(heading)+180
    azimuth_rotation = -np.degrees(heading)+90 if np.degrees(heading)>=0 else -np.degrees(heading)-90

    # -----------------------------------------
    
    figreal = plt.figure(figsize=(10,10), tight_layout=True)
    plot = xs.real.plot(cmap=cmap, vmin=0, x='k_east', y='k_north', add_colorbar=False)
    ax = plt.gca()
    for r in [400,200,100]:   # for r in [400,200,100,50]: with busrt dim, here are made the circles
        circle = plt.Circle((0, 0), 2*np.pi/r, color='k', fill=False, linestyle='--', linewidth=0.5)
        ax.add_patch(circle)
        plt.text(-2*np.pi/r*np.cos(np.radians(90)),2*np.pi/r*np.sin(np.radians(90))+0.002,'{} m'.format(r), rotation=0.,horizontalalignment='left',verticalalignment='center')
        # plt.text(-np.sqrt(2)*np.pi/r-0.002,np.sqrt(2)*np.pi/r+0.002,'{} m'.format(r), rotation=45.,horizontalalignment='center',verticalalignment='center')
    for a in [-60,-30,30,60]:
        plt.plot([-0.2*np.cos(np.radians(a)), 0.2*np.cos(np.radians(a))],[-0.2*np.sin(np.radians(a)), 0.2*np.sin(np.radians(a))], color='k', linestyle='--', linewidth=0.5)
    plt.vlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)
    plt.hlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)
    xp = kmax*np.cos(heading)
    yp = kmax*np.sin(heading)
    plt.plot([-xp,xp], [yp,-yp], color='r') # range line
    plt.plot([yp,-yp], [xp,-xp], color='r') # azimuth line
    if cutoff:
        plt.plot(np.array([-xp,xp])+2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])+2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff upper line
        plt.plot(np.array([-xp,xp])-2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])-2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff lower line
    plt.axis('scaled')
    plt.xlim([-kmax,kmax])
    plt.ylim([-kmax,kmax])
    plt.text(0.9*xp,-0.9*yp+0.006,'range', color='r', rotation=range_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # range name
    plt.text(0.85*yp+0.006,0.85*xp,'azimuth', color='r', rotation=azimuth_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # azimuth name
    if cutoff and np.isfinite(cutoff):
        plt.text(-0.96*kmax,-0.96*kmax,'cutoff : {} m'.format(cutoff))
    if incidence and np.isfinite(incidence):
        plt.text(-0.96*kmax,0.93*kmax,'incidence : {} deg'.format(incidence))
    if lon and np.isfinite(lon):
        plt.text(0.46*kmax,0.93*kmax,'longitude : {} deg'.format(lon))
    if lat and np.isfinite(lat):
        plt.text(0.46*kmax,0.86*kmax,'latitude : {} deg'.format(lat))
    if tau and np.isfinite(tau):
        plt.text(-0.96*kmax,0.86*kmax,'tau : {} s'.format(tau))
    cbar = plt.colorbar(plot, ax=ax, label='Real cross-spectra intensity (u.a)', extend='min', aspect=22, shrink=0.73) # colorbar
    plt.title('X-spectrum (real)', fontsize=17)
    plt.show()

    # -----------------------------------------

    figimag = plt.figure(figsize=(10,10), tight_layout=True)
    plot = xs.imag.plot(cmap=PuOr, x='k_east', y='k_north', add_colorbar=False)
    ax = plt.gca()
    # for r in [400,200,100,50]:
    for r in [400,200,100]:
        circle = plt.Circle((0, 0), 2*np.pi/r, color='k', fill=False, linestyle='--', linewidth=0.5)
        ax.add_patch(circle)
        plt.text(-np.sqrt(2)*np.pi/r-0.002,np.sqrt(2)*np.pi/r+0.002,'{} m'.format(r), rotation=45.,horizontalalignment='center',verticalalignment='center')
    for a in [-60,-30,30,60]:
        plt.plot([-0.2*np.cos(np.radians(a)), 0.2*np.cos(np.radians(a))],[-0.2*np.sin(np.radians(a)), 0.2*np.sin(np.radians(a))], color='k', linestyle='--', linewidth=0.5)
    plt.vlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)
    plt.hlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)

    plt.plot([-xp,xp], [yp,-yp], color='r') # range line
    plt.plot([yp,-yp], [xp,-xp], color='r') # azimuth line
    if cutoff:
        plt.plot(np.array([-xp,xp])+2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])+2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff upper line
        plt.plot(np.array([-xp,xp])-2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])-2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff lower line

    plt.axis('scaled')
    plt.xlim([-kmax,kmax])
    plt.ylim([-kmax,kmax])
    # plt.grid()
    plt.text(0.9*xp,-0.9*yp+0.006,'range', color='r', rotation=range_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # range name
    plt.text(0.85*yp+0.006,0.85*xp,'azimuth', color='r', rotation=azimuth_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # azimuth name
    if cutoff and np.isfinite(cutoff):
        plt.text(-0.96*kmax,-0.96*kmax,'cutoff : {} m'.format(cutoff))
    if incidence and np.isfinite(incidence):
        plt.text(-0.96*kmax,0.93*kmax,'incidence : {} deg'.format(incidence))
    if lon and np.isfinite(lon):
        plt.text(0.46*kmax,0.93*kmax,'longitude : {} deg'.format(lon))
    if lat and np.isfinite(lat):
        plt.text(0.46*kmax,0.86*kmax,'latitude : {} deg'.format(lat))
    if tau and np.isfinite(tau):
        plt.text(-0.96*kmax,0.86*kmax,'tau : {} s'.format(tau))
    cbar = plt.colorbar(plot, ax=ax, label='Imaginary cross-spectra intensity (u.a)', aspect=22, shrink=0.73) # colorbar
    plt.title('X-spectrum (imaginary)', fontsize=17)
    plt.show()

    return figreal, figimag