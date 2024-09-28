import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_wallpaper(l1b_path, group='intraburst', output_dir='.', png_name='wallpaper_xspectra.png'):
    """
    Plot the output of make_wallpaper()
    
    Args:
        l1b_path (str) : absolute path to L1B product
        group (str, optional) : 'intraburst' or 'interburst'
        output_dir (str) : output directory of the plot
        png_name (str) : root name of the figures
    """
    filename, ext = os.path.splitext(png_name)
    real_path = os.path.join(output_dir, filename+'_'+group+'_real'+ext)
    imag_path = os.path.join(output_dir, filename+'_'+group+'_imag'+ext)
    xsplot, orbit_pass = make_wallpaper(l1b_path, group = group)
    if orbit_pass.upper()=='ASCENDING':
        xsplot = xsplot.stack(trange=['tile_sample','range'])
        # xsplot = xsplot.assign_coords({'burst':np.flip(xsplot['burst'].data), 'tile_line':np.flip(xsplot['tile_line'].data)}).sortby('burst', 'tile_line')
        xsplot = xsplot.assign_coords({'tile_line':np.flip(xsplot['tile_line'].data)}).sortby('tile_line')
        # xsplot = xsplot.stack(tazimuth=['burst','tile_line','azimuth'])
        xsplot = xsplot.stack(tazimuth=['tile_line','azimuth'])
    elif orbit_pass.upper()=='DESCENDING':
        xsplot = xsplot.assign_coords({'tile_sample':np.flip(xsplot['tile_sample'].data)}).sortby('tile_sample')
        xsplot = xsplot.stack(trange=['tile_sample','range'])
        # xsplot = xsplot.stack(tazimuth=['burst','tile_line','azimuth'])
        xsplot = xsplot.stack(tazimuth=['tile_line','azimuth'])
    else:
        raise ValueError('Unknown orbit pass: {}'.format(orbit_pass))

    plt.figure(figsize=(320,180))
    plt.imshow(xsplot['xspectra_real'].transpose('tazimuth','trange','rgb').data)
    plt.savefig(real_path, dpi='figure',bbox_inches='tight')
    plt.close()
    plt.figure(figsize=(320,180))
    plt.imshow(xsplot['xspectra_imag'].transpose('tazimuth','trange','rgb').data)
    plt.savefig(imag_path, dpi='figure',bbox_inches='tight')
    plt.close()



def make_wallpaper(l1b_path, group):
    """
    Do a plot af all Cross-spectra in the provided L1B product.
    Return a xarray.Dataset with same dimensions as L1B plus three new dimensions ('azimuth','range','rgb') corresponding to
    pixel number in ('azimuth', 'range') direction of the plots and red-green-blue values
    
    Args:
        l1b_path (str) : absolute path to L1B product
        group (str) : 'intraburst' or 'interburst'
    Return:
        xarray.Dataset : dataset of intraburst cross-spectra plots
        str : orbit pass (ascending or descending)
    """
    import datatree
    from xsarslc.processing.xspectra import symmetrize_xspectrum
    from xsarslc.tools import xndindex
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

    l1b = datatree.open_datatree(l1b_path)
    if group == 'intraburst':
        xss = l1b[group]['xspectra_2tau_Re']+1j*l1b[group]['xspectra_2tau_Im']
    elif group == 'interburst':
        xss = l1b[group]['xspectra_Re']+1j*l1b[group]['xspectra_Im']
    else:
        raise ValueError('Unknown group: {}'.format(group))
    xss = xss.assign_coords({'k_rg':xss['k_rg'].mean(dim=set(xss['k_rg'].dims)-set(['freq_sample']), keep_attrs=True)}).swap_dims({'freq_sample':'k_rg', 'freq_line':'k_az'})
    if '2tau' in xss.dims:
        xss = xss.squeeze(dim='2tau')
    xss = symmetrize_xspectrum(xss)
    

    xsplot = list()
    # for mytile in xndindex({'burst':xss.sizes['burst'], 'tile_line':xss.sizes['tile_line'], 'tile_sample':xss.sizes['tile_sample']}):
    for mytile in xndindex({'tile_line':xss.sizes['tile_line'], 'tile_sample':xss.sizes['tile_sample']}):

        xs = xss[mytile]
        heading = np.radians(l1b[group]['ground_heading'][mytile].data)
        incidence = np.round(l1b[group]['incidence'][mytile].item(),2)
        tau = np.round(l1b[group]['tau'][mytile].item(),3)
        cutoff = np.rint(l1b[group].ds[mytile]['azimuth_cutoff'].data)
        cutoff = int(cutoff) if np.isfinite(cutoff) else cutoff # handle nan cutoff
        lon = np.round(l1b[group].ds[mytile]['longitude'].item(), 2)
        lat = np.round(l1b[group].ds[mytile]['latitude'].item(), 2)
        heading = np.degrees(np.arctan2(np.sin(heading), np.cos(heading)))
        figreal, figimag = xs_figures(xs, heading = heading, incidence = incidence, tau = tau, cutoff = cutoff, lon=lon, lat=lat)

        canvas = FigureCanvas(figreal)
        width, height = figreal.get_size_inches() * figreal.get_dpi()
        canvas.draw()
        imagereal = np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(int(height), int(width), 3)[146:854,102:810,:]
        plt.close()

        canvas = FigureCanvas(figimag)
        width, height = figimag.get_size_inches() * figimag.get_dpi()
        canvas.draw()
        imageimag = np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(int(height), int(width), 3)[146:854,102:810,:]
        plt.close()

        xsreal = xr.DataArray(imagereal, dims=('azimuth','range','rgb'), name='xspectra_real').assign_coords(mytile)
        xsimag = xr.DataArray(imageimag, dims=('azimuth','range','rgb'), name='xspectra_imag').assign_coords(mytile)
        xs = xr.merge([xsreal, xsimag])
        xsplot.append(xs)
    # xsplot = xr.combine_by_coords([xs.expand_dims(['burst', 'tile_line','tile_sample']) for xs in xsplot])
    xsplot = xr.combine_by_coords([xs.expand_dims(['tile_line','tile_sample']) for xs in xsplot])
    return xsplot, l1b[group].attrs['orbit_pass']



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
    
    figreal = plt.figure(figsize=(10,10), tight_layout=True)
    xs.real.plot(cmap=cmap, vmin=0, x='k_east', y='k_north')
    ax = plt.gca()
    
    # for r in [400,200,100,50]:
    for r in [400,200,100]:
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
    plt.title('X-spectrum (real)')
    plt.close()

    # -----------------------------------------

    figimag = plt.figure(figsize=(10,10), tight_layout=True)
    xs.imag.plot(cmap=PuOr, x='k_east', y='k_north')
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
    plt.title('X-spectrum (imaginary)')
    plt.close()

    return figreal, figimag
    



    # def make_wallpaper(l1b_path):
    # """
    # """
    # from xsarslc.processing.xspectra import symmetrize_xspectrum
    # from xsarslc.tools import xndindex
    # from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

    # l1b = datatree.open_datatree(l1b_path)
    # xs2tau = l1b['intraburst']['xspectra_2tau_Re']+1j*l1b['intraburst']['xspectra_2tau_Im']
    # xs2tau = xs2tau.assign_coords({'k_rg':xs2tau['k_rg'].mean(dim=set(xs2tau['k_rg'].dims)-set(['freq_sample']), keep_attrs=True)}).swap_dims({'freq_sample':'k_rg', 'freq_line':'k_az'})
    # xs2tau = symmetrize_xspectrum(xs2tau).squeeze(dim='2tau')

    # xsplot = list()
    # for mytile in xndindex({'burst':xs2tau.sizes['burst'], 'tile_line':xs2tau.sizes['tile_line'], 'tile_sample':xs2tau.sizes['tile_sample']}):

    #     heading = np.radians(l1b['intraburst']['ground_heading'][mytile].data)
    #     incidence = np.round(l1b['intraburst']['incidence'][mytile].item(),2)
    #     tau = np.round(l1b['intraburst']['tau'][mytile].item(),3)
    #     cutoff = int(l1b['intraburst'].ds[mytile]['azimuth_cutoff'].data)
    #     lon = np.round(l1b['intraburst'].ds[mytile]['longitude'].item(), 2)
    #     lat = np.round(l1b['intraburst'].ds[mytile]['latitude'].item(), 2)
    #     keast = xs2tau['k_rg']*np.cos(heading)+xs2tau['k_az']*np.sin(heading)
    #     knorth = xs2tau['k_az']*np.cos(heading)-xs2tau['k_rg']*np.sin(heading)
    #     keast.attrs.update({'long_name':'wavenumber in East direction', 'units':'rad/m'})
    #     knorth.attrs.update({'long_name':'wavenumber in North direction', 'units':'rad/m'})
    #     xs2tau = xs2tau.assign_coords({'k_east':keast,'k_north':knorth})
    #     heading=np.arctan2(np.sin(heading), np.cos(heading))
    #     range_rotation = -np.degrees(heading) if np.abs(np.degrees(heading))<=90 else -np.degrees(heading)+180
    #     azimuth_rotation = -np.degrees(heading)+90 if np.degrees(heading)>=0 else -np.degrees(heading)-90

    #     figreal = plt.figure(figsize=(10,10), tight_layout=True)
    #     xs2tau[mytile].real.plot(cmap=cmap, vmin=0, x='k_east', y='k_north')
    #     ax = plt.gca()
    #     for r in [400,200,100,50]:
    #         circle = plt.Circle((0, 0), 2*np.pi/r, color='k', fill=False, linestyle='--', linewidth=0.5)
    #         ax.add_patch(circle)
    #         plt.text(-np.sqrt(2)*np.pi/r-0.002,np.sqrt(2)*np.pi/r+0.002,'{} m'.format(r), rotation=45.,horizontalalignment='center',verticalalignment='center')
    #     for a in [-60,-30,30,60]:
    #         plt.plot([-0.2*np.cos(np.radians(a)), 0.2*np.cos(np.radians(a))],[-0.2*np.sin(np.radians(a)), 0.2*np.sin(np.radians(a))], color='k', linestyle='--', linewidth=0.5)
    #     plt.vlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)
    #     plt.hlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)
    #     xp = 0.14*np.cos(heading)
    #     yp = 0.14*np.sin(heading)
    #     plt.plot([-xp,xp], [yp,-yp], color='r') # range line
    #     plt.plot([yp,-yp], [xp,-xp], color='r') # azimuth line
    #     plt.plot(np.array([-xp,xp])+2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])+2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff upper line
    #     plt.plot(np.array([-xp,xp])-2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])-2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff lower line
    #     plt.axis('scaled')
    #     plt.xlim([-0.14,0.14])
    #     plt.ylim([-0.14,0.14])

    #     plt.text(0.9*xp,-0.9*yp+0.006,'range', color='r', rotation=range_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # range name
    #     plt.text(0.85*yp+0.006,0.85*xp,'azimuth', color='r', rotation=azimuth_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # azimuth name
    #     plt.text(-0.135,-0.135,'cut-off : {} m'.format(cutoff))
    #     plt.text(-0.135,0.13,'incidence : {} deg'.format(incidence))
    #     plt.text(0.065,0.13,'longitude : {} deg'.format(lon))
    #     plt.text(0.065,0.12,'latitude : {} deg'.format(lat))
    #     plt.text(-0.135,0.12,'tau : {} s'.format(tau))
    #     plt.title('X-spectrum (real)')

    #     canvas = FigureCanvas(figreal)
    #     width, height = figreal.get_size_inches() * figreal.get_dpi()
    #     canvas.draw()
    #     imagereal = np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(int(height), int(width), 3)[146:854,102:810,:]
    #     plt.close()
        
    #     # -----------------------------------------

    #     figimag = plt.figure(figsize=(10,10), tight_layout=True)
    #     xs2tau[mytile].imag.plot(cmap=PuOr, x='k_east', y='k_north')
    #     ax = plt.gca()
    #     for r in [400,200,100,50]:
    #         circle = plt.Circle((0, 0), 2*np.pi/r, color='k', fill=False, linestyle='--', linewidth=0.5)
    #         ax.add_patch(circle)
    #         plt.text(-np.sqrt(2)*np.pi/r-0.002,np.sqrt(2)*np.pi/r+0.002,'{} m'.format(r), rotation=45.,horizontalalignment='center',verticalalignment='center')
    #     for a in [-60,-30,30,60]:
    #         plt.plot([-0.2*np.cos(np.radians(a)), 0.2*np.cos(np.radians(a))],[-0.2*np.sin(np.radians(a)), 0.2*np.sin(np.radians(a))], color='k', linestyle='--', linewidth=0.5)
    #     plt.vlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)
    #     plt.hlines(0,-0.2,0.2, color='k', linestyle='--', linewidth=0.5)

    #     plt.plot([-xp,xp], [yp,-yp], color='r') # range line
    #     plt.plot([yp,-yp], [xp,-xp], color='r') # azimuth line
    #     plt.plot(np.array([-xp,xp])+2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])+2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff upper line
    #     plt.plot(np.array([-xp,xp])-2*np.pi/cutoff*np.sin(heading), np.array([yp,-yp])-2*np.pi/cutoff*np.cos(heading), color='k', linestyle='--') # cutoff lower line

    #     plt.axis('scaled')
    #     plt.xlim([-0.14,0.14])
    #     plt.ylim([-0.14,0.14])
    #     # plt.grid()
    #     plt.text(0.9*xp,-0.9*yp+0.006,'range', color='r', rotation=range_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # range name
    #     plt.text(0.85*yp+0.006,0.85*xp,'azimuth', color='r', rotation=azimuth_rotation, fontsize=15,horizontalalignment='center',verticalalignment='center') # azimuth name
    #     plt.text(-0.135,-0.135,'cuto-off : {} m'.format(cutoff))
    #     plt.text(-0.135,0.13,'incidence : {} deg'.format(incidence))
    #     plt.text(-0.135,0.12,'tau : {} s'.format(tau))
    #     plt.text(0.065,0.13,'longitude : {} deg'.format(lon))
    #     plt.text(0.065,0.12,'latitude : {} deg'.format(lat))
    #     plt.title('X-spectrum (imaginary)')

    #     canvas = FigureCanvas(figimag)
    #     width, height = figimag.get_size_inches() * figimag.get_dpi()
    #     canvas.draw()
    #     imageimag = np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(int(height), int(width), 3)[146:854,102:810,:]
    #     plt.close()

    #     xsreal = xr.DataArray(imagereal, dims=('azimuth','range','rgb'), name='xspectra_real').assign_coords(mytile)
    #     xsimag = xr.DataArray(imageimag, dims=('azimuth','range','rgb'), name='xspectra_imag').assign_coords(mytile)
    #     xs = xr.merge([xsreal, xsimag])
    #     xsplot.append(xs)
    # xsplot = xr.combine_by_coords([xs.expand_dims(['burst', 'tile_line','tile_sample']) for xs in xsplot])
    # return xsplot, l1b['intraburst'].attrs['orbit_pass']