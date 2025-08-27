from matplotlib import pyplot as plt
import numpy as np

def fig_MACS_row_incidence_multi_azi_wrt_windspeed(df,satellite='S1A+B',part='Re',burstkind='intraburst',polarization='vv',lambda_val='50',ymax=0.4):
    incidences = [34,39,44]
    #azicases = {'upwind':0,'downwind':180,'crosswind':90}
    az_wdir = df["wdir_az_scat"]
    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))
    #az_up = az_wdir[(az_wdir < 0.5) | (az_wdir > 359.5)];
    az_up = (az_wdir < 0.5) | (az_wdir > 359.5)
    #az_down = az_wdir[(az_wdir <= 180.5) & (az_wdir > 179.5)];
    az_down = (az_wdir <= 180.5) & (az_wdir > 179.5)
    #az_cross = az_wdir[(abs(az_wdir - 90) < 0.5) | (abs(az_wdir - 270) < 0.5)];
    az_cross = (abs(az_wdir - 90) < 0.5) | (abs(az_wdir - 270) < 0.5)
    azicases = {'upwind':az_up,'downwind':az_down,'crosswind':az_cross}
    deltaws = 0.7
    windspeed_vect = np.arange(0,28,deltaws)
    fig, ax = plt.subplots(len(incidences),1, figsize=(15, 15))
    fig.suptitle(satellite+' | IMACS versus wind speed | %s %s' % (burstkind,polarization), fontsize=15)
    for ii,inc in enumerate(incidences):
        print('incidence',inc)
        cond_inc = (abs(df['incidence']-inc)<=0.5)
        ax[ii].axhline(y=0,c='k',alpha=.7)
        #print('cond_inc',cond_inc.sum())
        for aa,azik in enumerate(azicases):
            az_sub = df.loc[df.index.intersection(azicases[azik].index)]
            #print(len(az_sub.index))

            inc_sub = df.loc[df.index.intersection(cond_inc.index)]
            print(len(inc_sub.index))
            #sub = az_sub.loc[az_sub.index.intersection(inc_sub.index)]
            sub = df[cond_inc & azicases[azik]]
            #print('sub',type(sub),len(sub))
            mean_imacs = []
            std_imacs = []
            for ww in windspeed_vect:
                condwindspeed = (abs(df['Wspeed'][sub.index]-ww)<deltaws)
                #print('condwindspeed',condwindspeed.sum())
                subw = sub[condwindspeed]
                #subw = sub.loc[sub.index.intersection(condwindspeed.index)]
                #print('subw',subw)
                mean_imacs.append(df[varname][subw.index].mean())
                std_imacs.append(df[varname][subw.index].std())
            mean_imacs = np.array(mean_imacs)
            std_imacs = np.array(std_imacs)
            #ax[ii].plot(df['Wspeed'][sub.index],df[varname][sub.index],'.',label=azik)
            line, = ax[ii].plot(windspeed_vect,mean_imacs,'.-',lw=3,markersize=5,label='obs '+azik+':%s'%len(sub.index),alpha=0.9)
            color = line.get_color()
            ax[ii].fill_between(windspeed_vect,mean_imacs-std_imacs,mean_imacs+std_imacs,alpha=0.4,color=color)
        ax[ii].legend(fontsize=12)
        ax[ii].grid(True)
        ax[ii].set_ylabel('IMACS $\lambda$=%sm []'%lambda_val)
        ax[ii].set_xlabel('Wind speed [m/s]')
        ax[ii].set_ylim(-ymax,ymax)
        ax[ii].set_title('incidence: %1.1f° +/- 0.5°'%inc)
    return fig,ax,incidences
  
