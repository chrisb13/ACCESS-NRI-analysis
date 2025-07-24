#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  ACCESS-NRI, Australia’s climate simulator 
#                 Building 69, 5 Liversidge St., Australian National University ACT 
#   Contact: chris.bull@anu.edu.au
#   Date created: Fri, 23 Jul 2025 10:00:12
#   Machine created on: SB2Vbox
#

"""
Quick script to plot out salinity differences between the different runs with modified bathymetry
"""
import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib import colors as c
import datetime
import modify_med_restart as medtool
import glob
from matplotlib import gridspec
import matplotlib.pyplot as plt
import collections
import itertools

def mk_filled_array_from_mask(array,oceanpts):
    array_fill=np.ma.masked_where(array==True,np.zeros(np.shape(oceanpts))) 
    array_fill.fill_value=1
    array_fill=array_fill.T
    return array_fill.filled()

def mk_salty_tseries(ifiles):
    #hack to make this run a bit faster for testing!
    fls=xr.open_mfdataset(ifiles[0:2])
    fls=xr.open_mfdataset(ifiles)

    data=fls['surface_salt'].values*wholemedf
    data[data == 0] = np.nan
    medmean=np.nanmean(np.nanmean(data,axis=-1),axis=-1)

    data=fls['surface_salt'].values*westmedf
    data[data == 0] = np.nan
    westmean=np.nanmean(np.nanmean(data,axis=-1),axis=-1)

    data=fls['surface_salt'].values*eastmedf
    data[data == 0] = np.nan
    eastmean=np.nanmean(np.nanmean(data,axis=-1),axis=-1)

    return {'fs':fls,'medmean':medmean,'westmean':westmean,'eastmean':eastmean}

def pmedtseries(pltax,pdict,color,label):
    pltax.plot(pdict['medmean'],label= label+' whole med',color=color,alpha=0.7,lw=2)
    pltax.plot(pdict['westmean'],label=label+' west',linestyle='dashed',color=color,alpha=0.6,lw=2)
    pltax.plot(pdict['eastmean'],label=label+' east',linestyle='dotted',color=color,alpha=0.6,lw=2)
    return

if __name__ == "__main__": 
    ifol='/home/561/cyb561/esm16/'
    infile=ifol+'grid_spec.nc'
    assert(os.path.exists(infile)),"netCDF file does not exist!"
    ifile=xr.open_dataset(infile)
    x,y,pme,wholemed,westmed,eastmed,oceanpts=medtool.mk_med_mask(ifile)

    wholemedf=mk_filled_array_from_mask(wholemed,oceanpts)
    westmedf =mk_filled_array_from_mask(westmed,oceanpts)
    eastmedf =mk_filled_array_from_mask(eastmed,oceanpts)

    
    ##########################################
    #  define input folders for experiments  #
    ##########################################
    
    ifiles=sorted(glob.glob('/scratch/p66/jxs599/access-esm/archive/JuneSpinUp-JuneSpinUp-bfaa9c5b/output0*/ocean/ocean-2d-surface_salt-1monthly-mean-ym_*.nc'))
    assert(ifiles!=[]),"glob didn't find anything!"
    JuneSpinUp=mk_salty_tseries(ifiles[0:30])

    ifiles=sorted(glob.glob('/scratch/tm70/cyb561/access-esm/archive/cntrlR_01-cntrlR_01-09d5ffea/output*/ocean/ocean-2d-surface_salt-1monthly-mean-ym_*.nc'))

    assert(ifiles!=[]),"glob didn't find anything!"
    cntrlr_01=mk_salty_tseries(ifiles)

    ifiles=sorted(glob.glob('/scratch/tm70/cyb561/access-esm/archive/sicilyR_01-sicilyR_01-189d463d/output*/ocean/ocean-2d-surface_salt-1monthly-mean-ym_*.nc'))
    sicilyr_01=mk_salty_tseries(ifiles)

    ifiles=sorted(glob.glob('/scratch/tm70/cyb561/access-esm/archive/sardiniaR_01-sardiniaR_01-02261d49/output*/ocean/ocean-2d-surface_salt-1monthly-mean-ym_*.nc'))
    sardiniar_01=mk_salty_tseries(ifiles)

    #for spatial plot
    pme_cntrlr01=cntrlr_01['fs']['surface_salt']
    pme_sicilyr_01=sicilyr_01['fs']['surface_salt']-cntrlr_01['fs']['surface_salt']
    pme_sardiniar_01=sardiniar_01['fs']['surface_salt']-cntrlr_01['fs']['surface_salt']

    #################################
    #  spatial plots on their own   #
    #################################
    
    #plt.close('all')
    #fig=plt.figure()
    #ax=fig.add_subplot(1, 3,1)

    #cs1=ax.pcolormesh(y,x,pme_cntrlr01[tstep,:].T)
    #ax.set_title('cntrlr_01')
    #plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    #ax.set_xlim([-10,40])
    #ax.set_ylim([25,49])

    #ax=fig.add_subplot(1, 3,2)
    #cs1=ax.pcolormesh(y,x,pme_sicilyr_01[tstep,:].T,cmap='seismic',vmin=-0.4,vmax=0.4)
    #ax.set_title('pme_sicilyr_01 - cntrlr_01')
    #plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    #ax.set_xlim([-10,40])
    #ax.set_ylim([25,49])

    #ax=fig.add_subplot(1, 3,3)
    #cs1=ax.pcolormesh(y,x,pme_sardiniar_01[tstep,:].T,cmap='seismic',vmin=-0.4,vmax=0.4)
    #ax.set_title('pme_sardiniar_01 - cntrlr_01')
    #plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    #ax.set_xlim([-10,40])
    #ax.set_ylim([25,49])

    #plt.show()

    #############################
    #  time series on it's own  #
    #############################

    #plt.close('all')
    #fig=plt.figure()
    #ax=fig.add_subplot(1, 1,1)
    #pmedtseries(ax,cntrlr_01,'green','cntrlR_0')
    #pmedtseries(ax,JuneSpinUp,'purple','JuneSpinUp')
    ##pmedtseries(ax,sicilyr_01,'blue','sicilyR_01')
    ##pmedtseries(ax,sardiniar_01,'red','sardiniarR_01')
    #ax.set_title("Comparing: cntrlR_0, JuneSpinUp med' salinity")
    #ax.legend()
    #ax.grid()
    #ax.set_xlabel('Months')
    #ax.set_ylabel('Surface salinity')
    #plt.show()

    #tstep=30

    print("Let's make a movie!!")
    print("Let's make a movie!!")
    print("Let's make a movie!!")

    for tstep in np.arange(360):
        print(tstep)
        plt.close('all')
        row=4
        col=3

        fig=plt.figure(figsize=(5.0*col,3*row))
        gs = gridspec.GridSpec(row, col,height_ratios=[1,1,.00005,.1],width_ratios=[1]*col,hspace=.3,wspace=0.25)

        #timeseries
        ax = plt.subplot(gs[0,0:3])
        # vertical line
        ax.vlines(x=tstep, ymin=37.7, ymax=38.2,linewidth=2, color='black', zorder=1,alpha=.8)
        pmedtseries(ax,cntrlr_01,'green','cntrlR_0')
        pmedtseries(ax,sicilyr_01,'blue','sicilyR_01')
        pmedtseries(ax,sardiniar_01,'red','sardiniarR_01')
        ax.legend(loc='upper left',ncol=3)
        ax.grid()
        ax.set_xlabel('Months')
        ax.set_ylabel('Average Surface Salinity')

        #spatial maps
        ax = plt.subplot(gs[1,0])
        cs1_cntrlr01=ax.pcolormesh(y,x,pme_cntrlr01[tstep,:].T,vmin=36.5, vmax=40)
        ax.set_title('cntrlr_01')
        #plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
        ax.set_xlim([-10,40])
        ax.set_ylim([25,49])

        ax = plt.subplot(gs[1,1])
        cs1_sicilyr=ax.pcolormesh(y,x,pme_sicilyr_01[tstep,:].T,cmap='seismic',vmin=-0.6,vmax=0.6)
        ax.set_title('pme_sicilyr_01 - cntrlr_01')
        #plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
        ax.set_xlim([-10,40])
        ax.set_ylim([25,49])

        ax = plt.subplot(gs[1,2])
        cs1_sardiniar=ax.pcolormesh(y,x,pme_sardiniar_01[tstep,:].T,cmap='seismic',vmin=-0.6,vmax=0.6)
        ax.set_title('pme_sardiniar_01 - cntrlr_01')
        #plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
        ax.set_xlim([-10,40])
        ax.set_ylim([25,49])

        #colourbars
        #ax = plt.subplot(gs[3,0])

        cbar=plt.colorbar(cs1_cntrlr01,cax=plt.subplot(gs[3,0]),orientation='horizontal')
        cbar=plt.colorbar(cs1_sicilyr,cax=plt.subplot(gs[3,1]),orientation='horizontal')
        cbar=plt.colorbar(cs1_sardiniar,cax=plt.subplot(gs[3,2]),orientation='horizontal')
        #plt.show()

        fig.savefig('./medsurface_salinity_movie_'+str(tstep).zfill(5)+'.png',dpi=300,bbox_inches='tight')


