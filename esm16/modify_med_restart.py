#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  ACCESS-NRI, Australia’s climate simulator 
#                 Building 69, 5 Liversidge St., Australian National University ACT 
#   Contact: chris.bull@anu.edu.au
#   Date created: Fri, 18 Jul 2025 06:50:12
#   Machine created on: SB2Vbox
#

"""
Quick script to modify the restart files for esm1.6 mediternean sea problem
"""
import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib import colors as c
import datetime

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx

def pbathy(paxis):
    cs1=paxis.pcolormesh(y,x,pme,vmin=0,vmax=5000,alpha=0.5)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    paxis.set_xlim([-10,40])
    paxis.set_ylim([25,49])
    paxis.grid(True)
    return

def mkmask(bathyarr):
    return np.ma.masked_where(bathyarr==0,bathyarr) 

def mk_med_mask(ifile):
    yT,xT=ifile['grid_y_T'].values,ifile['grid_x_T'].values
    x,y=np.meshgrid(yT,xT)
    pme=mkmask(ifile['depth_t'].T)

    #pme[100:150,200:220]=-10000000000

    yTs={'30':find_nearest(yT,30),'46':find_nearest(yT,46)}
    xTs={'-8':find_nearest(xT,-8),'38':find_nearest(xT,38)}
    yTidx=[yTs['30'][1],yTs['46'][1]]
    xTidx=[xTs['-8'][1],xTs['38'][1]]
    xT_0=find_nearest(xT,0)
    #where does y go to 25 and x to zero?
    oceanpts=~pme.mask

    #grab the whole mediterranean
    wholemed=np.zeros(np.shape(oceanpts),dtype=bool)
    wholemed[xTidx[0]:xTidx[1],yTidx[0]:yTidx[1]]=True
    wholemed[xTidx[0]:xTidx[0]+9,yTidx[1]-3:yTidx[1]]=False #remove top-left corner that grabs the Atlantic
    wholemed=wholemed&oceanpts

    westmed=np.zeros(np.shape(oceanpts),dtype=bool)
    westmed[xTidx[0]:xTidx[0]+19,yTidx[0]:yTidx[1]]=True
    westmed[xTidx[0]:xTidx[0]+9,yTidx[1]-3:yTidx[1]]=False #remove top-left corner that grabs the Atlantic

    #western med funny shaped inlet
    westmed[xTidx[0]+19:xTidx[0]+24,yTidx[1]-9:yTidx[1]-3]=True
    westmed[xTidx[0]+24-1,yTidx[1]-3-1]=False   #we go too far so we have to undo this point in the eastern med
    westmed=westmed&oceanpts

    eastmed=np.zeros(np.shape(oceanpts),dtype=bool)
    eastmed=wholemed&~westmed
    return x,y,pme,wholemed, westmed,eastmed,oceanpts

if __name__ == "__main__": 

    ifol='/home/561/cyb561/esm16/'
    infile=ifol+'grid_spec.nc'
    assert(os.path.exists(infile)),"netCDF file does not exist!"
    ifile=xr.open_dataset(infile)

    plt.close('all')
    fig=plt.figure()
    ax=fig.add_subplot(1, 3,1)

    x,y,pme,wholemed,westmed,eastmed,oceanpts=mk_med_mask(ifile)

    pbathy(ax)
    ifile['depth_t'].T
    pme=mkmask(ifile['depth_t'].T)
    pme[wholemed]=-10000000000
    pme=np.ma.masked_where(pme!=-10000000000,pme) 
    cMap = c.ListedColormap(['r','b','m'])
    cs1=ax.pcolormesh(y,x,pme.mask,cmap=cMap,alpha=0.5)
    ax.set_title('Whole med')

    ax=fig.add_subplot(1, 3,2)
    pbathy(ax)
    pme=mkmask(ifile['depth_t'].T)
    pme[eastmed]=-10000000000
    pme=np.ma.masked_where(pme!=-10000000000,pme) 
    cMap = c.ListedColormap(['r','b','m'])
    cs1=ax.pcolormesh(y,x,pme.mask,cmap=cMap,alpha=0.5)
    ax.set_title('East med')

    ax=fig.add_subplot(1, 3,3)
    pbathy(ax)
    pme=mkmask(ifile['depth_t'].T)
    pme[westmed]=-10000000000
    pme=np.ma.masked_where(pme!=-10000000000,pme) 
    cMap = c.ListedColormap(['r','b','m'])
    cs1=ax.pcolormesh(y,x,pme.mask,cmap=cMap,alpha=0.5)
    ax.set_title('West med')

    #plt.show()

    before_long=ifol+'ocean_temp_salt.res.nc-19501231'
    after_long=ifol+'ocean_temp_salt.res.nc-10001231'
    ifile_before=xr.open_dataset(before_long)
    ifile_after=xr.open_dataset(after_long)
    #source
    #/g/data/tm70/cyb561/esm16_bathy_tests/2025.06.12/ocean/ocean_temp_salt.res.nc

    #https://access-nri.zulipchat.com/#narrow/channel/470332-ocean-seaice/topic/Mediterranean.20sea.20in.20esm1.2E6/near/528205171
    #based on Spencer Wong
    #You could clone from this tag https://github.com/ACCESS-NRI/access-esm1.6-configs/tree/20250612-spinup-dev-preindustrial%2Bconcentrations to get the config as at the start of JuneSpinup
    rfile=xr.open_dataset(ifol+'ocean_temp_salt.res.nc')


    plt.close('all')
    fig=plt.figure()
    ax=fig.add_subplot(1, 3,1)

    #these are a bit weird so just using the old values....
    #yT,xT=ifile_before['yaxis_1'].values,ifile_before['xaxis_1'].values
    #x,y=np.meshgrid(yT,xT)

    pme=mkmask(ifile_before['salt'][0,0,:,:].T)
    cs1=ax.pcolormesh(y,x,pme,vmin=30,vmax=45,alpha=1)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    ax.set_xlim([-10,40])
    ax.set_ylim([25,49])
    ax.set_title('ocean_temp_salt.res.nc-19501231')

    ax=fig.add_subplot(1, 3,2)
    pme=mkmask(ifile_after['salt'][0,0,:,:].T)
    cs1=ax.pcolormesh(y,x,pme,vmin=30,vmax=45,alpha=1)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    ax.set_xlim([-10,40])
    ax.set_ylim([25,49])
    ax.set_title('ocean_temp_salt.res.nc-19501231')

    ax=fig.add_subplot(1, 3,3)
    pme=mkmask(rfile['salt'][0,0,:,:].T)
    cs1=ax.pcolormesh(y,x,pme,vmin=30,vmax=45,alpha=1)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    ax.set_xlim([-10,40])
    ax.set_ylim([25,49])
    ax.set_title('rfile')

    #plt.show()

    wholemed_fill=np.ma.masked_where(wholemed==True,np.zeros(np.shape(oceanpts))) 
    wholemed_fill.fill_value=1
    wholemed_fill=wholemed_fill.T

    wholenmed_fill=np.ma.masked_where(wholemed==False,np.zeros(np.shape(oceanpts))) 
    wholenmed_fill.fill_value=1
    wholenmed_fill=wholenmed_fill.T

    #remove the old salt
    salt_none=rfile['salt'][0,:,:,:]*wholenmed_fill.filled()

    #put in the new salt
    salt=(ifile_before['salt'][0,:,:,:]*wholemed_fill.filled())+salt_none
    #print(np.sum(salt))
    rfile['salt'][0,:,:,:]=salt
    #print(np.sum(rfile['salt'][0,:,:,:]))


    #check that it worked..
    plt.close('all')
    fig=plt.figure()
    ax=fig.add_subplot(1, 3,1)

    pme=mkmask(salt_none[0,:].T)
    cs1=ax.pcolormesh(y,x,pme,vmin=30,vmax=45,alpha=1)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    ax.set_xlim([-10,40])
    ax.set_ylim([25,49])
    ax.set_title('did the salt get removed?')

    ax=fig.add_subplot(1, 3,2)
    pme=mkmask(salt[0,:].T)
    cs1=ax.pcolormesh(y,x,pme,vmin=30,vmax=45,alpha=1)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    ax.set_xlim([-10,40])
    ax.set_ylim([25,49])
    ax.set_title('has the salt been inserted in the salt field?')

    ax=fig.add_subplot(1, 3,3)
    pme=mkmask(rfile['salt'][0,0,:,:].T)
    cs1=ax.pcolormesh(y,x,pme,vmin=30,vmax=45,alpha=1)
    plt.colorbar(cs1,cax=make_axes_locatable(ax).append_axes("bottom", size="5%", pad=0.25),orientation='horizontal')
    ax.set_xlim([-10,40])
    ax.set_ylim([25,49])
    ax.set_title('has the salt been assigned to the old xarray?')

    plt.show()

    ##we need to write out salt
    rfile.attrs['history']=   r'Created by /home/561/cyb561/esm16/modify_med_restart.py'
    rfile.attrs['Comment']=   r'This file was created to reset the salinity back to what was in /g/data/access/projects/access/data/ACCESS_CMIP5/restart/hPI-C01a/ocn/ocean_temp_salt.res.nc-19501231 for the Mediterranean sea in esm1.6, see discussion https://access-nri.zulipchat.com/#narrow/channel/470332-ocean-seaice/topic/Mediterranean.20sea.20in.20esm1.2E6/near/528205419.'
    rfile.attrs['contact']=   r'chris.bull@anu.edu.au'
    rfile.attrs['creation_date']=   str(datetime.datetime.now())

    rfile.to_netcdf('./ocean_temp_salt.res_medreset.nc')
    rfile.close()
    print('./ocean_temp_salt.res_medreset.nc written')

