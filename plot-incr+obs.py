import getopt
import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs
from cartopy import config
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

import netCDF4

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

  def plot(self, lons, lats, data=[]):
   #ax.coastlines(resolution='110m')
   #ax.gridlines()

    nrows = len(data)
    ncols = 1

   #set up the plot
    proj = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw=dict(projection=proj),
                            figsize=(11,8.5))
 
   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      axs[i].set_global()

      pvar = data[i]

     #print('Plot No. ', i)
     #print('\tpvar.shape = ', pvar.shape)

      cyclic_data, cyclic_lons = add_cyclic_point(pvar, coord=lons)

      cs=axs[i].contourf(cyclic_lons, lats, cyclic_data, transform=proj,
                         levels=self.clevs, extend=self.extend,
                         alpha=self.alpha, cmap=self.cmapname)
     #               cmap=self.cmapname, extend='both')

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
      axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

      if(self.add_obs_marker and len(self.obslonarray) > 0):
       #adding dotes:
       #dotes = self.basemap.plot(x, y, 'bo', markersize=12)
       #dotes = self.basemap.plot(x, y, 'bo', markersize=6)
        marker = 'x'
        size = 3
        color = 'green'
        dotes = axs[i].scatter(self.obslonarray[i], self.obslatarray[i],
                               marker=marker, s=size, color=color)

      axs[i].set_title(self.runname[i])

   #Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.8,
                        wspace=0.02, hspace=0.02)

   #Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.85])

   #Draw the colorbar
    cbar=fig.colorbar(cs, cax=cbar_ax, pad=self.pad, ticks=self.cblevs,
                      orientation='vertical')

    cbar.set_label(self.label, rotation=90)

   #Add a big title at the top
    plt.suptitle(self.title)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(self.imagename is None):
        imagename = 't_aspect.png'
      else:
        imagename = self.imagename
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

  def set_default(self):
    self.imagename = 'sample.png'

    self.runname = ['Halo', 'RR', 'RR - Halo']

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
    self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
   #self.cmapname = 'jet'

    self.clevs = np.arange(-0.2, 0.21, 0.01)
    self.cblevs = np.arange(-0.2, 0.3, 0.1)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Unit (C)'
    self.title = 'Temperature Increment'

    self.obslon = []
    self.obslat = []
    self.add_obs_marker = False

  def set_label(self, label='Unit (C)'):
    self.label = label

  def set_title(self, title='Temperature Increment'):
    self.title = title

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_imagename(self, imagename):
    self.imagename = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

  def set_obs_lonlat(self, lonarray, latarray):
    self.obslonarray = lonarray
    self.obslatarray = latarray

  def switch_marker_on(self):
    self.add_obs_marker = True

  def switch_marker_off(self):
    self.add_obs_marker = False

def get_obs_latlon_from_file(filename):
  ncfile = netCDF4.Dataset(filename, 'r')
  ncgroup = ncfile['MetaData']
  lat = ncgroup.variables['latitude'][:]
  lon = ncgroup.variables['longitude'][:]
  ncfile.close()

  lon = np.where(lon > 0, lon, lon+360.0)

  return lat, lon

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                                                'file1=', 'file2='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--file1'):
      file1 = a
    elif o in ('--file2'):
      file2 = a
    else:
      assert False, 'unhandled option'

#--------------------------------------------------------------------------------
  dir1 = '/work2/noaa/gsienkf/weihuang/jedi/per_core_timing/run/solver_halo_aircraft/run_80.40t1n_36p'
  dir2 = '/work2/noaa/gsienkf/weihuang/jedi/per_core_timing/run/solver_RoundRobin_aircraft/run_80.40t1n_36p'

  file1 = '%s/analysis/increment/xainc.20200110_030000z.nc4' %(dir1)
  file2 = '%s/analysis/increment/xainc.20200110_030000z.nc4' %(dir2)

#--------------------------------------------------------------------------------
  obslonarry = []
  obslatarry = []
  obs1 = '%s/ioda_v2_data/aircraft_tsen_obs_2020011006_0000.nc4' %(dir1)
  obslat, obslon = get_obs_latlon_from_file(obs1)
  obslonarry.append(obslon)
  obslatarry.append(obslat)
  obs2 = '%s/ioda_v2_data/aircraft_tsen_obs_2020011006_0000.nc4' %(dir2)
  obslat, obslon = get_obs_latlon_from_file(obs2)
  obslonarry.append(obslon)
  obslatarry.append(obslat)

#--------------------------------------------------------------------------------

  gp = GeneratePlot(debug=debug, output=output)

  nc1 = netCDF4.Dataset(file1)
  nc2 = netCDF4.Dataset(file2)
  lats = nc1.variables['lat'][:]
  lons = nc1.variables['lon'][:]

#-----------------------------------------------------------------------------------------
  gp.set_obs_lonlat(obslonarry, obslatarry)
  gp.switch_marker_on()

#-----------------------------------------------------------------------------------------
  clevs = np.arange(-1.0, 1.01, 0.01)
  cblevs = np.arange(-1.0, 1.1, 0.1)

  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#-----------------------------------------------------------------------------------------
 #varlist1 = ['T', 'ua', 'va', 'sphum', 'delp', 'DZ', 'o3mr']
 #varlist2 = ['T', 'ua', 'va', 'sphum', 'delp', 'DZ', 'o3mr']
  varlist1 = ['T', 'delp', 'sphum']
  varlist2 = ['T', 'delp', 'sphum']

  unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)',
              'Unit (kg/kg)', 'Unit (Pa', 'Unit (m', 'Unit (ppm)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist1)):
    var1 = nc1.variables[varlist1[n]][0,:,:,:]
    var2 = nc2.variables[varlist2[n]][0,:,:,:]

    nlev, nlat, nlon = var1.shape
    print('var1.shape = ', var1.shape)
    print('var2.shape = ', var2.shape)

    gp.set_label(unitlist[n])

    for lev in range(5, nlev, 10):
      v1 = var1[lev,:,:]
      v2 = var2[lev,:,:]

      data = [v1, v2]

      title = '%s at Level %d' %(varlist1[n], lev)
      gp.set_title(title)

      print('Plotting ', title)
      print('\tv1.shape = ', v1.shape)
      print('\tv2.shape = ', v2.shape)
 
      print('\tv1.max: %f, v1.min: %f' %(np.max(v1), np.min(v1)))
      print('\tv2.max: %f, v2.min: %f' %(np.max(v2), np.min(v2)))

      imagename = '%s_lev_%3.3d.png' %(varlist1[n], lev)
      gp.set_imagename(imagename)

      gp.plot(lons, lats, data=data)

#-----------------------------------------------------------------------------------------
  nc1.close()
  nc2.close()

