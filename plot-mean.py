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

from netCDF4 import Dataset as netcdf_dataset

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

  def plot(self, lons, lats, data):
   #ax.coastlines(resolution='110m')
   #ax.gridlines()

    nrows = 1
    ncols = 1

   #set up the plot
    proj = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw=dict(projection=proj),
                            figsize=(11,8.5))
    axs.set_global()
    pvar = data

    cyclic_data, cyclic_lons = add_cyclic_point(pvar, coord=lons)

    cs=axs.contourf(cyclic_lons, lats, cyclic_data, transform=proj,
                    levels=self.clevs, extend=self.extend,
                    alpha=self.alpha, cmap=self.cmapname)

    axs.set_extent([-180, 180, -90, 90], crs=proj)
    axs.coastlines(resolution='auto', color='k')
    axs.gridlines(color='lightgrey', linestyle='-', draw_labels=True)
    axs.set_title(self.runname[0])

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
   #self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
    self.cmapname = 'rainbow'
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

  def set_runname(self, runname):
    self.runname = runname

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
 #datadir = '/work2/noaa/gsienkf/weihuang/jedi/develop/build/fv3-jedi/test/Data/analysis/letkf/gfs/mts'
  datadir = '/work2/noaa/gsienkf/weihuang/jedi/develop/build/fv3-jedi/test'

  filename = 'letkf.mid.00020201215_000000z.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                                                'filename=', 'datadir='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--datadir'):
      datadir = a
    elif o in ('--filename'):
      filename = a
    else:
      assert False, 'unhandled option'

  gp = GeneratePlot(debug=debug, output=output)

  fullpath = '%s/%s' %(datadir, filename)
  ncfile = netcdf_dataset(fullpath)

  lats = ncfile.variables['lat'][:]
  lons = ncfile.variables['lon'][:]

#-----------------------------------------------------------------------------------------
  clevs = np.arange(200.0, 306.0, 1.0)
  cblevs = np.arange(200.0, 310.0, 5.0)

  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#-----------------------------------------------------------------------------------------
  varlist = ['T', 'ua', 'va', 'sphum', 'delp', 'DZ', 'o3mr']
  unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)',
              'Unit (kg/kg)', 'Unit (Pa', 'Unit (m', 'Unit (ppm)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist)):
    var = ncfile.variables[varlist[n]][0,:, :, :]
    gp.set_runname([varlist[n]])

    nlev, nlat, nlon = var.shape
    print('var.shape = ', var.shape)

    gp.set_label(unitlist[n])

    for lev in range(5, nlev, 10):
      data = var[lev,:,:]

      title = '%s at Level %d' %(varlist[n], lev)
      gp.set_title(title)

      print('Plotting ', title)
      print('\tdata.shape = ', data.shape)
 
      print('\tdata.max: %f, data.min: %f' %(np.max(data), np.min(data)))

      imagename = '%s_lev_%3.3d.png' %(varlist[n], lev)
      gp.set_imagename(imagename)

      gp.plot(lons, lats, data)

#-----------------------------------------------------------------------------------------
  ncfile.close()

