#=========================================================================
import os
import sys
import types
import getopt

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from netCDF4 import Dataset
#from matplotlib import cm

class FileReader():
  def __init__(self, debug=0, datafiles=[], gridspecfiles=[]):
    self.debug = debug

    if(self.debug):
      print('debug = ', debug)

    self.setDataFiles(datafiles=datafiles)
    self.setGridSpecFiles(gridspecfiles=gridspecfiles)
    self.has_snd_file = 0
    self.snd_files = []

  def setDataFiles(self, datafiles=[]):
    self.datafiles = datafiles

  def setSecondFiles(self, files):
    self.has_snd_file = 1
    self.snd_files = files

  def setGridSpecFiles(self, gridspecfiles=[]):
    self.gridspecfiles = gridspecfiles

  def read3Dvar(self,datafiles,varname,ntime=0):
    """
    read FV3 cubed sphere 3D data.

    datafiles : list of data filenames for each tile
    varname : var name to read from data files

    returns data array"""
    data = None
    nt = len(self.datafiles)
    for it in range(nt):
      datafile = datafiles[it]
      print('reading ',datafile)
      nc = Dataset(datafile)
      arr = nc.variables[varname][ntime,:,:,:]
      nz, ny, nx = arr.shape

      print('arr.ndim=', arr.ndim)
      print('arr.shape=', arr.shape)
      print('arr.size=', arr.size)

      if(data is None):
        data = np.zeros((nt, nz, ny, nx))

      data[it,:,:,:] = arr[:,:,:]
      nc.close()

    return data

  def readTileInfo(self,datafiles,varname):
    var1d = []
    nt = len(self.datafiles)
    for it in range(nt):
      datafile = datafiles[it]
      nc = Dataset(datafile)
      arr = nc.variables[varname][:,:]
      ny, nx = arr.shape

     #print('arr.ndim=', arr.ndim)
     #print('arr.shape=', arr.shape)
     #print('arr.size=', arr.size)

      varc = np.zeros(((ny-1),(nx-1)))
      varc[0:ny-1,0:nx-1] = 0.25*(arr[0:ny-1,0:nx-1] + arr[0:ny-1,1:nx] + arr[1:ny,0:nx-1] + arr[1:ny,1:nx])

      varc1d = np.reshape(varc, ((nx-1)*(ny-1),))
      var1d.extend(varc1d)

      nc.close()

    return var1d

  def get_original_data(self, varname):
    lat1d, lon1d = self.get_GridSpec_latlon()
    
    print('len(lon1d) = ', len(lon1d))
    print('len(lat1d) = ', len(lat1d))

   #varname = 'T'
    if(self.has_snd_file):
      var1 = self.read3Dvar(self.datafiles, varname)
      var2 = self.read3Dvar(self.snd_files, varname)
      var = var2 - var1
    else:
      var = self.read3Dvar(self.datafiles, varname)

    print('var.ndim=', var.ndim)
    print('var.shape=', var.shape)
    print('var.size=', var.size)

    return lat1d, lon1d, var

  def get_latlon_tile(self, varname, nlon=360, nlat=181, method='linear'):
    lat1d, lon1d = self.readGridSpecFiles()

    var1d = self.readTileInfo(self.datafiles, varname)

    olons,olats,latlon_var = self.interp_to_latlon(lon1d, lat1d, var1d, nlon=nlon, nlat=nlat, method=method)

    return latlon_var

  def get_GridSpec_latlon(self):
    '''
    gridspecfiles : list of grid_spec filenames for each tile.
    '''

    nc = 0
    nt = len(self.gridspecfiles)
    lon1d = []
    lat1d = []
    for gridspecfile in self.gridspecfiles:
      print('reading ',gridspecfile)
      nc = Dataset(gridspecfile)
      lons = nc.variables['x'][:]
      lats = nc.variables['y'][:]

     #print('lons.ndim=', lons.ndim)
     #print('lons.shape=', lons.shape)
     #print('lons.size=', lons.size)

      ny, nx = lons.shape
      latc = np.zeros(((ny-1),(nx-1)))
      lonc = np.zeros(((ny-1),(nx-1)))
      latc[0:ny-1,0:nx-1] = 0.25*(lats[0:ny-1,0:nx-1] + lats[0:ny-1,1:nx] + lats[1:ny,0:nx-1] + lats[1:ny,1:nx])
      lonc[0:ny-1,0:nx-1] = 0.25*(lons[0:ny-1,0:nx-1] + lons[0:ny-1,1:nx] + lons[1:ny,0:nx-1] + lons[1:ny,1:nx])

      print('lonc.ndim=', lonc.ndim)
      print('lonc.shape=', lonc.shape)
      print('lonc.size=', lonc.size)

     #lonc1d = np.reshape(lonc, ((nx-1)*(ny-1),))
     #latc1d = np.reshape(latc, ((nx-1)*(ny-1),))

      lon1d.append(lonc)
      lat1d.append(latc)

      nc.close()

    return lat1d, lon1d

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C192'

  datadir = '/work/noaa/da/menetrie/StaticBTraining/c384/bump_gsi/dirac_cor_local'

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%s/C192_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s/20200131.000000.dirac_SABER.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)

  fr = FileReader(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)

 #varname = 't'
 #varname = 'chi'
  varname = 'psi'
  lat, lon, var = fr.get_original_data(varname)

  print('var.ndim = ', var.ndim)
  print('var.shape = ', var.shape)

  nt, nz, ny, nx = var.shape
  nlen = int(nx*ny)
  precision = 1

 #nq = int(nx/4)
 #nq = 100
  nq = 160
 #nq = 180

 #cmapname = coolwarm, bwr, rainbow, jet, seismi
  cmapname = 'rainbow'
  clevs = np.arange(0.05, 1.05, 0.05)
  cblevs = np.arange(0.1, 1.3, 0.2)

  for it in range(nt):
    imgname = 'tile_%d.png' %(it)
    title = 'Temperature on tile %d' %(it)

   #plon = lon[it]
   #plat = lat[it]

    for ilev in range(0, nz, 5):
     #pvar = var[it,ilev,:,:]
      pvar = var[it,ilev,nq:-nq,nq:-nq]

     #if(np.max(pvar) > 1.0e-10):
      if(np.max(pvar) > 0.2):
        title = '%s min: %8.5f, max: %8.3f on tile: %d, at level: %d' %(varname, np.min(pvar), np.max(pvar), it, ilev)
        print(title)

        fig, ax = plt.subplots()
        contplot = ax.contourf(pvar, levels=clevs, extend='both',
                                     alpha=0.5, cmap=cmapname)
        ax.clabel(contplot, inline=True, fontsize=10)
        ax.set_title(title)

        cb = fig.colorbar(contplot, orientation='horizontal',
                          pad=0.1, ticks=cblevs)

        cb.set_label(label=varname, size='large', weight='bold')

        cb.ax.tick_params(labelsize='medium')
        if(precision == 0):
          cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in cblevs], minor=False)
        elif(precision == 1):
          cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in cblevs], minor=False)
        elif(precision == 2):
          cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in cblevs], minor=False)
        else:
          cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in cblevs], minor=False)

        if(output):
          image_name='%s_tile%d_level%d.png' %(varname, it, ilev)
          plt.tight_layout()
          print('Saving image as ', image_name)
          plt.savefig(image_name)
        else:
          plt.show()

