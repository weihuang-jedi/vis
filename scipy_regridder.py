import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata as regridder

class RegridFV3():
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
   #if(len(self.gridspecfiles)):
   #  self.readGridSpecFiles()
   #else:
   #  print('No gridspecfiles specified. exit.')
   #  sys.exit(-1)

  def interp_to_latlon(self, lon_1d, lat_1d, var_1d, nlon=360, nlat=181, method='linear'):
    '''
    Interpolate a variable on cube-sphere grid (such as FV3) to LatLon grid
    '''
    dlon = 360.0/nlon
    dlat = 180.0/(nlat - 1)
   #Create a lat-lon uniform grid
    out_lon = np.arange(0.0, 360.0, dlon)
    out_lat = np.arange(-90.0, 91.0, dlat)
    lon, lat = np.meshgrid(out_lon, out_lat)

   #print('out_lon.size = ', out_lon.size)
   #print('lon.size = ', lon.size)

   #print('lon_1d.ndim = ', lon_1d.ndim)
   #print('lon_1d.size = ', lon_1d.size)
   #print('lon_1d.shape = ', lon_1d.shape)

   #print('var_1d.ndim = ', var_1d.ndim)
   #print('var_1d.size = ', var_1d.size)
   #print('var_1d.shape = ', var_1d.shape)

    # Interpolate from cube to lat-lon grid
    out_var = regridder((lon_1d,lat_1d), var_1d,
                        (lon,lat), method=method)

    nlen = int(nlon*nlat)
   #print('nlen = ', nlen)
   #print('lon.size = ', lon.size)
    olon = np.reshape(lon, (nlen, ))
    olat = np.reshape(lat, (nlen, ))
    ovar = np.reshape(out_var, (nlen, ))

    olat = olat[~np.isnan(ovar)]
    olon = olon[~np.isnan(ovar)]
    ovar = ovar[~np.isnan(ovar)]

   #Fill in extrapolated values with nearest neighbor
    out_var = regridder((olon,olat), ovar,
                        (lon,lat), method='nearest')

   #print('out_var.ndim=', out_var.ndim)
   #print('out_var.shape=', out_var.shape)
   #print('out_var.size=', out_var.size)

    return lon, lat, out_var

  def readGridSpecFiles(self):
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

     #print('lonc.ndim=', lonc.ndim)
     #print('lonc.shape=', lonc.shape)
     #print('lonc.size=', lonc.size)

      lonc1d = np.reshape(lonc, ((nx-1)*(ny-1),))
      latc1d = np.reshape(latc, ((nx-1)*(ny-1),))

     #print('len(lonc1d) = ', len(lonc1d))

      lon1d.extend(lonc1d)
      lat1d.extend(latc1d)

     #print('len(lon1d) = ', len(lon1d))

      nc.close()

   #print('len(lon1d) = ', len(lon1d))
    return lat1d, lon1d

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

      arr[np.isnan(arr)] = 0
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

  def get_level(self, data, level=0):
    if(3 == data.ndim):
      nz, ny, nx = data.shape
      var2d = data[level,:,:]
      var1d = np.reshape(var2d, (ny*nx, ))
    else:
      var1d = []
      nt, nz, ny, nx = data.shape
      for it in range(nt):
        var = data[it,level,:,:]
        var = var.reshape((ny*nx))
        var1d.extend(var)

   #print('len(var1d) = ', len(var1d))
    return var1d

  def get_latlon_data(self, varname, nlon=360, nlat=181, method='linear'):
    lat1d, lon1d = self.readGridSpecFiles()
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

    latlon_var = self.interp2latlon_data(lon1d, lat1d, var, nlon=nlon, nlat=nlat, method=method)

    return latlon_var

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

  def interp2latlon_data(self, lon1d, lat1d, var, nlon=360, nlat=181, method='linear'):
    print('var.ndim = ', var.ndim)
    if(2 == var.ndim):
      ny, nx = var.shape

      var1d = np.reshape(var, (ny*nx, ))
      olons,olats,latlon_var = self.interp_to_latlon(lon1d, lat1d, var1d, nlon=nlon, nlat=nlat, method=method)
    else:
      if(3 == var.ndim):
        nz, ny, nx = var.shape
      else:
        nt, nz, ny, nx = var.shape

      latlon_var = np.zeros((nz, int(nlat), int(nlon)), dtype=float)

      for level in range(nz):
        var1d = self.get_level(var, level=level)
        if(self.debug):
          print('processing level ', level)
        olons,olats,ovar = self.interp_to_latlon(lon1d, lat1d, var1d, nlon=nlon, nlat=nlat, method=method)
        latlon_var[level,:,:] = ovar[:,:]

    return latlon_var

#----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1

  datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs/increment/'
 #griddir = '/work/noaa/gsienkf/weihuang/tools/UFS-RNR-tools/JEDI.FV3-increments/grid/C96/'
  griddir = '/work/noaa/gsienkf/weihuang/tools/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
   #datafiles.append('fv3_history.tile%s.nc'%ntile)
    datafile = '%s20210109.000000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)
   #gridspecfiles.append('grid_spec.tile%s.nc'%ntile)
   #gridfile = '%sC96_grid.tile%s.nc' %(griddir, ntile)
    gridfile = '%sC48_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

  rf = RegridFV3(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)
  lat1d, lon1d = rf.readGridSpecFiles()

  varname = 'T'
  var = rf.read3Dvar(varname)

 #print('var.ndim=', var.ndim)
 #print('var.shape=', var.shape)
 #print('var.size=', var.size)

  var1d = rf.get_level(var, level=30)

  print('len(lon1d) = ', len(lon1d))
  print('len(lat1d) = ', len(lat1d))
  print('len(var1d) = ', len(var1d))

  olons,olats,pvar = rf.interp_to_latlon(lon1d, lat1d, var1d, nlon=360, nlat=181, method='linear')

 #make plot on output mesh
  m = Basemap(lon_0=180)
  m.drawcoastlines()
  m.drawmapboundary()
 #m.contourf(olons,olats,pvar,15)
  m.contourf(olons,olats,pvar,5)
  m.colorbar()
  plt.show()

