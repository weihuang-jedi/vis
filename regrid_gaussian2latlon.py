import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata as regridder

class RegridGaussian2LatLon():
  def __init__(self, debug=0, nlon=360, nlat=181, method='linear'):
    self.debug = debug
    self.nlon = nlon
    self.nlat = nlat
    self.method = method

    if(self.debug):
      print('debug = ', debug)

    dlon = 360.0/nlon
    dlat = 180.0/(nlat - 1)
   #Create a lat-lon uniform grid
    out_lon = np.arange(0.0, 360.0, dlon)
    out_lat = np.arange(-90.0, 91.0, dlat)
    self.lon, self.lat = np.meshgrid(out_lon, out_lat)

   #print('self.lon.size = ', self.lon.size)
   #print('self.lon.shape = ', self.lon.shape)

  def get_latlon(self):
    return self.lat, self.lon

  def interp_to_latlon(self, lon_1d, lat_1d, var_1d):
    '''
    Interpolate a variable on cube-sphere grid (such as FV3) to LatLon grid
    '''
   #print('var_1d.ndim = ', var_1d.ndim)
   #print('var_1d.size = ', var_1d.size)
   #print('var_1d.shape = ', var_1d.shape)

    # Interpolate from cube to lat-lon grid
    out_var = regridder((lon_1d,lat_1d), var_1d,
                        (self.lon, self.lat), method=self.method)

    nlen = int(self.nlon*self.nlat)
   #print('nlen = ', nlen)
   #print('self.lon.size = ', self.lon.size)
    olon = np.reshape(self.lon, (nlen, ))
    olat = np.reshape(self.lat, (nlen, ))
    ovar = np.reshape(out_var, (nlen, ))

   #print('olon.size = ', olon.size)
   #print('olon.shape = ', olon.shape)

   #print('olat.size = ', olat.size)
   #print('olat.shape = ', olat.shape)

    olat = olat[~np.isnan(ovar)]
    olon = olon[~np.isnan(ovar)]
    ovar = ovar[~np.isnan(ovar)]

   #Fill in extrapolated values with nearest neighbor
    out_var = regridder((olon,olat), ovar,
                        (self.lon,self.lat), method='nearest')

   #print('out_var.ndim=', out_var.ndim)
   #print('out_var.shape=', out_var.shape)
   #print('out_var.size=', out_var.size)

    return out_var

  def read3Dvar(self,filename,varname,ntime=0):
    ncfile = netCDF4.Dataset(filename, 'r')
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
    data = ncfile.variables[varname][ntime,:,:,:]
    ncfile.close()

    mlat = lat.size
    mlon = lon.size

    lat2d = np.zeros((mlat, mlon))
    lon2d = np.zeros((mlat, mlon))

    for i in range(mlon):
      lat2d[:,i] = lat[:]
    for j in range(mlat):
      lon2d[j,:] = lon[:]

    lat1d = lat2d.flatten()
    lon1d = lon2d.flatten()

    return lat1d, lon1d, data

  def get_level(self, data, level=0):
    nz, ny, nx = data.shape
    var2d = data[level,:,:]
    var1d = np.reshape(var2d, (ny*nx, ))

   #print('len(var1d) = ', len(var1d))
    return var1d

  def interp2latlon(self, lon1d, lat1d, var):
    nz, ny, nx = var.shape

   #print('var.shape=', var.shape)

    latlon_var = np.zeros((nz, int(self.nlat), int(self.nlon)), dtype=float)

    for level in range(nz):
      var1d = self.get_level(var, level=level)
      if(self.debug):
        print('processing level ', level)
      ovar = self.interp_to_latlon(lon1d, lat1d, var1d)
      latlon_var[level,:,:] = ovar[:,:]

    return latlon_var

#----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1

  datadir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.40t1n_36p/analysis/increment'
  filename = 'xainc.20200110_030000z.nc4'

  fullname = '%s/%s' %(datadir, filename)

#------------------------------------------------------------------------------
  nlon = 360
  nlat = nlon/2 + 1

  rg = RegridGaussian2LatLon(debug=debug, nlon=nlon, nlat=nlat, method='linear')

  varname = 'T'
  lat1d, lon1d, invar = rg.read3Dvar(fullname,varname)

  var = rg.interp2latlon(lon1d, lat1d, invar)

  print('var.ndim=', var.ndim)
  print('var.shape=', var.shape)
  print('var.size=', var.size)

 #pvar = rg.get_level(var, level=30)
  pvar = var[50,:,:]
  lat,lon = rg.get_latlon()

 #make plot on output mesh
  m = Basemap(lon_0=180)
  m.drawcoastlines()
  m.drawmapboundary()
 #m.contourf(lon,lat,pvar,15)
  m.contourf(lon,lat,pvar,5)
  m.colorbar()
  plt.show()

