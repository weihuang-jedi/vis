#########################################################################
#$Id: bld.py 28 2021-01-21 15:10:31Z whuang $
#$Revision: 28 $
#$HeadURL: file:///Users/whuang/.wei_svn_repository/trunk/jedi-build-tools/bld.py $
#$Date: 2021-01-21 08:10:31 -0700 (Thu, 21 Jan 2021) $
#$Author: whuang $
#########################################################################

import getopt
import os, sys
import types
import time
import datetime
import subprocess
import netCDF4

import matplotlib.pyplot as plt
import numpy as np

import matplotlib
import matplotlib.pyplot

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

def cmdout(command):
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
    ostr = result.stdout
    return ostr.strip()

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

    self.precision = 1

  def set_precision(self, precision=1):
    self.precision = precision

  def set_grid(self, lat, lon):
    self.lat = np.array(lat)
    self.lon = np.array(lon)

  def set_default(self):
    self.image_name = 'sample.png'

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
   #self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
    self.cmapname = 'jet'

    self.obslat = []
    self.obslon = []

   #self.clevs = np.arange(-0.2, 0.21, 0.01)
   #self.cblevs = np.arange(-0.2, 0.3, 0.1)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Time (sec)'
    self.title = 'Time (sec)'

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_imagename(self, imagename):
    self.image_name = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

  def set_label(self, label):
    self.label = label

  def set_title(self, title):
    self.title = title

  def build_basemap(self):
    basemap_dict = {'resolution': 'c', 'projection': 'cyl',
                    'llcrnrlat': -90.0, 'llcrnrlon': 0.0,
                    'urcrnrlat':  90.0, 'urcrnrlon': 360.0}
    basemap_dict['lat_0'] = 0.0
    basemap_dict['lon_0'] = 180.0

    basemap = Basemap(**basemap_dict)

    return basemap

  def create_image(self, plt_obj, savename):
    msg = ('Saving image as %s.' % savename)
    print(msg)
    kwargs = {'transparent': True, 'dpi': 500}
    plt_obj.savefig(savename, **kwargs)

  def display(self, output=False, image_name=None):
    if(output):
      if(image_name is None):
        image_name=self.image_name
      self.plt.tight_layout()
      kwargs = {'plt_obj': self.plt, 'savename': image_name}
      self.create_image(**kwargs)
    else:
      self.plt.show()

  def plot(self, pvar):
    self.basemap = self.build_basemap()

    self.plt = matplotlib.pyplot
    try:
      self.plt.close('all')
      self.plt.clf()
    except Exception:
      pass

    self.fig = self.plt.figure()
    self.ax = self.plt.subplot()

    msg = ('plot variable min: %s, max: %s' % (np.min(pvar), np.max(pvar)))
    print(msg)

    (self.x, self.y) = self.basemap(self.lon, self.lat)

    v1d = np.array(pvar)

    print('self.x.shape = ', self.x.shape)
    print('self.y.shape = ', self.y.shape)
    print('v1d.shape = ', v1d.shape)

   #contfill = self.basemap.contourf(self.x, self.y, v1d, tri=True,
   #                                 levels=self.clevs, extend=self.extend,
   #                                 alpha=self.alpha, cmap=self.cmapname)
    contfill = self.plt.tricontourf(self.x, self.y, v1d,
                                    alpha=self.alpha, cmap=self.cmapname)

    cb = self.fig.colorbar(contfill, orientation=self.orientation,
                           pad=self.pad)

    cb.set_label(label=self.label, size=self.size, weight=self.weight)

    cb.ax.tick_params(labelsize=self.labelsize)
   #if(self.precision == 0):
   #  cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in self.cblevs], minor=False)
   #elif(self.precision == 1):
   #  cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in self.cblevs], minor=False)
   #elif(self.precision == 2):
   #  cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in self.cblevs], minor=False)
   #else:
   #  cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in self.cblevs], minor=False)

    self.ax.set_title(self.title)

    self.plot_coast_lat_lon_line()

    self.display(output=self.output, image_name=self.image_name)

  def plot_coast_lat_lon_line(self):
   #https://matplotlib.org/basemap/users/geography.html
   #map.drawmapboundary(fill_color='aqua')
   #map.fillcontinents(color='#cc9955', lake_color='aqua')
   #map.drawcounties()
   #map.drawstates(color='0.5')

   #draw coastlines
    color = 'black'
    linewidth = 0.5
    self.basemap.drawcoastlines(color=color, linewidth=linewidth)

   #draw parallels
    color = 'green'
    linewidth = 0.5
    fontsize = 8
    dashes = [10, 10]
    circles = np.arange(-90,90,30)
    self.basemap.drawparallels(np.arange(-90,90,30),labels=[1,1,0,1],
                               color=color, linewidth=linewidth,
                               dashes=dashes, fontsize=fontsize)

   #draw meridians
    color = 'green'
    linewidth = 0.5
    fontsize = 8
    dashes = [10, 10]
    meridians = np.arange(0,360,30)
    self.basemap.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1],
                               color=color, linewidth=linewidth,
                               dashes=dashes, fontsize=fontsize)

#------------------------------------------------------------------
def get_grid_latlon(gridsize):
  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/%s' %(gridsize)
  lon1d = []
  lat1d = []

  for ntile in range(1,7,1):
    gridspecfile = '%s/%s_grid.tile%s.nc' %(griddir, gridsize, ntile)
    print('reading ',gridspecfile)

    if(not os.path.exists(gridspecfile)):
      print('filename: %s does not exist, stop' %(gridspecfile))
      sys.exit(-1)

    ncfl = netCDF4.Dataset(gridspecfile)
    lons = ncfl.variables['x'][:]
    lats = ncfl.variables['y'][:]

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

    ncfl.close()

  return lon1d, lat1d

#------------------------------------------------------------------
def get_fv3_var(datadir, varname, prefix=None, nt=0):
  var = []

  for ntile in range(1,7,1):
    if(prefix is None):
      filename = '%s/fv_core.res.tile%d.nc' %(datadir, ntile)
    else:
      filename = '%s/%s.fv_core.res.tile%d.nc' %(datadir, prefix, ntile)

    print('filename: ', filename)
    if(not os.path.exists(filename)):
      print('filename: %s does not exist, stop' %(filename))
      sys.exit(-1)

    ncfl = netCDF4.Dataset(filename)
    val = ncfl.variables[varname][nt,:,:,:]
    ncfl.close()

    var.append(val)

  print('len(var) = ', len(var))
  print('var[0].shape = ', var[0].shape)

  return var

#------------------------------------------------------------------
def get_var1d_at_level(var, lvl):
  ntile = 6
  nlev, nlat, nlon = var[0].shape

  print('len(var) = ', len(var))
  print('var[0].shape = ', var[0].shape)

  var1d = []
  for n in range(ntile):
   #v1d = var[n][lvl,:,:].flatten()
    v1d = np.reshape(var[n][lvl,:,:], (nlat*nlon,))
    var1d.extend(v1d)

  print('nlat*nlon = ', nlat*nlon)
  print('len(var1d) = ', len(var1d))

  arr1d = np.array(var1d)

  return arr1d
  
#------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
 #casename = 'sondes'
  casename = 'aircraft'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'casename='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--casename'):
      casename = a
    else:
      print('option: ', a)
      assert False, 'unhandled option'

  gp = GeneratePlot(debug=debug, output=output)

  lon1d, lat1d = get_grid_latlon('C48')
  gp.set_grid(lat1d, lon1d)

  basedir = '/work2/noaa/gsienkf/weihuang/jedi/case_study'
  casedir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/develop_code'
  runcase = 'run_80.40t1n_36p'
  dir1 = '%s/%s/%s/analysis/increment' %(basedir, casename, runcase)
  dir2 = '%s/%s/%s/analysis/increment' %(casedir, casename, runcase)

  varname = 'T'

  basevar = get_fv3_var(dir1, varname, prefix='20200110.030000', nt=0)
  casevar = get_fv3_var(dir2, varname, prefix='20200110.030000', nt=0)

  ntile = 6
  nlev, nlat, nlon = basevar[0].shape

  print('len(basevar) = ', len(basevar))
  print('basevar[0].shape = ', basevar[0].shape)

  levs = [0, 10, 20, 30, 40, 50, 60]
  for lvl in levs:
    basevar1d = get_var1d_at_level(basevar, lvl)
    casevar1d = get_var1d_at_level(casevar, lvl)

    val1d = casevar1d - basevar1d

    print('%s diff at lvl: %d, min: %f, max: %f' %(varname, lvl, np.min(val1d), np.max(val1d)))

    name = '%s-diff_%s_level_%d.png' %(casename, varname, lvl)
    gp.set_imagename(name)
    title = '%s diff %s at level %d' %(casename, varname, lvl)
    gp.set_title(title)
    gp.plot(val1d)

