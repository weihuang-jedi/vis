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
""" TileData """
class TileData:
  """ Constructor """
  def __init__(self, debug=0, output=0, griddir=None, gridtype='C48'):
    """ Initialize class attributes """
    self.debug = debug
    self.output = output
    self.griddir = griddir
    self.gridtype = gridtype

    if(griddi is None):
      print('griddir not defined. Exit.')
      sys.exit(-1)

    self.lon = []
    self.lat = []
    for n in range(6):
      nt = n + 1
      datafile = '%s/%s/%s_oro_data.tile%d.nc' %(griddir, gridtype, gridtype, nt)
      if(os.path.exists(datafile)):
       #print('File No %d: %s' %(nt, datafile))
        self.datalist.append(datafile)

      ncf = netCDF4.Dataset(datafile)
      lon = ncf.variables['geolon'][:,:]
      lat = ncf.variables['geolat'][:,:]
      
     #print('lon.ndim=', lon.ndim)
     #print('lon.shape=', lon.shape)
     #print('lon.size=', lon.size)

      ny, nx = lon.shape

      lonc1d = np.reshape(lon, (nx*ny,))
      latc1d = np.reshape(lat, (nx*ny,))

      self.lon.extend(lonc1d)
      self.lat.extend(latc1d)

      ncf.close()

    print('len(self.lon) = ', len(self.lon))
    print('len(self.lat) = ', len(self.lat))

  def get_var(self, varname):
    data = []
    for datafile in self.datalist:
      ncf = netCDF4.Dataset(datafile)
      var = ncf.variables[varname][:,:]
    
     #print('lon.ndim=', lon.ndim)
     #print('lon.shape=', lon.shape)
     #print('lon.size=', lon.size)

      ny, nx = lon.shape

      lonc1d = np.reshape(lon, (nx*ny,))
      latc1d = np.reshape(lat, (nx*ny,))
      var1d = np.reshape(lat, (nx*ny,))

      self.lon.extend(lonc1d)
      self.lat.extend(latc1d)
      self.var.extend(var1d)
      ncf.close()

    print('len(self.lon) = ', len(self.lon))
    print('len(self.lat) = ', len(self.lat))
    print('len(self.orog) = ', len(self.orog))

    return np.array(self.lat), np.array(self.lon), np.array(self.orog)

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  gridtype = 'C96'
  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'gridtype=', 'griddir='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--griddir'):
      griddir = a
    elif o in ('--gridtype'):
      gridtype = a
    else:
      assert False, 'unhandled option'

  pt = TileData(debug=debug, griddir=griddir, gridtype=gridtype)
  pt.process()
  lat, lon, orog = dt.get_data()

  gp = GeneratePlot(debug=debug, output=output)
  gp.set_grid(lat, lon)

  imgname = 'ps.png'
  gp.set_imagename(imgname)
  title = 'Surface Pressure (Unit: hPa)'
  gp.set_title(title)
  gp.plot(ps)

  imgname = 'orog.png'
  gp.set_imagename(imgname)
  title = 'Orograph (Unit: m)'
  gp.set_title(title)
  gp.plot(orog)

