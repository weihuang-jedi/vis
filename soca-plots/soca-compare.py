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
if __name__== '__main__':
  debug = 1
  output = 0

  gridfile = '../regrid/grids/ocn_2014_01.nc'
  ncg = netCDF4.Dataset(gridfile, 'r')
  lon = ncg.variables['geolon'][:,:]
  lat = ncg.variables['geolat'][:,:]
  ncg.close()

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    else:
      print('option: ', a)
      assert False, 'unhandled option'

  gp = GeneratePlot(debug=debug, output=output)

  lat1d = lat.flatten()
  lon1d = lon.flatten()
  lon1d = np.where(lon1d > 0, lon1d, lon1d+360.0)
  gp.set_grid(lat1d, lon1d)

 #dir1 = '/work2/noaa/gsienkf/weihuang/jedi/run.soca/soca_solver.40t2n_80p'
 #dir2 = '/work2/noaa/gsienkf/weihuang/ufs/soca/new-soca-solver/soca_solver.40t2n_80p'
  dir1 = '/work2/noaa/gsienkf/weihuang/jedi/run.soca/soca_solver.36t2n_72p'
 #dir2 = '/work2/noaa/gsienkf/weihuang/jedi/run.soca/solver.36t2n_72p'
  dir2 = '/work2/noaa/gsienkf/weihuang/jedi/run.soca/solver.36t2n_72p-soca_observer_output'

  varname = 'Temp'
 #varname = 'Salt'

  file1 = '%s/ocn.LETKF.an.2015-12-01T12:00:00Z.nc' %(dir1)
  if(os.path.exists(file1)):
    pass
  else:
    print('file1: %s does not exist, stop' %(file1))
   #continue
    sys.exit(-1)

  file2 = '%s/ocn.LETKF.an.2015-12-01T12:00:00Z.nc' %(dir2)
  if(os.path.exists(file2)):
    pass
  else:
    print('file2: %s does not exist, stop' %(file2))
   #continue
    sys.exit(-1)

  print('file1: ', file1)
  print('file2: ', file2)

  nc1 = netCDF4.Dataset(file1, 'r')
  nc2 = netCDF4.Dataset(file2, 'r')
  temp1 = nc1.variables[varname][0,:,:,:]
  temp2 = nc2.variables[varname][0,:,:,:]
  nc1.close()
  nc2.close()

  difftemp = temp2 - temp1

  print('%s diff min: %f, max: %f' %(varname, np.min(difftemp), np.max(difftemp)))

  nlay, nlat, nlon = temp1.shape
 #print('temp1.shape: ', temp1.shape)

  plotlvls = [0, 2]
 #case = 'solver-full-letkf'
  case = 'solver-onepass'
  for n in range(nlay):
    val = temp2[n,:,:] - temp1[n,:,:]
    val1d = val.flatten()
   #print('val.shape: ', val.shape)
    print('%s diff at level %d min: %f, max: %f' %(varname, n, np.min(val1d), np.max(val1d)))
    if n in plotlvls:
      name = '%s-diff_%s_level_%d.png' %(case, varname, n)
      gp.set_imagename(name)
      title = '%s diff %s at level %d' %(case, varname, n)
      gp.set_title(title)
      gp.plot(val1d)

