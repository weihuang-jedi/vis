#=========================================================================
import os
import sys
import types
import getopt

import numpy as np
import matplotlib
import matplotlib.pyplot

import matplotlib.cm as cm

#from matplotlib.colors import Normalize
#from matplotlib.ticker import MultipleLocator

from matplotlib import colors
from matplotlib.ticker import PercentFormatter

import netCDF4

from readIODA2Obs import ReadIODA2Obs

#=========================================================================
class ScatterPlotsFor2Runs():
  def __init__(self, debug=0):
    self.debug = debug

    self.set_default()

    self.precision = 1

  def set_precision(self, precision=1):
    self.precision = precision

  def set_default(self):
    self.image_name = 'sample.png'

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
    self.cmapname = 'bwr'
   #self.cmapname = 'rainbow'
   #self.cmapname = 'seismic'

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.xlabel = 'MTS-JEDI'
    self.ylabel = 'MTS-GSI'

    self.title = '%s vs %s Scatter Plot' %(self.xlabel, self.ylabel)

  def set_imagename(self, imagename):
    self.image_name = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

  def set_label(self, xlabel, ylabel, varname, datestr=None, obstype=None):
    self.xlabel = xlabel
    self.ylabel = ylabel

    self.title = '%s vs %s: %s Scatter Plot' %(xlabel, ylabel, varname)
    imagename = '%s_%s_%s' %(xlabel, ylabel, varname)
    if(datestr is not None):
      self.title = '%s, %s' %(self.title, datestr)
      imagename = '%s_%s' %(imagename, datestr)
    if(obstype is not None):
      self.title = '%s, %s' %(self.title, obstype)
      imagename = '%s_%s' %(imagename, obstype)
    self.set_imagename(imagename)

  def set_title(self, title):
    self.title = title

  def create_image(self, plt_obj, savename):
    msg = ('Saving image as %s.' % savename)
    print(msg)
    kwargs = {'transparent': True, 'dpi': 500}
    plt_obj.savefig(savename, **kwargs)

  def display(self, image_name=None):
    if(image_name is None):
      image_name=self.image_name
    self.plt.tight_layout()
    kwargs = {'plt_obj': self.plt, 'savename': image_name}
    self.create_image(**kwargs)

    if(self.debug):
      self.plt.show()

  def scatter_plot(self, x, y, varname):
    self.plt = matplotlib.pyplot
    try:
      self.plt.close('all')
      self.plt.clf()
    except Exception:
      pass

    self.fig = self.plt.figure()
    self.ax = self.plt.subplot()

    scatterplot = self.plt.scatter(x, y, c='blue', s=2, alpha=self.alpha)

    self.ax.set_title(self.title)

    if(varname == 'surface_pressure'):
      self.plt.xlim((-7, 7))
      self.plt.ylim((-7, 7))
    elif(varname == 'airTemperature' or varname == 'virtualTemperature'):
      self.plt.xlim((-10, 10))
      self.plt.ylim((-10, 10))
    elif(varname == 'eastward_wind' or varname == 'northward_wind'):
      self.plt.xlim((-10, 10))
      self.plt.ylim((-10, 10))
    elif(varname == 'specific_humidity'):
      self.plt.xlim((-5, 5))
      self.plt.ylim((-5, 5))
    else:
      self.plt.xlim((-10, 10))
      self.plt.ylim((-10, 10))

    self.plt.xlabel(self.xlabel, fontsize=14)
    self.plt.ylabel(self.ylabel, fontsize=14)
    self.plt.grid(True)

    ax = self.plt.gca()
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    self.display(image_name=self.image_name)

def reorder(latx, lonx, prsx, laty, lony, prsy, varyin):
  varyout = varyin[:]
  nv = len(latx)
  idx = [i for i in range(nv)]
  dlt = 0.0000001
  for n in range(nv):
    k = -1
    i = 0
    while (i < len(idx)):
      if(abs(latx[n]-laty[i]) < dlt):
        if(abs(lonx[n]-lony[i]) < dlt):
          if(abs(prsx[n]-prsy[i]) < dlt):
            varyout[n] = varyin[i]
            k = i
            i = len(idx)
      i += 1;
    if(k >= 0):
      del idx[k]
  return varyout

#----------------------------------------------------------------------
if __name__ == '__main__':
  debug = 0
  topdir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data'
 #obstype = 'aircraft'
 #obstype = 'scatwind'
  obstype = 'sfcship'
  datestr = '2020121500'

  grplist = ['GsiHofX','hofx0','hofx_y_mean_xb0','hofx_y_mean_xb1']
 #grplist = ['hofx0_1','hofx0_2','hofx0_3','hofx0_4','hofx0_5','hofx0_6','hofx0_7','hofx0_8','hofx0_9','hofx0_10']
 #grplist = ['hofxm0_1_1','hofxm0_1_2','hofxm0_1_3','hofxm0_1_4','hofxm0_1_5','hofxm0_1_6','hofxm0_1_7','hofxm0_1_8','hofxm0_1_9','hofxm0_1_10']
 #varname = 'airTemperature'
 #varname = 'windEastward'
  varname = 'stationPressure'
  obsname = '/ObsValue/%s' %(varname)

#-----------------------------------------------------------------------
  xflist = ['nl', 'nl.getkf']
  yflist = ['lo', 'lo.getkf']
  enlist = ['letkf-gfs', 'lgetkf-geos']

  xslist = ['NonLinear', 'NonLinearGETKF']
  yslist = ['LinObs', 'LinObsGETKF']

#-----------------------------------------------------------------------
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'xf=', 'yf='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--xf'):
      xf = a
    elif o in ('--yf'):
      yf = a
    else:
      assert False, 'unhandled option'

#-----------------------------------------------------------------------
  print('debug = ', debug)
 #print('xflist = ', xflist)
 #print('yflist = ', yflist)
 #print('xslist = ', xslist)
 #print('yslist = ', yslist)

  sp42r = ScatterPlotsFor2Runs(debug=debug)

  for n in range(len(xflist)):
    if(0 == n):
      continue

    xdatadir = '%s/hofx.%s' %(topdir, xflist[n])
    xf = '%s/%s_%s_%s_s.nc4' %(xdatadir, obstype, enlist[n], datestr)
    ydatadir = '%s/hofx.%s' %(topdir, yflist[n])
    yf = '%s/%s_%s_%s_s.nc4' %(ydatadir, obstype, enlist[n], datestr)

    if(os.path.exists(xf)):
      print('xf = %s' %(xf))
    else:
      print('xf = %s, does not exist. Exit.' %(xf))

    if(os.path.exists(yf)):
      print('yf = %s' %(yf))
    else:
      print('yf = %s, does not exist. Exit.' %(yf))

    xlabel = xslist[n]
    ylabel = yslist[n]

    ncxf = ReadIODA2Obs(debug=debug, filename=xf)
    latx, lonx = ncxf.get_latlon()
    prsx = ncxf.get_var('/MetaData/pressure')
    obsx = ncxf.get_var(obsname)

    ncyf = ReadIODA2Obs(debug=debug, filename=yf)
    laty, lony = ncyf.get_latlon()
    prsy = ncyf.get_var('/MetaData/pressure')
   #hgty = ncyf.get_var('/MetaData/heightOfSurface')
    obsy = ncyf.get_var(obsname)

    for grpname in grplist:
      grpvarname = '/%s/%s' %(grpname, varname)
      basx = ncxf.get_var(grpvarname)
      basy = ncyf.get_var(grpvarname)

      varx = obsx - basx
      vary = obsy - basy

     #varx = 0.01*(obsx - basx)
     #vary = 0.01*(obsy - basy)

     #lenx = len(latx)
     #leny = len(laty)

     #print('lenx = %d, leny = %d' %(lenx, leny))
      print('grp: %s var: %s min: %f max: %f' %(grpname, varname, np.min(varx), np.max(varx)))

     #newvary = reorder(latx, lonx, prsx, laty, lony, prsy, vary)
      labelname = '_%s_%s' %(grpname, varname)
      sp42r.set_label(xlabel, ylabel, labelname, datestr=datestr, obstype=obstype)
      sp42r.scatter_plot(varx, vary, labelname)

