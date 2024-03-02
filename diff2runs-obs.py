#=========================================================================
import os
import sys
import types
import getopt
import netCDF4

import numpy as np

from readIODA2Obs import ReadIODA2Obs

#=========================================================================
class Diff2Obs():
  def __init__(self, debug=0):
    self.debug = debug
    self.precision = 1

  def set_precision(self, precision=1):
    self.precision = precision

  def reorder(self, latx, lonx, prsx, laty, lony, prsy, varyin):
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

  def process(self, xf, yf, ndim=1):
    ncxf = ReadIODA2Obs(debug=debug, filename=xf)
    latx, lonx = ncxf.get_latlon()

    ncyf = ReadIODA2Obs(debug=debug, filename=yf)
   #laty, lony = ncyf.get_latlon()

    xgrplist = ncxf.get_grouplist()
    ygrplist = ncyf.get_grouplist()

   #print('xgrplist = ', xgrplist)
   #print('ygrplist = ', ygrplist)

    delt = 1.0e-6

    for xgrp in xgrplist:
      if(xgrp in ygrplist):
       #print('xgrp = ', xgrp)
        xvarlist = ncxf.get_variablelist_in_group(xgrp)
        yvarlist = ncyf.get_variablelist_in_group(xgrp)
       #print('xvarlist = ', xvarlist)
       #print('yvarlist = ', yvarlist)
        for xvar in xvarlist:
          if(xvar in yvarlist):
           #print('xvar = ', xvar)
            if(xvar not in ['stationIdentification', 'variable_names']):
              grpvarname = '/%s/%s' %(xgrp, xvar)
             #print('grp: %s, var: %s' %(xgrp, xvar))
              if(ndim > 1):
                varx = ncxf.get_var_2d(grpvarname)
                vary = ncyf.get_var_2d(grpvarname)
              else:
                varx = ncxf.get_var(grpvarname)
                vary = ncyf.get_var(grpvarname)
              dvar = varx - vary
              dmin = np.min(dvar)
              dmax = np.max(dvar)
             #print('grp: %s, var: %s, min: %e, max: %e' %(xgrp, xvar, dmin, dmax))
              if(dmin < -delt or dmax > delt):
                print('grp: %s, var: %s, min: %e, max: %e' %(xgrp, xvar, dmin, dmax))

#----------------------------------------------------------------------
if __name__ == '__main__':
  debug = 0
  topdir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data'
  datestr = '2020121500'
 #obslist = ['aircraft', 'scatwind', 'sfc']
  obslist = ['sfcship']
 #obslist = ['amsua_n19']

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

  d2f = Diff2Obs(debug=debug)

#-----------------------------------------------------------------------
  for obstype in obslist:
    for n in range(len(xflist)):
      xdatadir = '%s/hofx.%s' %(topdir, xflist[n])
      xfile = '%s/%s_%s_%s_s.nc4' %(xdatadir, obstype, enlist[n], datestr)
      ydatadir = '%s/hofx.%s' %(topdir, yflist[n])
      yfile = '%s/%s_%s_%s_s.nc4' %(ydatadir, obstype, enlist[n], datestr)

      if(os.path.exists(xfile)):
        print('xfile = %s' %(xfile))
      else:
        print('xfile = %s, does not exist. Exit.' %(xfile))
        sys.exit(-1)

      if(os.path.exists(yfile)):
        print('yfile = %s' %(yfile))
      else:
        print('yfile = %s, does not exist. Exit.' %(yfile))
        sys.exit(-1)

      d2f.process(xfile, yfile)

#-----------------------------------------------------------------------

