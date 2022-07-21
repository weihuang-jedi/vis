#=========================================================================
import os
import sys
import types
import getopt

import numpy as np
import matplotlib
import matplotlib.pyplot

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from genplot import GeneratePlot as genplot
from scipy_regridder import RegridFV3 as regridder

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  prefix = 'ori'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'prefix='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--prefix'):
      prefix = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('prefix = ', output)

 #griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C96/'
  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

 #dirname = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/manual-obs'
 #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/manual-obs/analysis.getkf.80members.36procs.ori/increment/'
 #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/manual-obs/analysis.getkf.80members.36procs.new/increment/'
  dirname = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/surfpres'
 #dirname = '/work/noaa/gsienkf/weihuang/jedi/surf'
  datadir = '%s/analysis.getkf.80members.36procs.%s/increment/' %(dirname, prefix)

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s20210109.000000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)

  rg = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)

  nlon = 360
  nlat = nlon/2 + 1
  varname = 'T'
  var = rg.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

  print('var.ndim = ', var.ndim)
  print('var.shape = ', var.shape)

 #------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)

  gp.set_label('Temperature (K)')

 #clevs = np.arange(-0.01, 0.011, 0.001)
 #cblevs = np.arange(-0.01, 0.015, 0.005)

  clevs = np.arange(-0.5, 0.51, 0.01)
  cblevs = np.arange(-0.5, 0.6, 0.1)
  gp.set_precision(precision=1)

  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

  image_prefix = '%s_JEDI_sondes_PSonly_singleobs' %(prefix)
  title_prefix = '%s JEDI sondes PSonly Single Obs' %(prefix)

 #------------------------------------------------------------------------------
  lons = [170]
  for lon in lons:
    pvar = var[:,:,lon]
    title = '%s Temperature at longitude %d' %(title_prefix, lon)
    imgname = '%s_logp_lon_%d.png' %(image_prefix, lon)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot_meridional_section_logp(pvar)
    imgname = '%s_level_lon_%d.png' %(image_prefix, lon)
    gp.set_imagename(imgname)
    gp.plot_meridional_section(pvar)

 #------------------------------------------------------------------------------
 #levs = [0, 10, 13, 16, 24, 29, 36, 42, 52, 63]
  levs = [0, 1, 2, 61, 62, 63]
  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(image_prefix, lev)
    title = '%s Temperature at level %d' %(title_prefix, lev)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot(pvar)

 #------------------------------------------------------------------------------
  lats = [-41, -22, 0, 22, 41]
  for lat in lats:
    pvar = var[:,90+lat,:]
    title = '%s Temperature at latitude %d' %(title_prefix, lat)
    imgname = '%s_logp_lat_%d.png' %(image_prefix, lat)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot_zonal_section_logp(pvar)
    imgname = '%s_level_lat_%d.png' %(image_prefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section(pvar)

