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
  plotdiff = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'plotdiff='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--plotdiff'):
      plotdiff = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('plotdiff = ', plotdiff)

 #griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C96/'
  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

 #datadir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.5members.36procs/increment/'
  datadir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs/increment/'
  if(plotdiff):
    snd_dir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.all/increment/'

  datafiles = []
  snd_files = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s20200110.030000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)
    if(plotdiff):
      snd_file = '%s20200110.030000.fv_core.res.tile%s.nc' %(snd_dir, ntile)
      snd_files.append(snd_file)

  rg = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)
  if(plotdiff):
    rg.setSecondFiles(snd_files)

  nlon = 360
  nlat = nlon/2 + 1
  varname = 't'
  var = rg.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

  print('var.ndim = ', var.ndim)
  print('var.shape = ', var.shape)
  nz, ny, nx = var.shape

  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)

  gp.set_label('Temperature (K)')

  imgpreix = 'uvTq'

  title_preix = 'JEDI'

 #levs = [30, 55]
  levs = [10, 23, 30, 52, 55, 60, 63]

 #for lev in range(nz):
  for lev in levs:
    pvar = var[lev,:,:]

    print('var %s min %8.5f max %8.5f at level %d' %(varname, np.min(pvar), np.max(pvar), lev))

    if((np.min(pvar) > -1.0e-2) and (np.max(pvar) < 1.0e-2)):
      continue

    if(plotdiff):
      imgname = 'diff_getkf_sondes_lev_%d.png' %(lev)
      title = 'Diff GETKF Sondes Temperature at level %d' %(lev)
    else:
      imgname = 'sondes_lev_%d.png' %(lev)
      title = 'GETKF Sondes Temperature at level %d' %(lev)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot(pvar)

 #lons = [100, 115, 125, 140, 175, 190, 225]
  lons = [40, 105, 170, 270, 300]

  for lon in lons:
    pvar = var[:,:,lon]
    if(plotdiff):
      imgname = 'diff_getkf_sondes_lon_%d.png' %(lon)
      title = 'Diff GETKF Sondes Temperature at longitude %d' %(lon)
    else:
      imgname = 'getkf_sondes_lon_%d.png' %(lon)
      title = 'GETKF Sondes Temperature at longitude %d' %(lon)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot_meridional_section(pvar)

 #lats = [-35, -20, 45, 55]
  lats = [-30, 0, 45, 70]
  for lat in lats:
    pvar = var[:,90+lat,:]
    if(plotdiff):
      imgname = 'diff_getkf_sondes_lat_%d.png' %(lat)
      title = 'Diff GETKF Sondes Temperature at latitude %d' %(lat)
    else:
      imgname = 'getkf_sondes_lat_%d.png' %(lat)
      title = 'GETKF Sondes Temperature at latitude %d' %(lat)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot_zonal_section(pvar)

