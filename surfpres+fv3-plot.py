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
from readIODA2Obs import ReadIODA2Obs

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uvOnly = 0
  plotdiff = 0
  addobs = 0
  uselogp = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'uvOnly=',
                             'plotdiff=', 'addobs=', 'uselogp='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uvOnly'):
      uvOnly = int(a)
    elif o in ('--plotdiff'):
      plotdiff = int(a)
    elif o in ('--addobs'):
      addobs = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('plotdiff = ', plotdiff)
  print('uvOnly = ', uvOnly)

#------------------------------------------------------------------------------
 #griddir = '/work/noaa/gsienkf/weihuang/tools/UFS-RNR-tools/JEDI.FV3-increments/grid/C96/'
  griddir = '/work/noaa/gsienkf/weihuang/tools/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

  if(uvOnly):
    datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.uvOnly/increment/'
    if(plotdiff):
      snd_dir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.all/increment/'
  else:
   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.all/increment/'
   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.uvTq/increment/'
   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs/increment/'
   #datadir = '/work/noaa/gsienkf/weihuang/jedi/surface/analysis.getkf.5members.36procs/increment/'
   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.psOnly/increment/'
    datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/amsua/analysis.getkf.5members.36procs/increment/'

  datafiles = []
  snd_files = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s20210109.000000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)
    if(plotdiff):
      snd_file = '%s20210109.000000.fv_core.res.tile%s.nc' %(snd_dir, ntile)
      snd_files.append(snd_file)

  rg = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)
  if(plotdiff):
    rg.setSecondFiles(snd_files)

  nlon = 360
  nlat = nlon/2 + 1
  varname = 'T'
  var = rg.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

 #print('var.ndim = ', var.ndim)
 #print('var.shape = ', var.shape)

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)
  clevs = np.arange(-0.5, 0.51, 0.01)
  cblevs = np.arange(-0.5, 0.6, 0.1)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#------------------------------------------------------------------------------
  if(addobs):
    filename = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/ioda_v2_data/obs/ncdiag.oper.ob.PT6H.sondes.2021-01-08T21:00:00Z.nc4'
    rio = ReadIODA2Obs(debug=debug, filename=filename)
   #lat, lon = rio.get_latlon()
    lat, lon = rio.get_latlon4var(varname='/ObsValue/surface_pressure')

   #print('lat = ', lat)
   #print('lon = ', lon)

    gp.set_obs_latlon(obslat=lat, obslon=lon)

#------------------------------------------------------------------------------
  gp.set_label('Temperature (K)')

  if(uvOnly):
    if(plotdiff):
      image_prefix = 'diff_uvOnly_jedi_sondes'
      title_preix = 'Diff uvOnly JEDI Sondes Temperature at'
    else:
      image_prefix = 'uvOnly'
      title_preix = 'uvOnly JEDI Sondes Temperature at'
  else:
   #image_prefix = 'uvTq_jedi_sondes'
   #title_preix = 'uvTq JEDI Sondes Temperature at
    image_prefix = 'jedi_surface_Pressure'
    title_preix = 'JEDI Surface Pressure (use surface data) Temperature at'

 #levs = [30, 55]
  levs = [1, 10, 23, 30, 52, 62]

  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(image_prefix, lev)
    title = '%s level %d' %(title_preix, lev)
    gp.set_imagename(imgname)
    gp.set_title(title)
    if(addobs):
     #gp.plot(pvar, addmark=1, marker='x', size=3, color='green')
      gp.plot(pvar, addmark=1, marker='x', size=1, color='green')
    else:
      gp.plot(pvar)

 #lons = [100, 115, 125, 140, 175, 190, 225]
  lons = [40, 105, 170, 270, 300]

  for lon in lons:
    pvar = var[:,:,lon]
    imgname = '%s_lon_%d.png' %(image_prefix, lon)
    title = '%s longitude %d' %(title_preix, lon)
    gp.set_imagename(imgname)
    gp.set_title(title)
    if(uselogp):
      gp.plot_meridional_section_logp(pvar)
    else:
      gp.plot_meridional_section(pvar)

 #lats = [-35, -20, 45, 55]
  lats = [-30, 0, 45, 70]
  for lat in lats:
    pvar = var[:,90+lat,:]
    imgname = '%s_lat_%d.png' %(image_prefix, lat)
    title = '%s latitude %d' %(title_preix, lat)
    gp.set_imagename(imgname)
    gp.set_title(title)
    if(uselogp):
      gp.plot_zonal_section_logp(pvar)
    else:
      gp.plot_zonal_section(pvar)

