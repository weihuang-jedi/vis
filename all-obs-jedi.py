#=========================================================================
import os
import sys
import types
import getopt
import netCDF4
import matplotlib

import numpy as np
import matplotlib.pyplot

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from genplot import GeneratePlot as genplot
from readIODA2Obs import ReadIODA2Obs

#=========================================================================
class PlotObsOverlayIncrement():
  def __init__(self, debug=0, output=0, filename=None):
    self.debug = debug
    self.output = output
    self.filename = filename

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.filename = ', self.filename)

  def get_var(self, varname):
   #print('varname =', varname)

    ncfile = netCDF4.Dataset(self.filename, 'r')
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
   #var = ncfile.variables[varname][itime, :, :, :]
    var = ncfile.variables[varname][:, :, :]
    ncfile.close()

    var[np.isnan(var)] = 0.0

    if(self.debug):
      msg = ('var range for variable %s: (%s, %s).' % (varname, var.min(), var.max()))
      print(msg)

    return lon, lat, var

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uselogp = 1
  varname = 'T'
  casename = 'surf'
  gridfile = '/work2/noaa/gsienkf/weihuang/jedi/case_study/vis/regrid/fv3latlon.nc'

  obsdir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/ioda_v2_data'
  filelist = ['aircraft_tsen_obs_2020011006.nc4', 'amsua_n19_obs_2020121500_m.nc4',
              'iasi_metop-b_obs_2020121500_m.nc4', 'satwind_obs_2020011006.nc4',
              'scatwind_obs_2020011006.nc4', 'sfc_ps_obs_2020011006.nc4',
              'sfcship_tsen_obs_2020011006.nc4', 'sondes_tsen_obs_2020011006.nc4',
              'vadwind_obs_2020011006.nc4', 'windprof_obs_2020011006.nc4']

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'uselogp=',
                             'gridfile=', 'obsfile=', 'casename=', 'varname='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
    elif o in ('--casename'):
      casename = a
    elif o in ('--varname'):
      varname = a
    elif o in ('--gridfile'):
      gridfile = a
    elif o in ('--obsfile'):
      obsfile = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

#------------------------------------------------------------------------------
  po = PlotObsOverlayIncrement(debug=debug, output=output, filename=gridfile)
  lon1d, lat1d, var = po.get_var(varname)

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat1d, lon=lon1d)
  robs = ReadIODA2Obs(debug=debug, filename=None)
  obslist = []
  for file in filelist:
    filename = '%s/%s' %(obsdir, file)
    obslat, obslon = robs.get_latlon_from_file(filename)
    obs = [obslon, obslat]
    obslist.append(obs)

#------------------------------------------------------------------------------
  gp.set_label('Temperature (K)')

  imageprefix = '%s_%s_increment' %(casename, varname)
  titleprefix = '%s %s Increment at' %(casename, varname)

#------------------------------------------------------------------------------
 #clevs = np.arange(-0.5, 0.51, 0.01)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
 #clevs = np.arange(-0.2, 0.21, 0.01)
 #cblevs = np.arange(-0.2, 0.3, 0.1)
 #clevs = np.arange(-2., 2.1, 0.1)
 #cblevs = np.arange(-2., 3., 1.)
  clevs = np.arange(-5.0, 5.1, 0.1)
  cblevs = np.arange(-5.0, 6.0, 1.0)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

  levs = [30, 40, 50, 60]

  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(imageprefix, lev)
    title = '%s level %d' %(titleprefix, lev)
    gp.set_imagename(imgname)
    gp.set_title(title)
    gp.plot_with_obs(pvar, obslist)

  sys.exit(-1)

#------------------------------------------------------------------------------
  lons = [40, 105, 170, 270, 300]

  for lon in lons:
    pvar = var[:,:,lon]
    title = '%s longitude %d' %(titleprefix, lon)
    gp.set_title(title)

    imgname = '%s_lon_%d_logp.png' %(imageprefix, lon)
    gp.set_imagename(imgname)
    gp.plot_meridional_section_logp(pvar)

    imgname = '%s_lon_%d_level.png' %(imageprefix, lon)
    gp.set_imagename(imgname)
    gp.plot_meridional_section(pvar)

#------------------------------------------------------------------------------
  lats = [-30, 0, 45, 70]

  for lat in lats:
    pvar = var[:,90+lat,:]
    gp.set_title(title)
    title = '%s latitude %d' %(titleprefix, lat)

    imgname = '%s_lat_%d_logp.png' %(imageprefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section_logp(pvar)

    imgname = '%s_lat_%d_level.png' %(imageprefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section(pvar)

