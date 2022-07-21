#=========================================================================
import os
import sys
import time
import types
import getopt
import netCDF4
import matplotlib

import numpy as np
import matplotlib.pyplot

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from genplot import GeneratePlot as genplot
from scipy_regridder import RegridFV3 as regridder
from readIODA2Obs import ReadIODA2Obs
from readGSIobs import ReadGSIobs

#=========================================================================
class PlotGaussian():
  def __init__(self, debug=0, output=0, bkg=None, anl=None):
    self.debug = debug
    self.output = output
    self.bkg = bkg
    self.anl = anl

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.bkg = ', self.bkg)
      print('self.anl = ', self.anl)

    self.has_snd_file = 0
    self.snd_file = None

  def set_snd_file(self, filename):
    self.has_snd_file = 1
    self.snd_file = filename

  def get_vardims(self, filename, varname):
    if(self.debug > 1):
      print('varname = ', varname)
      print('filename = ', filename)
    ncfile = netCDF4.Dataset(filename, 'r')
    dimids = ncfile.variables[varname].dimensions

    if(self.debug > 1):
      print('dimids = ', dimids)

    self.nx = 0
    self.ny = 0
    self.nz = 0
    self.ntime = 0

    print('dimids = ', dimids)

    for dimname in dimids:
      if(self.debug > 1):
        print('dimname:', dimname)
      if(dimname == 'time'):
        self.ntime = len(ncfile.dimensions[dimname])
      elif(dimname == 'pfull'):
        self.nz = len(ncfile.dimensions[dimname])
      elif(dimname == 'grid_yt'):
        self.ny = len(ncfile.dimensions[dimname])
      elif(dimname == 'grid_xt'):
        self.nx = len(ncfile.dimensions[dimname])

    if(self.debug):
      print('self.nx = ', self.nx)
      print('self.ny = ', self.ny)
      print('self.nz = ', self.nz)
      print('self.ntime = ', self.ntime)

  def get_var(self, varname):
    print('varname =', varname)

    self.get_vardims(self.bkg, varname)

    lat = np.zeros((self.ny, self.nx))
    lon = np.zeros((self.ny, self.nx))

    anl = np.zeros((self.nz, self.ny, self.nx))
    bkg = np.zeros((self.nz, self.ny, self.nx))

    anl_file = netCDF4.Dataset(self.anl, 'r')
    lat = anl_file.variables['lat'][:, :]
    lon = anl_file.variables['lon'][:, :]
    anl = anl_file.variables[varname][0, :, :, :]
    anl_file.close()

    self.lats = lat.flatten()
    self.lons = lon.flatten()

    bkg_file = netCDF4.Dataset(self.bkg, 'r')
    bkg = bkg_file.variables[varname][0, :, :, :]
    bkg_file.close()

    if(self.debug):
      msg = ('analy range for variable %s: (%s, %s).' % (varname, anl.min(), anl.max()))
      print(msg)
      msg = ('bkgrd range for variable %s: (%s, %s).' % (varname, bkg.min(), bkg.max()))
      print(msg)

    incr = anl - bkg

   #print('incr = ', incr)

    return self.lons, self.lats, incr

  def get_diff(self, varname):
    print('varname =', varname)

    self.get_vardims(self.anl, varname)

    lat = np.zeros((self.ny, self.nx))
    lon = np.zeros((self.ny, self.nx))

    fst = np.zeros((self.nz, self.ny, self.nx))
    snd = np.zeros((self.nz, self.ny, self.nx))

    fst_file = netCDF4.Dataset(self.anl, 'r')
    lat = fst_file.variables['lat'][:, :]
    lon = fst_file.variables['lon'][:, :]
    fst = fst_file.variables[varname][0, :, :, :]
    fst_file.close()

    self.lats = lat.flatten()
    self.lons = lon.flatten()

    snd_file = netCDF4.Dataset(self.snd_file, 'r')
    snd = snd_file.variables[varname][0, :, :, :]
    snd_file.close()

    if(self.debug):
      msg = ('fst range for variable %s: (%s, %s).' % (varname, fst.min(), fst.max()))
      print(msg)
      msg = ('snd range for variable %s: (%s, %s).' % (varname, snd.min(), snd.max()))
      print(msg)

    diff = snd - fst

   #print('incr = ', incr)

    return self.lons, self.lats, diff


#=========================================================================
class SplitObs():
  def __init__(self, debug=0, jedifile=None, gsifile=None, varname=None):
    self.debug = debug

    rio = ReadIODA2Obs(debug=debug, filename=jedifile)
    if(varname is None):
      self.jedilat, self.jedilon = rio.get_latlon()
    else:
      self.jedilat, self.jedilon = rio.get_latlon4var(varname=varname)

    rgo = ReadGSIobs(debug=debug, filename=gsifile)
    self.gsilat, self.gsilon = rgo.get_latlon()

    self.delt = 0.000001

    if(self.debug):
      print('debug = ', debug)
      print('len(self.jedilat) = ', len(self.jedilat))
      print('len(self.gsilat) = ', len(self.gsilat))

    if(os.path.exists('jedi_only_obs.txt')):
      self.jedionlylat, self.jedionlylon = self.readdata('jedi_only_obs.txt')
      self.gsionlylat, self.gsionlylon = self.readdata('gsi_only_obs.txt')
      self.comlat, self.comlon = self.readdata('common_obs.txt')
    else:
      self.splitobs()

  def splitobs(self):
    self.jedimask = []
    self.gsimask = []
    self.comlat = []
    self.comlon = []
    self.gsiindex = np.linspace(0, len(self.gsilat), len(self.gsilat), endpoint=False, dtype=int)
    print('self.gsiindex = ', self.gsiindex)

    nlat = len(self.jedilat)
    i = 0
    xb = time.time()
    while(i < nlat):
      found, idx = self.findobs(lat=self.jedilat[i], lon=self.jedilon[i])
      if(found):
        self.comlat.append(self.jedilat[i])
        self.comlon.append(self.jedilon[i])
      
        self.jedimask.append(i)
        self.gsimask.append(idx)
       #print('deleted jedi %d, gsi %d, mask size: %d' %(i, idx, len(self.jedimask)))
      i += 1;
      xe = time.time()
      print('Loop on No. %d, len(self.gsimask)=%d, use time: %f' %(i, len(self.gsimask), xe - xb))
      xb = xe

    print('len(self.jedimask) = ', len(self.jedimask))
    print('len(self.gsimask) = ', len(self.gsimask))
    print('len(self.comlat) = ', len(self.comlat))

    self.jedionlylat = np.delete(self.jedilat, self.jedimask)
    self.jedionlylon = np.delete(self.jedilon, self.jedimask)

    self.gsionlylat = np.delete(self.gsilat, self.gsimask)
    self.gsionlylon = np.delete(self.gsilon, self.gsimask)

    self.writedata(self.jedionlylat, self.jedionlylon, 'jedi_only_obs.txt')
    self.writedata(self.gsionlylat, self.gsionlylon, 'gsi_only_obs.txt')
    self.writedata(self.comlat, self.comlon, 'common_obs.txt')

  def findobs(self, lat=100.0, lon=400.0):
    loopidx = np.delete(self.gsiindex, self.gsimask)
   #print('loopidx = ', loopidx)
    for n in loopidx:
      if(abs(lat - self.gsilat[n]) > self.delt):
        continue
      if(abs(lon - self.gsilon[n]) > self.delt):
        continue
      return 1, n
    return 0, 0

  def get_ori_jedi_latlon(self):
    return self.jedilat, self.jedilon

  def get_ori_gsi_latlon(self):
    return self.gsilat, self.gsilon

  def get_jedi_latlon(self):
    return self.jedionlylat, self.jedionlylon

  def get_gsi_latlon(self):
    return self.gsionlylat, self.gsionlylon

  def get_common_latlon(self):
    return self.comlat, self.comlon

  def writedata(self, lat, lon, flnm):
    OF = open(flnm, 'w')
    for n in range(len(lat)):
      OF.write("%f, %f\n" %(lat[n], lon[n]))
    OF.close()

  def readdata(self, flnm):
    lat = []
    lon = []
    INF = open(flnm, 'r')
    lines = INF.readlines()
    for line in lines:
      item = line.split(', ')
      lat.append(float(item[0]))
      lon.append(float(item[1]))
    INF.close()

    return lat, lon

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uvOnly = 0
  addobs = 1
  uselogp = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'uvOnly=', 'addobs=', 'uselogp='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uvOnly'):
      uvOnly = int(a)
    elif o in ('--addobs'):
      addobs = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('uvOnly = ', uvOnly)
  print('uvOnly = ', uvOnly)

#=======================================================================================================================
  if(addobs):
    jediobsfile = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/ioda_v2_data/obs/ncdiag.oper.ob.PT6H.sondes.2021-01-08T21:00:00Z.nc4'
   #gsiobsfile = '/work/noaa/gsienkf/weihuang/jedi/vis_tools/visfv3/jeff-all-obs/diag_conv_t_ges.2021010900_ensmean.nc4'
   #gsiobsfile = '/work/noaa/gsienkf/weihuang/jedi/vis_tools/visfv3/jeff-all-obs/diag_conv_uv_ges.2021010900_ensmean.nc4'
    gsiobsfile = 'jeff-runs/PSonly/diag_conv_ps_ges.2021010900_ensmean.nc4'

   #so = SplitObs(debug=debug, jedifile=jediobsfile, gsifile=gsiobsfile)
    so = SplitObs(debug=debug, jedifile=jediobsfile, gsifile=gsiobsfile, varname='/ObsValue/surface_pressure')

#=======================================================================================================================
  fv3_griddir = '/work/noaa/gsienkf/weihuang/tools/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

  if(uvOnly):
    gsi_bkg = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/jeff-gsi-run/uv-only/sfg_2021010900_fhr06_ensmean'
    gsi_anl = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/jeff-gsi-run/uv-only/sanl_2021010900_fhr06_ensmean.uv'

    datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.uvOnly/increment/'

    imgprefix = 'uvOnly_jedi-gsi_sondes'
    title_prefix = 'uvOnly JEDI-GSI Sondes Temperature at'
  else:
    gsi_bkg = 'jeff-runs/PSonly/sfg_2021010900_fhr06_ensmean'
    gsi_anl = 'jeff-runs/PSonly/sanl_2021010900_fhr06_ensmean'

   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.all/increment/'
   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.uvTq/increment/'
    datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.psOnly/increment/'

   #imgprefix = 'diff_uvTq_jedi-gsi_sondes'
   #title_prefix = 'uvTq JEDI-GSI Sondes Temperature at'
    imgprefix = 'diff_psOnly_jedi-gsi_sondes'
    title_prefix = 'PS only JEDI-GSI Sondes Temperature at'

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(fv3_griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s20210109.000000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)

#=======================================================================================================================
  nlon = 360
  nlat = nlon/2 + 1
  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

#=======================================================================================================================
  pg = PlotGaussian(debug=debug, output=output, bkg=gsi_bkg, anl=gsi_anl)
  lon1d, lat1d, var1d = pg.get_var('tmp')

  regrid_gsi = regridder(debug=debug, datafiles=[], gridspecfiles=[])
  gsi_var = regrid_gsi.interp2latlon_data(lon1d, lat1d, var1d, nlon=nlon, nlat=nlat, method='linear')

#=======================================================================================================================
  regrid_jedi = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)

  varname = 'T'
  jedi_var = regrid_jedi.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

 #print('var.ndim = ', var.ndim)
 #print('var.shape = ', var.shape)

#=======================================================================================================================
  var = jedi_var - gsi_var

#=======================================================================================================================
  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)
  clevs = np.arange(-0.5, 0.51, 0.01)
  cblevs = np.arange(-0.5, 0.6, 0.1)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

  gp.set_label('Temperature (K)')

#=======================================================================================================================
  levs = [1, 10, 23, 30, 52, 62]

  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(imgprefix, lev)
    title = '%s level %d' %(title_prefix, lev)
    gp.set_imagename(imgname)
    gp.set_title(title)

    if(addobs):
      gp.plot_without_marker(pvar)
      jedilat, jedilon = so.get_jedi_latlon()
      print('len(jedilat) = ', len(jedilat))
      gp.set_obs_latlon(obslat=jedilat, obslon=jedilon)
     #gp.add_obs_marker(marker='x', size=3, color='green')
      gp.add_obs_marker(marker='x', size=2, color='green')
      comlat, comlon = so.get_common_latlon()
      print('len(comlat) = ', len(comlat))
      gp.set_obs_latlon(obslat=comlat, obslon=comlon)
     #gp.add_obs_marker(marker='o', size=3, color='yellow')
      gp.add_obs_marker(marker='o', size=1, color='yellow')
      gsilat, gsilon = so.get_gsi_latlon()
      print('len(gsilat) = ', len(gsilat))
      gp.set_obs_latlon(obslat=gsilat, obslon=gsilon)
     #gp.add_obs_marker(marker='+', size=3, color='cyan')
      gp.add_obs_marker(marker='+', size=0.5, color='cyan')
      gp.display(output=output)
    else:
      gp.plot(pvar)

#=======================================================================================================================
  lons = [40, 105, 170, 270, 300]

  for lon in lons:
    pvar = var[:,:,lon]
    imgname = '%s_lon_%d.png' %(imgprefix, lon)
    title = '%s longitude %d' %(title_prefix, lon)
    gp.set_imagename(imgname)
    gp.set_title(title)
    if(uselogp):
      gp.plot_meridional_section_logp(pvar)
    else:
      gp.plot_meridional_section(pvar)

#=======================================================================================================================
  lats = [-30, 0, 45, 70]

  for lat in lats:
    pvar = var[:,90+lat,:]
    imgname = '%s_lat_%d.png' %(imgprefix, lat)
    title = '%s latitude %d' %(title_prefix, lat)
    gp.set_imagename(imgname)
    gp.set_title(title)
    if(uselogp):
      gp.plot_zonal_section_logp(pvar)
    else:
      gp.plot_zonal_section(pvar)

