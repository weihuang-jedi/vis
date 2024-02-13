import getopt
import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as ccrs
from cartopy import config
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker

import netCDF4

import numpy.ma as ma

from readIODA2Obs import ReadIODA2Obs

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.set_default()

#--------------------------------------------------------------------------------
  def plot(self, gridlon, gridlat, data=[], obsvar=[]):
    nrows = len(data)
    ncols = 1

   #set up the plot
    proj = ccrs.PlateCarree()
    proj2 = ccrs.Mollweide()

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw=dict(projection=proj),
                            figsize=(11,8.5))
 
   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

   #cslist = []

    for i in range(len(axs)):
      axs[i].set_global()

      pvar = data[i]

     #print('Plot No. ', i)
     #print('\tpvar.shape = ', pvar.shape)
     #print('\tgridlon.shape = ', gridlon.shape)
     #print('\tgridlat.shape = ', gridlat.shape)

      cyclic_data, cyclic_lon = add_cyclic_point(pvar, coord=gridlon)

      cs=axs[i].contourf(cyclic_lon, gridlat, cyclic_data, transform=proj,
                         cmap=self.cmapname)
     #cs=axs[i].contourf(cyclic_lon, gridlat, cyclic_data, transform=proj,
     #                   levels=self.clevs, extend=self.extend,
     #                   alpha=self.alpha, cmap=self.cmapname)
     #cs=axs[i].contourf(gridlon, gridlat, pvar, transform=proj,
     #                   levels=self.clevs, extend=self.extend,
     #                   alpha=self.alpha, cmap=self.cmapname)

     #cslist.append(cs)

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
     #ax.coastlines(resolution='110m')
      axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

      if(len(obsvar) < 1):
        msz = 3
        colors = 'cyan'
      else:
        print('obsvar min: %f, max: %f' %(np.min(obsvar), np.max(obsvar)))
        msz = np.zeros((len(obsvar),), dtype=int)
        colors = np.zeros((len(obsvar),), dtype=str)
        for n in range(len(obsvar)):
         #print('obsvar %d: %f' %(n, obsvar[n]))
          if(ma.is_masked(obsvar[n])):
            msz[n] = 1
            colors[n] = 'cycan'
          else:
            msz[n] = int(500.0*abs(obsvar[n]))+1
           #print('Size %d: %d' %(n, msz[n]))
            if(obsvar[n] < 0.0):
              colors[n] = 'blue'
            else:
              colors[n] = 'red'
         
     #adding marker:
      dotes = axs[i].scatter(self.obslon, self.obslat, s=msz, c=colors)

      axs[i].set_title(self.runname[i])

   #Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.8,
                        wspace=0.02, hspace=0.02)

   #Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.85])

   #Draw the colorbar
   #cbar=fig.colorbar(cslist[0], cax=cbar_ax, pad=self.pad, ticks=self.cblevs,
    cbar=fig.colorbar(cs, cax=cbar_ax, pad=self.pad, ticks=self.cblevs,
                      orientation='vertical', extend='both')

   #cbar.ax.get_yaxis().set_ticks([])
   #for lev in self.cblevs:
   #  txt = '%4.2f' %(lev)
   #  cbar.ax.text(.5, (lev + 0.5)/len(self.cblevs), txt, ha='center', va='center')

   #cbar.ax.get_yaxis().labelpad = 15
   #cbar.ax.set_ylabel('# of contacts', rotation=270)
    cbar.set_label(self.label, rotation=90)

   #Add a big title at the top
    plt.suptitle(self.title)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(self.imagename is None):
        imagename = 'sample.png'
      else:
        imagename = self.imagename
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

#--------------------------------------------------------------------------------
  def set_default(self):
    self.imagename = 'sample.png'

    self.runname = ['Linear Observer', 'NonLinear', 'Linear Observer - NonLinear']

   #cmapname = coolwarm, bwr, rainbow, jet, seismic
    self.cmapname = 'bwr'
   #self.cmapname = 'coolwarm'
   #self.cmapname = 'rainbow'
   #self.cmapname = 'jet'

    self.clevs = np.arange(-0.2, 0.21, 0.01)
    self.cblevs = np.arange(-0.2, 0.3, 0.1)

    self.extend = 'both'
    self.alpha = 0.5
    self.pad = 0.1
    self.orientation = 'horizontal'
    self.size = 'large'
    self.weight = 'bold'
    self.labelsize = 'medium'

    self.label = 'Unit (C)'
    self.title = 'Temperature Increment'

    self.obslon = []
    self.obslat = []
    self.add_obs_marker = False

#--------------------------------------------------------------------------------
  def set_label(self, label='Unit (C)'):
    self.label = label

#--------------------------------------------------------------------------------
  def set_title(self, title='Temperature Increment'):
    self.title = title

#--------------------------------------------------------------------------------
  def set_clevs(self, clevs=[]):
    self.clevs = clevs

#--------------------------------------------------------------------------------
  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

#--------------------------------------------------------------------------------
  def set_imagename(self, imagename):
    self.imagename = imagename

#--------------------------------------------------------------------------------
  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

#--------------------------------------------------------------------------------
  def set_obs_lonlat(self, obslon, obslat):
    self.obslon = obslon
    self.obslat = obslat

#--------------------------------------------------------------------------------
  def switch_marker_on(self):
    self.add_obs_marker = True

#--------------------------------------------------------------------------------
  def switch_marker_off(self):
    self.add_obs_marker = False

#--------------------------------------------------------------------------------
  def set_runname(self, runname = ['Linear Observer', 'NonLinear', 'Linear Observer - NonLinear']):
    self.runname = runname

#===================================================================================
class PlotObsOnMap():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output

    self.gp = GeneratePlot(debug=debug, output=output)

#-----------------------------------------------------------------------------------------
  def set_default(self):
    self.gridvarlist1 = ['T', 'ua', 'va', 'sphum', 'DELP', 'DZ', 'o3mr']
    self.gridvarlist2 = ['T', 'ua', 'va', 'sphum', 'DELP', 'DZ', 'o3mr']
    self.unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)',
                     'Unit (kg/kg)', 'Unit (Pa', 'Unit (m', 'Unit (ppm)']

   #self.gridvarlist1 = ['ua']
   #self.gridvarlist2 = ['ua']
   #self.unitlist = ['Unit (m/s)']

    self.clevs = np.arange(-1.0, 1.01, 0.01)
    self.cblevs = np.arange(-1.0, 1.1, 0.1)

    self.gp.set_clevs(clevs=self.clevs)
    self.gp.set_cblevs(cblevs=self.cblevs)

    self.gp.set_runname(runname=['Linear Observer', 'NonLinear', 'Linear Observer - NonLinear'])
    self.gp.switch_marker_on()

#--------------------------------------------------------------------------------
  def process(self, gridbase=None, gridcase=None,
                    obsbase=None, obscase=None,
                    grplist=None, varlist=None):
    self.set_default()

    if(os.path.exists(gridbase)):
      self.ncbase = netCDF4.Dataset(gridbase)
    else:
      print('File %s does not exist. Stop' %(gridbase))
      sys.exit(-1)

    if(os.path.exists(gridcase)):
      self.nccase = netCDF4.Dataset(gridcase)
    else:
      print('File %s does not exist. Stop' %(gridcase))
      sys.exit(-1)

    if(os.path.exists(gridbase)):
      self.nccase = netCDF4.Dataset(gridcase)
    else:
      print('File %s does not exist. Stop' %(gridbase))
      sys.exit(-1)

    if(os.path.exists(obsbase)):
      self.baseobs = ReadIODA2Obs(debug=debug, filename=obsbase)
    else:
      print('File %s does not exist. Stop' %(obsbase))
      sys.exit(-1)

    if(os.path.exists(obscase)):
      self.caseobs = ReadIODA2Obs(debug=debug, filename=obscase)
    else:
      print('File %s does not exist. Stop' %(obscase))
      sys.exit(-1)

    obslat, obslon = self.baseobs.get_latlon()
    obslon = np.where(obslon > 0, obslon, obslon+360.0)

    self.gp.set_obs_lonlat(obslon, obslat)

    self.variables = self.ncbase.variables

    self.gridlat = self.ncbase.variables['lat'][:]
    self.gridlon = self.ncbase.variables['lon'][:]

    for grp in grplist:
      self.grp = grp
      for nv in range(len(varlist)):
        self.varbase = varlist[nv]
        varname = '/%s/%s' %(grp, self.varbase)
        print('varname: ', varname)

        if('unknown' == varlist[nv]):
          continue

        base_var = self.baseobs.get_var(varname)
        case_var = self.caseobs.get_var(varname)
        diffvar = case_var - base_var
        self.plot(nv, diffvar)

    self.baseobs.close()
    self.caseobs.close()

    self.ncbase.close()
    self.nccase.close()

#-----------------------------------------------------------------------------------------
  def set_runname(self, runname=['Linear Observer', 'NonLinear', 'Linear Observer - NonLinear']):
    self.gp.set_runname(runname=runname)

#-----------------------------------------------------------------------------------------
  def switch_marker_on(self):
    self.gp.switch_marker_on()

#-----------------------------------------------------------------------------------------
  def set_obstype(self, obstype):
    self.obstype = obstype

#-----------------------------------------------------------------------------------------
  def set_enskind(self, enskind):
    self.enskind = enskind

#-----------------------------------------------------------------------------------------
  def plot(self, nv, diffvar):
    n = nv
    while(n >= len(self.gridvarlist1)):
      n -= len(self.gridvarlist1)

    grd1 = self.ncbase.variables[self.gridvarlist1[n]][0,:,:,:]
    grd2 = self.nccase.variables[self.gridvarlist2[n]][0,:,:,:]

    nlev, nlat, nlon = grd1.shape
   #print('grd1.shape = ', grd1.shape)
   #print('grd2.shape = ', grd2.shape)

    self.gp.set_label(self.unitlist[n])

    titlelabel = '%s %s %s' %(self.obstype, self.enskind, self.varbase)
    plotlabel = '%s_%s_%s' %(self.obstype, self.enskind, self.varbase)

    for lev in range(95, nlev, 10):
      v1 = grd1[lev,:,:]
      v2 = grd2[lev,:,:]
      dv = v2 - v1

      data = [v1, v2, dv]

      title = '%s %s at Level %d' %(titlelabel, self.gridvarlist1[n], lev)
      self.gp.set_title(title)

      print('Plotting ', title)
      print('\tv1.min: %f, v1.max: %f' %(np.min(v1), np.max(v1)))
      print('\tv2.min: %f, v2.max: %f' %(np.min(v2), np.max(v2)))

      imagename = '%s_%s_lev_%3.3d.png' %(plotlabel, self.gridvarlist1[n], lev)
      self.gp.set_imagename(imagename)

      self.gp.plot(self.gridlon, self.gridlat, data=data, obsvar=diffvar)

#===================================================================================
if __name__== '__main__':
  debug = 1
  output = 0
  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                                                'file1=', 'file2='])
  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--file1'):
      file1 = a
    elif o in ('--file2'):
      file2 = a
    else:
      assert False, 'unhandled option'

#--------------------------------------------------------------------------------
  runname = ['Linear Observer', 'NonLinear', 'Linear Observer - NonLinear']

  poom = PlotObsOnMap(debug=debug, output=output)

  poom.set_runname(runname=runname)
  poom.switch_marker_on()

#--------------------------------------------------------------------------------
  workdir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test'
  griddir = '%s/Data/analysis/letkf/gfs' %(workdir)
  obs_dir = '%s/Data' %(workdir)
  obsdir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data'

 #gridbase = '%s/mem000.nl/xainc.20201215_000000z.nc4' %(griddir)
 #gridcase = '%s/mem000.lo/xainc.20201215_000000z.nc4' %(griddir)

  gridbase = '%s/mem000.nl.getkf/xainc.20201215_000000z.nc4' %(griddir)
  gridcase = '%s/mem000.lo.getkf/xainc.20201215_000000z.nc4' %(griddir)

#--------------------------------------------------------------------------------
  nlexps = ['nl', 'nl.getkf']
  loexps = ['lo', 'lo.getkf']
  enslist = ['letkf-gfs', 'lgetkf-geos']

  obslist = ['aircraft', 'scatwind', 'sfc']

  varbaselist = [['airTemperature', 'windEastward', 'windNorthward'],
                 ['unknown', 'windEastward', 'windNorthward'],
                 ['stationPressure']]

  grplist = ['ombg']

#--------------------------------------------------------------------------------
  for ne in range(len(nlexps)):
    gridbase = '%s/mem000.%s/xainc.20201215_000000z.nc4' %(griddir, nlexps[ne])
    gridcase = '%s/mem000.%s/xainc.20201215_000000z.nc4' %(griddir, loexps[ne])

    print('gridbase: ', gridbase)
    print('gridcase: ', gridcase)

    enskind = enslist[ne]
    poom.set_enskind(enskind)

    for no in range(len(obslist)):
      obstype = obslist[no]
      poom.set_obstype(obstype)

      obsname = '%s_%s_2020121500_m.nc4' %(obstype, enslist[ne])
      basefile = '%s/hofx.%s/%s' %(obs_dir, nlexps[ne], obsname)
      casefile = '%s/hofx.%s/%s' %(obs_dir, loexps[ne], obsname)

      print('basefile: ', basefile)
      print('casefile: ', casefile)

      poom.process(gridbase=gridbase, gridcase=gridcase,
                   obsbase=basefile, obscase=casefile,
                   grplist=grplist, varlist=varbaselist[no])

