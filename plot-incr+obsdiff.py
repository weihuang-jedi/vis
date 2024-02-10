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

  def plot(self, gridlon, gridlat, data=[]):
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

    for i in range(len(axs)):
      axs[i].set_global()

      pvar = data[i]

     #print('Plot No. ', i)
     #print('\tpvar.shape = ', pvar.shape)

     #print('lat = ', gridlat)
     #print('lon = ', gridlon)

      cyclic_data, cyclic_lon = add_cyclic_point(pvar, coord=gridlon)

      cs=axs[i].contourf(cyclic_lon, gridlat, cyclic_data, transform=proj,
                         cmap=self.cmapname)
     #cs=axs[i].contourf(cyclic_lon, gridlat, cyclic_data, transform=proj,
     #                   levels=self.clevs, extend=self.extend,
     #                   alpha=self.alpha, cmap=self.cmapname)
     #cs=axs[i].contourf(gridlon, gridlat, pvar, transform=proj,
     #                   levels=self.clevs, extend=self.extend,
     #                   alpha=self.alpha, cmap=self.cmapname)

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
     #ax.coastlines(resolution='110m')
     #axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

     #adding dotes:
     #dotes = self.basemap.plot(x, y, 'bo', markersize=12)
     #dotes = self.basemap.plot(x, y, 'bo', markersize=6)

      marker = '.'
      if(self.svar is None):
        msz = 3
        colors = 'cyan'
      else:
        print('svar min: %f, max: %f' %(np.min(self.svar), np.max(self.svar)))
        msz = np.zeros((len(self.svar),))
        colors = np.zeros((len(self.svar),), dtype=str)
        for n in range(len(self.svar)):
         #print('self.svar %d: %f' %(n, self.svar[n]))
          if(ma.is_masked(self.svar[n])):
            msz[n] = 1
            colors[n] = 'cycan'
          else:
            msz[n] = int(500.0*abs(self.svar[n]))+1
           #print('Size %d: %d' %(n, msz[n]))
            if(self.svar[n] < 0.0):
              colors[n] = 'blue'
            else:
              colors[n] = 'red'
         
      dotes = axs[i].scatter(self.obslon, self.obslat, s=msz, c=colors)
     #dotes = axs[i].scatter(self.obslon, self.obslat,
     #                       marker=marker, s=msz, color=color)

      axs[i].set_title(self.runname[i])

   #Adjust the location of the subplots on the page to make room for the colorbar
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.05, right=0.8,
                        wspace=0.02, hspace=0.02)

   #Add a colorbar axis at the bottom of the graph
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.85])

   #Draw the colorbar
    cbar=fig.colorbar(cs, cax=cbar_ax, pad=self.pad, ticks=self.cblevs,
                      orientation='vertical')

    cbar.set_label(self.label, rotation=90)

   #Add a big title at the top
    plt.suptitle(self.title)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(self.imagename is None):
        imagename = 't_aspect.png'
      else:
        imagename = self.imagename
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

  def set_default(self):
    self.imagename = 'sample.png'

    self.runname = ['Halo', 'RR', 'RR - Halo']

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

  def set_label(self, label='Unit (C)'):
    self.label = label

  def set_title(self, title='Temperature Increment'):
    self.title = title

  def set_clevs(self, clevs=[]):
    self.clevs = clevs

  def set_cblevs(self, cblevs=[]):
    self.cblevs = cblevs

  def set_imagename(self, imagename):
    self.imagename = imagename

  def set_cmapname(self, cmapname):
    self.cmapname = cmapname

  def set_obs_lonlat(self, obslon, obslat, svar=None):
    self.obslon = obslon
    self.obslat = obslat
    self.svar = svar

  def switch_marker_on(self):
    self.add_obs_marker = True

  def switch_marker_off(self):
    self.add_obs_marker = False

#--------------------------------------------------------------------------------
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
  griddir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data/analysis/letkf/gfs'

 #gridbase = '%s/mem000.nl/xainc.20201215_000000z.nc4' %(griddir)
 #gridcase = '%s/mem000.lo/xainc.20201215_000000z.nc4' %(griddir)

  gridbase = '%s/mem000.nl.getkf/xainc.20201215_000000z.nc4' %(griddir)
  gridcase = '%s/mem000.lo.getkf/xainc.20201215_000000z.nc4' %(griddir)

#--------------------------------------------------------------------------------
  obsdir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data'
 #filename = 'scatwind_letkf-gfs_2020121500_m.nc4'
 #filename = 'aircraft_letkf-gfs_2020121500_m.nc4'
 #filename = 'aircraft_lgetkf-geos_2020121500_m.nc4'
  filename = 'scatwind_lgetkf-geos_2020121500_m.nc4'

#--------------------------------------------------------------------------------
 #basefile = '%s/hofx.nl/%s' %(obsdir, filename)
  basefile = '%s/hofx.nl.getkf/%s' %(obsdir, filename)
  baseio = ReadIODA2Obs(debug=debug, filename=basefile)
  obslat, obslon = baseio.get_latlon()
  obslon = np.where(obslon > 0, obslon, obslon+360.0)

 #casefile = '%s/hofx.lo/%s' %(obsdir, filename)
  casefile = '%s/hofx.lo.getkf/%s' %(obsdir, filename)
  caseio = ReadIODA2Obs(debug=debug, filename=casefile)

#--------------------------------------------------------------------------------
  varbase = 'windEastward'
  varname = '/ombg/%s' %(varbase)
 #varname = '/hofx_y_mean_xb0/%s' %(varbase)
 #varname = '/hofx0/%s' %(varbase)

 #for group in ['hofx_y_mean_xb0', 'hofx0', 'hofx1', 'hofx02', 'hofx0_1', 'hofx0_2', 'hofx0_3', 'hofx1_1', 'hofx1_2', 'hofx1_3', 'oman', 'ombg']:
 #for group in ['hofx_y_mean_xb0', 'hofx0', 'hofx0_1', 'hofx0_2', 'hofx0_3', 'hofx1_1', 'hofx1_2', 'hofx1_3', 'oman', 'ombg']:
  for group in ['hofx_y_mean_xb0', 'oman', 'ombg']:
    varname = '/%s/%s' %(group, varbase)
    base_var = baseio.get_var(varname)
    case_var = caseio.get_var(varname)

    print('varname: ', varname)
    print('len(case_var) = ', len(case_var))

    diffvar = case_var - base_var
   #print('\tbase_var:\n', base_var)
   #print('\tcase_var:\n', case_var)

   #print('\ttotal base_var.max: %f, base_var.min: %f' %(np.max(base_var), np.min(base_var)))
   #print('\ttotal case_var.max: %f, case_var.min: %f' %(np.max(case_var), np.min(case_var)))
    print('\tdiffvar.max: %f, diffvar.min: %f' %(np.max(diffvar), np.min(diffvar)))

    delt = 1.0e-10

    nd = 0
    n = 0

    for i in range(len(diffvar)):
        n += 1
        if(abs(diffvar[i]) > delt):
       #print('%d: c %f, b %f, diff: %f' %(n+1, case_var[i,channel], base_var[i,channel], diff[i,channel]))
          nd += 1

    if(nd):
      print('Len(var): %d, nd: %d, mindif: %f, maxdif: %f' %(case_var.size, nd, np.min(diffvar), np.max(diffvar)))

 #sys.exit(0)

#--------------------------------------------------------------------------------
  gp = GeneratePlot(debug=debug, output=output)

  ncbase = netCDF4.Dataset(gridbase)
  nccase = netCDF4.Dataset(gridcase)
  gridlat = ncbase.variables['lat'][:]
  gridlon = ncbase.variables['lon'][:]

#-----------------------------------------------------------------------------------------
  gp.set_obs_lonlat(obslon, obslat, svar=diffvar)
  gp.switch_marker_on()

#-----------------------------------------------------------------------------------------
  clevs = np.arange(-1.0, 1.01, 0.01)
  cblevs = np.arange(-1.0, 1.1, 0.1)

  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#-----------------------------------------------------------------------------------------
 #varlist1 = ['T', 'ua', 'va', 'sphum', 'delp', 'DZ', 'o3mr']
 #varlist2 = ['T', 'ua', 'va', 'sphum', 'delp', 'DZ', 'o3mr']
 #unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)',
 #            'Unit (kg/kg)', 'Unit (Pa', 'Unit (m', 'Unit (ppm)']

  varlist1 = ['ua']
  varlist2 = ['ua']
  unitlist = ['Unit (m/s)']

#-----------------------------------------------------------------------------------------
  for n in range(len(varlist1)):
    var1 = ncbase.variables[varlist1[n]][0,:,:,:]
    var2 = nccase.variables[varlist2[n]][0,:,:,:]

    nlev, nlat, nlon = var1.shape
    print('var1.shape = ', var1.shape)
    print('var2.shape = ', var2.shape)

    gp.set_label(unitlist[n])

    for lev in range(95, nlev, 10):
      v1 = var1[lev,:,:]
      v2 = var2[lev,:,:]
      dv = v2 - v1

      data = [v1, v2, dv]

      title = '%s at Level %d' %(varlist1[n], lev)
      gp.set_title(title)

      print('Plotting ', title)
      print('\tv1.shape = ', v1.shape)
      print('\tv2.shape = ', v2.shape)
 
      print('\tv1.max: %f, v1.min: %f' %(np.max(v1), np.min(v1)))
      print('\tv2.max: %f, v2.min: %f' %(np.max(v2), np.min(v2)))

      imagename = '%s_lev_%3.3d.png' %(varlist1[n], lev)
      gp.set_imagename(imagename)

      gp.plot(gridlon, gridlat, data=data)

#-----------------------------------------------------------------------------------------
  ncbase.close()
  nccase.close()

