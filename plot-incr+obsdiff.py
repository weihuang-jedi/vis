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

from modelVerticalpressure import ModelVerticalPressure

#=========================================================================
class GeneratePlot():
  def __init__(self, debug=0, output=0, akbkfile='akbk127.nc4'):
    self.debug = debug
    self.output = output

    self.mvp = ModelVerticalPressure(debug=debug, filename=akbkfile)
    prs = self.mvp.get_pressure()
   #print('len(prs) = ', len(prs))
   #for n in range(len(prs)):
   #  print('Level %d pressure %f' %(n, prs[n]))

    self.logp = self.mvp.get_logp()
    self.markpres = self.mvp.get_markpres()
    self.marklogp = self.mvp.get_marklogp()

    self.precision = 1

    self.set_default()

#--------------------------------------------------------------------------------
  def plot(self, lons, lats, data=[], obsvar=[]):
    if(self.debug):
      print('obsvar min: %f, max: %f' %(np.min(obsvar), np.max(obsvar)))

    nrows = len(data)
    ncols = 1

   #set up the plot
    proj = ccrs.PlateCarree()

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols,
                            subplot_kw=dict(projection=proj),
                            figsize=(11,8.5))
 
   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      pvar = data[i]

      cyclic_data, cyclic_lons = add_cyclic_point(pvar, coord=lons)

      if(i > 1):
        cyclic_data *= 10.0
        title = '%s magnified 10 time' %(self.runname[i])
      else:
        title = self.runname[i]

      cs=axs[i].contourf(cyclic_lons, lats, cyclic_data, transform=proj,
                         levels=self.clevs, extend=self.extend,
                         alpha=self.alpha, cmap=self.cmapname)

      axs[i].set_extent([-180, 180, -90, 90], crs=proj)
      axs[i].coastlines(resolution='auto', color='k')
      axs[i].gridlines(color='lightgrey', linestyle='-', draw_labels=True)

      axs[i].set_title(title)

      if(len(obsvar) < 1):
        msz = 3
        colors = 'cyan'
      else:
       #print('obsvar min: %f, max: %f' %(np.min(obsvar), np.max(obsvar)))
        msz = np.zeros((len(obsvar),), dtype=int)
        colors = np.zeros((len(obsvar),), dtype=str)
        for n in range(len(obsvar)):
          if(ma.is_masked(obsvar[n])):
            msz[n] = 1
            colors[n] = 'cycan'
          else:
            msz[n] = int(500.0*abs(obsvar[n]))+1
            if(obsvar[n] < 0.0):
              colors[n] = 'cyan'
            else:
              colors[n] = 'magenta'

     #adding marker:
     #dotes = axs[i].scatter(self.obslon, self.obslat, s=msz, c=colors)

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

#--------------------------------------------------------------------------------
  def plot_meridional_section(self, lat, lev, data, title, imagename, suptitle):
    nlev = len(lev)
    nrows = 1
    ncols = len(data)

    self.imagename = imagename

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols, figsize=(11,8.5))

   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      pvar = data[i]

     #print('pvar.shape = ', pvar.shape)
     #print('lat = ', lat)
     #print('lev = ', lev)

      self.title = title[i]

      contfill = axs[i].contourf(lat, -lev[::-1], pvar[::-1,:],
                                 levels=self.clevs, extend=self.extend,
                                 alpha=self.alpha, cmap=self.cmapname)

      axs[i].set_title(self.title)

      major_ticks_top=np.linspace(-90,90,7)
      axs[i].set_xticks(major_ticks_top)

      minor_ticks_top=np.linspace(-90,90,19)
      axs[i].set_xticks(minor_ticks_top,minor=True)

      intv = int(nlev/10)
      nbeg = -10*intv
      major_ticks_top=np.linspace(nbeg,0,intv+1)
      axs[i].set_yticks(major_ticks_top)

      axs[i].grid(b=True, which='major', color='green', linestyle='-', alpha=0.5)
      axs[i].grid(b=True, which='minor', color='green', linestyle='dotted', alpha=0.2)

      cb = plt.colorbar(contfill, ax=axs[i], orientation=self.orientation,
                        pad=self.pad, ticks=self.cblevs)

      cb.set_label(label=self.label, size=self.size, weight=self.weight)

      cb.ax.tick_params(labelsize=self.labelsize)
      if(self.precision == 0):
        cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 1):
        cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 2):
        cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in self.cblevs], minor=False)
      else:
        cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in self.cblevs], minor=False)

   #Add a big title at the top
    plt.suptitle(suptitle)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(imagename is None):
        imagename = 't_aspect.png'
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

#--------------------------------------------------------------------------------
  def plot_meridional_section_logp(self, lat, data, title, imagename, suptitle):
    nrows = 1
    ncols = len(data)

    self.imagename = imagename

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols, figsize=(11,8.5))

   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      pvar = data[i]

     #print('pvar.shape = ', pvar.shape)
     #print('lat = ', lat)
     #print('lev = ', lev)

      self.title = title[i]

      contfill = axs[i].contourf(lat, self.logp[::-1], pvar[::-1,:],
                                 levels=self.clevs, extend=self.extend,
                                 alpha=self.alpha, cmap=self.cmapname)
      axs[i].set_title(self.title)

      major_ticks_top=np.linspace(-90,90,7)
      axs[i].set_xticks(major_ticks_top)

      minor_ticks_top=np.linspace(-90,90,19)
      axs[i].set_xticks(minor_ticks_top,minor=True)

      axs[i].set_yticks(self.marklogp)
      axs[i].set_ylabel('Unit: hPa')

      yticklabels = []
      for p in self.markpres:
        lbl = '%d' %(int(p+0.1))
        yticklabels.append(lbl)
      axs[i].set_yticklabels(yticklabels)

     #axs[i].grid(b=True, which='major', color='green', linestyle='-', alpha=0.5)
     #axs[i].grid(b=True, axis='x', which='minor', color='green', linestyle='dotted', alpha=0.2)

      axs[i].grid(b=True, which='major', color='green', linestyle='-', alpha=0.5)
      axs[i].grid(b=True, which='minor', color='green', linestyle='dotted', alpha=0.2)

      cb = plt.colorbar(contfill, ax=axs[i], orientation=self.orientation,
                        pad=self.pad, ticks=self.cblevs)

      cb.set_label(label=self.label, size=self.size, weight=self.weight)

      cb.ax.tick_params(labelsize=self.labelsize)
      if(self.precision == 0):
        cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 1):
        cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 2):
        cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in self.cblevs], minor=False)
      else:
        cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in self.cblevs], minor=False)

   #Add a big title at the top
    plt.suptitle(suptitle)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(imagename is None):
        imagename = 't_aspect.png'
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

#--------------------------------------------------------------------------------
  def plot_zonal_section(self, lon, lev, data, title, imagename, suptitle):
    nlev = len(lev)
    nrows = 1
    ncols = len(data)

    self.imagename = imagename

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols, figsize=(11,8.5))

   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      pvar = data[i]

      self.title = title[i]

      contfill = axs[i].contourf(lon, -lev[::-1], pvar[::-1,:],
                                 levels=self.clevs, extend=self.extend,
                                 alpha=self.alpha, cmap=self.cmapname)

      axs[i].set_title(self.title)

      minor_ticks_top=np.linspace(0,360,13)
      axs[i].set_xticks(minor_ticks_top,minor=True)

      minor_ticks_top=np.linspace(0,360,37)
      axs[i].set_xticks(minor_ticks_top,minor=True)

      intv = int(nlev/10)
      nbeg = -10*intv
      major_ticks_top=np.linspace(nbeg,0,intv+1)
      axs[i].set_yticks(major_ticks_top)

      axs[i].grid(b=True, which='major', color='green', linestyle='-', alpha=0.5)
      axs[i].grid(b=True, which='minor', color='green', linestyle='dotted', alpha=0.2)

      cb = plt.colorbar(contfill, ax=axs[i], orientation=self.orientation,
                        pad=self.pad, ticks=self.cblevs)

      cb.set_label(label=self.label, size=self.size, weight=self.weight)

      cb.ax.tick_params(labelsize=self.labelsize)
      if(self.precision == 0):
        cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 1):
        cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 2):
        cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in self.cblevs], minor=False)
      else:
        cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in self.cblevs], minor=False)

   #Add a big title at the top
    plt.suptitle(suptitle)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(imagename is None):
        imagename = 't_aspect.png'
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

#--------------------------------------------------------------------------------
  def plot_zonal_section_logp(self, lon, data, title, imagename, suptitle):
    nrows = 1
    ncols = len(data)

    self.imagename = imagename

    fig, axs = plt.subplots(nrows=nrows,ncols=ncols, figsize=(11,8.5))

   #axs is a 2 dimensional array of `GeoAxes`. Flatten it into a 1-D array
    axs=axs.flatten()

    for i in range(len(axs)):
      pvar = data[i]

      self.title = title[i]

      contfill = axs[i].contourf(lon, self.logp[::-1], pvar[::-1,:],
                                 levels=self.clevs, extend=self.extend,
                                 alpha=self.alpha, cmap=self.cmapname)

      axs[i].set_title(self.title)

      minor_ticks_top=np.linspace(0,360,13)
      axs[i].set_xticks(minor_ticks_top,minor=True)

      minor_ticks_top=np.linspace(0,360,37)
      axs[i].set_xticks(minor_ticks_top,minor=True)

      axs[i].set_yticks(self.marklogp)
      axs[i].set_ylabel('Unit: hPa')

      yticklabels = []
      for p in self.markpres:
        lbl = '%d' %(int(p+0.1))
        yticklabels.append(lbl)
      axs[i].set_yticklabels(yticklabels)

      axs[i].grid(b=True, which='major', color='green', linestyle='-', alpha=0.5)
      axs[i].grid(b=True, which='minor', color='green', linestyle='dotted', alpha=0.2)

      cb = plt.colorbar(contfill, ax=axs[i], orientation=self.orientation,
                        pad=self.pad, ticks=self.cblevs)

      cb.set_label(label=self.label, size=self.size, weight=self.weight)

      cb.ax.tick_params(labelsize=self.labelsize)
      if(self.precision == 0):
        cb.ax.set_xticklabels(['{:.0f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 1):
        cb.ax.set_xticklabels(['{:.1f}'.format(x) for x in self.cblevs], minor=False)
      elif(self.precision == 2):
        cb.ax.set_xticklabels(['{:.2f}'.format(x) for x in self.cblevs], minor=False)
      else:
        cb.ax.set_xticklabels(['{:.3f}'.format(x) for x in self.cblevs], minor=False)

   #Add a big title at the top
    plt.suptitle(suptitle)

    fig.canvas.draw()
    plt.tight_layout()

    if(self.output):
      if(imagename is None):
        imagename = 't_aspect.png'
      plt.savefig(imagename)
      plt.close()
    else:
      plt.show()

#===================================================================================
class PlotObsOnMap():
  def __init__(self, debug=0, output=0):
    self.debug = debug
    self.output = output
    self.casename = 'unknown'

    self.gp = GeneratePlot(debug=debug, output=output)

#--------------------------------------------------------------------------------
  def set_casename(self, casename='unknown'):
    self.casename = casename

#-----------------------------------------------------------------------------------------
  def set_default(self):
    self.gridvarlist1 = ['T', 'ua', 'va', 'sphum', 'DELP', 'DZ', 'o3mr']
    self.gridvarlist2 = ['T', 'ua', 'va', 'sphum', 'DELP', 'DZ', 'o3mr']
    self.unitlist = ['Unit (C)', 'Unit (m/s)', 'Unit (m/s)',
                     'Unit (kg/kg)', 'Unit (Pa', 'Unit (m', 'Unit (ppm)']

   #self.gridvarlist1 = ['ua']
   #self.gridvarlist2 = ['ua']
   #self.unitlist = ['Unit (m/s)']

   #self.clevs = np.arange(-1.0, 1.01, 0.01)
   #self.cblevs = np.arange(-1.0, 1.1, 0.1)

    self.clevs = np.arange(-1.0, 1.02, 0.02)
    self.cblevs = np.arange(-1.0, 1.2, 0.2)

    self.gp.set_clevs(clevs=self.clevs)
    self.gp.set_cblevs(cblevs=self.cblevs)

    self.gp.set_runname(runname=['Linear Observer', 'NonLinear', 'Linear Observer - NonLinear'])
    self.gp.switch_marker_on()

#--------------------------------------------------------------------------------
  def process(self, gridbase=None, gridcase=None,
                    obsbase=None, obscase=None,
                    grplist=None, varlist=None, casename='unknown'):
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
    self.set_casename(casename=casename)

    self.variables = self.ncbase.variables

    self.gridlat = self.ncbase.variables['lat'][:]
    self.gridlon = self.ncbase.variables['lon'][:]

   #for n in range(len(self.gridlon)):
   #  print('No %d lon %f' %(n, self.gridlon[n]))
   #for n in range(len(self.gridlat)):
   #  print('No %d lat %f' %(n, self.gridlat[n]))

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

    lev = np.arange(0.0, float(nlev), 1.0)

    self.gp.set_label(self.unitlist[n])

   #for lon in [8, 15, 32]:
    for lon in []:
      v1 = grd1[:,:,lon]
      v2 = grd2[:,:,lon]
      dv = v2 - v1
      print('\tv1.min: %f, v1.max: %f' %(np.min(v1), np.max(v1)))
      print('\tv2.min: %f, v2.max: %f' %(np.min(v2), np.max(v2)))
      print('\tdv.min: %f, dv.max: %f' %(np.min(dv), np.max(dv)))

      data = [v1, v2, dv]

      suptitle = '%s Meridional Section at Lon %f' %(self.casename, self.gridlon[lon])
      title0 = 'NonLinear'
      title1 = 'Linear'
      title2 = 'Linear - NonLinear'

      title = [title0, title1, title2]

      imagename = '%s_meridional_section_at_lon_%f.png' %(self.casename, self.gridlon[lon])
      self.gp.plot_meridional_section(self.gridlat, lev, data, title, imagename, suptitle)
      imagename = '%s_meridional_logp_section_at_lon_%f.png' %(self.casename, self.gridlon[lon])
      self.gp.plot_meridional_section_logp(self.gridlat, data, title, imagename, suptitle)

   #for lat in [17, 21]:
    for lat in []:
      v1 = grd1[:,lat,:]
      v2 = grd2[:,lat,:]
      dv = v2 - v1
      print('\tv1.min: %f, v1.max: %f' %(np.min(v1), np.max(v1)))
      print('\tv2.min: %f, v2.max: %f' %(np.min(v2), np.max(v2)))
      print('\tdv.min: %f, dv.max: %f' %(np.min(dv), np.max(dv)))

      data = [v1, v2, dv]

      suptitle = '%s Zonal Section at Lat %f' %(self.casename, self.gridlon[lon])
      title0 = 'NonLinear'
      title1 = 'Linear'
      title2 = 'Linear - NonLinear'

      title = [title0, title1, title2]

      imagename = '%s_zonal_section_at_lat_%f.png' %(self.casename, self.gridlat[lat])
      self.gp.plot_zonal_section(self.gridlon, lev, data, title, imagename, suptitle)
      imagename = '%s_zonal_logp_section_at_lat_%f.png' %(self.casename, self.gridlat[lat])
      self.gp.plot_zonal_section_logp(self.gridlon, data, title, imagename, suptitle)

    titlelabel = '%s %s %s %s' %(self.casename, self.obstype, self.enskind, self.varbase)
    plotlabel = '%s_%s_%s_%s' %(self.casename, self.obstype, self.enskind, self.varbase)

    for lev in range(95, nlev, 10):
      v1 = grd1[lev,:,:]
      v2 = grd2[lev,:,:]
      dv = v2 - v1
      print('\tv1.min: %f, v1.max: %f' %(np.min(v1), np.max(v1)))
      print('\tv2.min: %f, v2.max: %f' %(np.min(v2), np.max(v2)))
      print('\tdv.min: %f, dv.max: %f' %(np.min(dv), np.max(dv)))

      data = [v1, v2, dv]

      if(diffvar.ndim > 1):
        for channel in range(len(diffvar[0])):
          title = 'channel %d %s %s at Level %d' %(channel, titlelabel, self.gridvarlist1[n], lev)
          self.gp.set_title(title)

          print('Plotting ', title)

          imagename = 'channel_%d_%s_%s_lev_%3.3d.png' %(channel, plotlabel, self.gridvarlist1[n], lev)
          self.gp.set_imagename(imagename)

          self.gp.plot(self.gridlon, self.gridlat, data=data, obsvar=diffvar[:, channel])
      else:
        title = '%s %s at Level %d' %(titlelabel, self.gridvarlist1[n], lev)
        self.gp.set_title(title)

        print('Plotting ', title)

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
 #nlexps = ['nl', 'nl.getkf']
 #loexps = ['lo', 'lo.getkf']
 #enslist = ['letkf-gfs', 'lgetkf-geos']

  nlexps = ['nl.getkf', 'nl']
  loexps = ['lo.getkf', 'lo']
  enslist = ['lgetkf-geos', 'letkf-gfs']

 #obslist = ['aircraft', 'scatwind', 'sfc', 'amsua_n19']

 #varbaselist = [['airTemperature', 'windEastward', 'windNorthward'],
 #               ['unknown', 'windEastward', 'windNorthward'],
 #               ['stationPressure'],
 #               ['brightnessTemperature']]
  obslist = ['sfcship']
  varbaselist = [['stationPressure']]

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

      obsname = '%s_%s_2020121500_s.nc4' %(obstype, enslist[ne])
      basefile = '%s/hofx.%s/%s' %(obs_dir, nlexps[ne], obsname)
      casefile = '%s/hofx.%s/%s' %(obs_dir, loexps[ne], obsname)

      print('basefile: ', basefile)
      print('casefile: ', casefile)

      poom.process(gridbase=gridbase, gridcase=gridcase,
                   obsbase=basefile, obscase=casefile,
                   grplist=grplist, varlist=varbaselist[no],
                   casename=enslist[ne])

