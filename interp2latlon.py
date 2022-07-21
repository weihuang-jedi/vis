import xarray as xr
import numpy as np
import xesmf as xe
import os, sys
import glob
import getopt

#-------------------------------------------------------------------------------------------
class Ocean2LatLon():
  def __init__(self, debug=0, wgtdir='grids/'):
    self.debug = debug
    self.wgtdir = wgtdir

    if(self.debug):
      print('debug = ', debug)
      print('wgtdir = ', wgtdir)

    self.latlonfilename = 'grids/ocn_2014_01.nc'

   #Read in tri-polar rot vars.
    self.grid_in = xr.open_dataset(self.latlonfilename)
    self.cos_rot = self.grid_in.cos_rot
    self.sin_rot = self.grid_in.sin_rot
    self.grid_in.close()

   #specify an default input/output resolution
    self.ires = 'mx100'
    self.ores='360x180'

    self.set_ires(ires=self.ires)

    self.has_wgtfile = 0

    self.gen_wgtfile(infilename=self.latlonfilename)

  def set_ires(self, ires='mx100'):
    self.ires = ires

  def set_ores(self, ores='360x180'):
    self.ores = ores

  def set_wgtdir(self, wgtdir='./'):
    self.wgtdir = wgtdir

  def gen_wgtfile(self, infilename=None):
    if(not os.path.isfile(infilename)):
      print('input file %s does not exist. Stop' %infilename)
      sys.exit(-1)

   #interpolation of tripolar t-points to 1-degree grid
    self.t2t_wgtfile = '%s%s.Ct.%s.Ct.bilinear.nc' %(self.wgtdir, self.ires, self.ores)

   #OPEN 1 degree history file to get input grid
    ds_in=xr.open_dataset(infilename)

   #rename lat and lons from grid files for ESMF interpolation
    output_t_grid=ds_in.rename({'geolon': 'lon', 'geolat': 'lat'})

    print('output_t_grid=', output_t_grid)

   #define target grid
    lon1d=np.arange(0.5,360.5,1.0)
    lat1d=np.arange(-89.5,90.5,1.0)

    lons, lats=np.meshgrid(lon1d,lat1d)
    da_out_lons=xr.DataArray(lons,dims=['nx','ny'])
    da_out_lats=xr.DataArray(lats,dims=['nx','ny'])
    ds_out_lons=da_out_lons.to_dataset(name='lon')
    ds_out_lats=da_out_lats.to_dataset(name='lat')
    ds_out=xr.merge([ds_out_lons, ds_out_lats])

    self.rg_tt = xe.Regridder(output_t_grid, ds_out, 'bilinear',
                              periodic=True, filename=self.t2t_wgtfile)

    if(not os.path.isfile(self.t2t_wgtfile)):
      print('weight file %s does not exist. Stop' %self.t2t_wgtfile)
      sys.exit(-1)

    self.has_wgtfile = 1
    print('gen_wgtfile done!')

  def process_file(self, infilename=None, outfilename=None):
    if(not os.path.isfile(infilename)):
      print('input file %s does not exist. Stop' %infilename)
      sys.exit(-1)
      
    if(outfilename is None or infilename == outfilename):
      path, flnm = os.path.split(infilename)
      if(flnm.find('.nc') > 0):
        namestr = '_%s.nc' %(self.ores)
        outfilename =flnm.replace('.nc', namestr)
      else:
        outfilename ='./%s_%s.nc' %(flnm, self.ores)

    if(not self.has_wgtfile):
      self.gen_wgtfile(infilename=self.latlonfilename)
      self.has_wgtfile = 1

   #output grid
    lon1d=np.arange(0.5,360.5,1.0)
    lat1d=np.arange(-89.5,90.5,1.0)
    lons,lats=np.meshgrid(lon1d,lat1d)

    da_out_lons=xr.DataArray(lons,dims=['nx','ny'])
    da_out_lats=xr.DataArray(lats,dims=['nx','ny'])
    ds_out_lons=da_out_lons.to_dataset(name='lon')
    ds_out_lats=da_out_lats.to_dataset(name='lat')

    grid_out=xr.merge([ds_out_lons,ds_out_lats])

    print('Processing file: ', infilename)

    ds_in=xr.open_dataset(infilename)
    ds_out=[]
   #print('ds_in.keys(): ', ds_in.keys())
    for i in list(ds_in.keys()):
      print('\tWorking on: ', i)
     #print('\t\tds_in[i].coords = ', ds_in[i].coords)
      if len(ds_in[i].coords) > 2:
        coords=ds_in[i].coords.to_index()
       #print('\tcoords.names = ', coords.names)
        pos='skip'

        if coords.names[0] != 'Time':
          print('\t\tcoords.names[:] = ', coords.names[:])
        else:
          if i=='Salt':
            pos='T'
          elif i=='Temp':
            pos='T'
          elif i=='ave_ssh':
            pos='T'
          elif i=='h':
            pos='T'
          else:
            pos='skip'

        if coords.names[1] == 'zaxis_1':  # 3-dimensional data
          if pos=='T':  # interplate to lat-lon grid
            interp_out= self.rg_tt(ds_in[i].values)
            da_out=xr.DataArray(interp_out,dims=['time','lay','lat','lon'])                    
            da_out.attrs['long_name']=ds_in[i].long_name
            da_out.attrs['units']=ds_in[i].units
            ds_out.append(da_out.to_dataset(name=i))
        
        else: # 2 dimension data
          if pos=='T':  # interplate to lat-lon grid
            interp_out= self.rg_tt(ds_in[i].values)
            da_out=xr.DataArray(interp_out,dims=['time','lat','lon'])                    
            da_out.attrs['long_name']=ds_in[i].long_name
            da_out.attrs['units']=ds_in[i].units
            ds_out.append(da_out.to_dataset(name=i))
  
    ds_out=xr.merge(ds_out)
    ds_out=ds_out.assign_coords(lon=('lon',lon1d))
    ds_out=ds_out.assign_coords(lat=('lat',lat1d))
   #ds_out=ds_out.assign_coords(lay=('lay',ds_in.Layer.values))
    ds_out=ds_out.assign_coords(lay=('lay',ds_in.zaxis_1.values))
    ds_out=ds_out.assign_coords(time=('time',ds_in.Time.values))
    ds_out['lon'].attrs['units']='degrees_east'
    ds_out['lon'].attrs['axis']='X'
    ds_out['lon'].attrs['standard_name']='longitude'
    ds_out['lat'].attrs['units']='degrees_north'
    ds_out['lat'].attrs['axis']='Y'
    ds_out['lat'].attrs['standard_name']='latitude'
    ds_out['lay'].attrs['units']='meters'
    ds_out['lay'].attrs['positive']='down'
    ds_out['lay'].attrs['axis']='Z'

    ds_out.to_netcdf(outfilename)
    ds_out.close()

    ds_in.close()

#------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  wgtdir = 'grids/'
  datadir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-ship/analysis/increment'
  outdir = './'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'wgtdir=', 'datadir=', 'outdir='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--wgtdir'):
      wgtdir = a
    elif o in ('--datadir'):
      datadir = a
    elif o in ('--outdir'):
      outdir = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('wgtdir = ', wgtdir)

  filename = '%s/xainc.20200110_030000z.nc4' %(datadir)
  files = []
  files.append(filename)

  o2ll = Ocean2LatLon(debug=debug, wgtdir=wgtdir)
  inres = 'mx100'
  outres='360x180'

  o2ll.set_ires(ires=inres)
  o2ll.set_ores(ores=outres)

  for infile in files:
    print('Processing:', infile)
    path, flnm = os.path.split(infile)
    if(flnm.find('.nc') > 0):
      namestr = '_%s.nc' %(outres)
      outfile = outdir + flnm.replace('.nc', namestr)
    else:
      outfile = '%s%s_%s.nc' %(outdir, flnm, outres)
    print('\toutput file:', outfile)

    o2ll.process_file(infilename=infile, outfilename=outfile)

