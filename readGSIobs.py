# =========================================================================
import os
import sys
import getopt
import numpy as np
import netCDF4

#=========================================================================
class ReadGSIobs():
  def __init__(self, debug=0, filename=None):
    self.debug = debug
    self.filename = filename

    if(self.debug):
      print('debug = ', debug)
      print('filename = ', filename)

  def set_filename(self, filename=None):
    self.filename = filename

  def get_latlon(self):
    ncfile = netCDF4.Dataset(self.filename, 'r')
    lat = ncfile.variables['Latitude'][:]
    lon = ncfile.variables['Longitude'][:]
    ncfile.close()

   #print('lat = ', lat)
   #print('lon = ', lon)

    return lat, lon

# ----
if __name__ == '__main__':
  debug = 1
  filename = '/work/noaa/gsienkf/weihuang/jedi/vis_tools/visfv3/jeff-all-obs/diag_conv_t_ges.2021010900_ensmean.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'yaml_file='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--filename'):
      filename = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('filename = ', filename)

  rgo = ReadGSIobs(debug=debug, filename=filename)
  lat, lon = rgo.get_latlon()
 
  print('lat = ', lat)
  print('lon = ', lon)

