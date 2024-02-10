#########################################################################
#$Id: bld.py 28 2021-01-21 15:10:31Z whuang $
#$Revision: 28 $
#$HeadURL: file:///Users/whuang/.wei_svn_repository/trunk/jedi-build-tools/bld.py $
#$Date: 2021-01-21 08:10:31 -0700 (Thu, 21 Jan 2021) $
#$Author: whuang $
#########################################################################

import getopt
import os, sys
import types
import time
import datetime
import subprocess
import netCDF4

import numpy as np


from readIODA2Obs import ReadIODA2Obs

#------------------------------------------------------------------
class Compare2ObsFiles():
  def __init__(self, debug=0):
    self.debug = debug

    self.delt = 1.0e-10

#------------------------------------------------------------------
  def setup(self, datadir=None, filename=None):
    if(self.debug):
      print('datadir = ', datadir)
      print('filename = ', filename)

    basefile = '%s/hofx.nl/%s' %(datadir, filename)
    if(os.path.exists(basefile)):
      self.baseio = ReadIODA2Obs(debug=debug, filename=basefile)
    else:
      print('basefile: <%s> does not exist. Stop.' %(basefile))
      sys.exit(-1)

    casefile = '%s/hofx.lo/%s' %(datadir, filename)
    if(os.path.exists(casefile)):
      self.caseio = ReadIODA2Obs(debug=debug, filename=casefile)
    else:
      print('casefile: <%s> does not exist. Stop.' %(casefile))
      sys.exit(-1)

#------------------------------------------------------------------
  def get_groups(self):
    basegrps = self.baseio.get_groups()
    casegrps = self.caseio.get_groups()
    grps = []
    for grp in casegrps:
      if(grp in basegrps):
        if(grp != 'MetaData'):
          grps.append(grp)

    return grps

#------------------------------------------------------------------
  def get_variables(self, grp):
    return self.caseio.groups[grp].variables

#------------------------------------------------------------------
  def compare_int_variable(self, varname):
    base_var = self.baseio.get_var(varname)
    case_var = self.caseio.get_var(varname)

   #print('len(case_var) = ', len(case_var))
   #print('case_var) = ', case_var)
   #print('len(base_var) = ', len(base_var))
   #print('base_var) = ', base_var)

    diff = case_var - base_var
    nd = 0
    for n in range(len(case_var)):
      if(diff[n] != 0):
        nd += 1

    if(nd):
      print('Difference of %s\n' %(varname))
      print('\tLen(var): %d, nd: %d, mindif: %f, maxdif: %f' %(len(case_var), nd, np.min(diff), np.max(diff)))
      for n in range(len(case_var)):
        if(diff[n] != 0):
          print('\tNo %d: c %d, b %d, diff: %d' %(n+1, case_var[n], base_var[n], diff[n]))
    else:
      print('var: %s has no difference' %(varname))

#------------------------------------------------------------------
  def compare_float_variable(varname):
    base_var = self.baseio.get_var(varname)
    case_var = self.caseio.get_var(varname)

   #print('len(case_var) = ', len(case_var))
   #print('case_var) = ', case_var)
   #print('len(base_var) = ', len(base_var))
   #print('base_var) = ', base_var)

    diff = case_var - base_var

    nd = 0
    for n in range(len(case_var)):
      if(abs(diff[n]) > self.delt):
        nd += 1

    if(nd):
      print('Difference of %s\n' %(varname))
      print('\tLen(var): %d, nd: %d, mindif: %f, maxdif: %f' %(len(case_var), nd, np.min(diff), np.max(diff)))
      for n in range(len(case_var)):
        if(abs(diff[n]) > self.delt):
          print('\tNo %d: c %f, b %f, diff: %f' %(n+1, case_var[n], base_var[n], diff[n]))
    else:
      print('var: %s has no difference' %(varname))
 
#------------------------------------------------------------------
  def compare_2d_int_variable(varname):
    base_var = self.baseio.get_2d_var(varname)
    case_var = self.caseio.get_2d_var(varname)

   #print('len(case_var) = ', len(case_var))
   #print('case_var) = ', case_var)
   #print('len(base_var) = ', len(base_var))
   #print('base_var) = ', base_var)

    diff = case_var - base_var
    nd = 0
    for j in range(len(case_var)):
      for i in range(len(case_var[0])):
        if(diff[j,i] != 0):
          nd += 1

    if(nd):
      print('Difference of %s\n' %(varname))
      print('\tLen(var): %d, nd: %d, mindif: %f, maxdif: %f' %(len(case_var), nd, np.min(diff), np.max(diff)))
      for j in range(len(case_var)):
        for i in range(len(case_var[0])):
          if(diff[j,i] != 0):
            print('\tJ: %d, I: %d: c %d, b %d, diff: %d' %(j, i, case_var[j,i], base_var[j,i], diff[j,i]))
    else:
      print('var: %s has no difference' %(varname))

#------------------------------------------------------------------
  def compare_2d_float_variable(datadir, filename, varname):
    base_var = self.baseio.get_2d_var(varname)
    case_var = self.caseio.get_2d_var(varname)

    diff = case_var - base_var

   #print('len(case_var) = ', len(case_var))
   #print('case_var) = ', case_var)
   #print('len(base_var) = ', len(base_var))
   #print('base_var) = ', base_var)

   #print('\tbase_var.max: %f, base_var.min: %f' %(np.max(base_var[:,0]), np.min(base_var[:,0])))
   #print('\tcase_var.max: %f, case_var.min: %f' %(np.max(case_var[:,0]), np.min(case_var[:,0])))
   #print('\diffvar.max: %f, diffvar.min: %f' %(np.max(diffvar), np.min(diffvar)))

    nd = 0
    for j in range(len(case_var)):
      for i in range(len(case_var[0])):
        if(abs(diff[j,i]) > delt):
          nd += 1

    if(nd):
      print('Difference of %s\n' %(varname))
      print('\tbase_var.max: %f, base_var.min: %f' %(np.max(base_var[:,0]), np.min(base_var[:,0])))
      print('\tcase_var.max: %f, case_var.min: %f' %(np.max(case_var[:,0]), np.min(case_var[:,0])))
      print('\tLen(var): %d, nd: %d, mindif: %f, maxdif: %f' %(len(case_var), nd, np.min(diff), np.max(diff)))
      for j in range(len(case_var)):
        for i in range(len(case_var[0])):
          if(abs(diff[j,i]) > delt):
            print('\tJ: %d, I: %d: c %d, b %d, diff: %d' %(j, i, case_var[j,i], base_var[j,i], diff[j,i]))
    else:
      print('var: %s has no difference' %(varname))

#------------------------------------------------------------------
if __name__== '__main__':
  debug = 0
  datadir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data'
 #filename = 'sfc_letkf-gfs_2020121500_m.nc4'
  filename = 'sfc_letkf-gfs_2020121500_m.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'filename='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--filename'):
      filename = a
    else:
      print('option: ', a)
      assert False, 'unhandled option'

#------------------------------------------------------------------
  cof = Compare2ObsFiles()

#------------------------------------------------------------------
  for filename in ['sfc_letkf-gfs_2020121500_m.nc4',
                   'aircraft_letkf-gfs_2020121500_m.nc4',
                   'scatwind_letkf-gfs_2020121500_m.nc4']:
    cof.setup(datadir=datadir, filename=filename)
    varlist = cof.get_variables('EffectiveQC0')
    for grp in ctof.get_groups():
      for var in varlist:
        varname = '/%s/%s' %(grp, var)
        if(grp in ['EffectiveQC0', 'EffectiveQC2', 'ObsType', 'PreQC', 'PreUseFlag']):
          cof.compare_int_variable(varname)
        else:
          cof.compare_float_variable(varname)

#------------------------------------------------------------------
  filename = 'amsua_n19_letkf-gfs_2020121500_m.nc4'
  cof.setup(datadir=datadir, filename=filename)
  varlist = cof.get_variables('EffectiveQC0')
  for grp in ctof.get_groups():
    for var in varlist:
      varname = '/%s/%s' %(grp, var)
      if(grp in ['EffectiveQC0', 'EffectiveQC2', 'ObsType', 'PreQC', 'PreUseFlag']):
        cof.compare_2d_int_variable(varname)
      else:
        cof.compare_2d_float_variable(varname)

