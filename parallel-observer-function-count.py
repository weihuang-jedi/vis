#########################################################################
#$Id: bld.py 28 2021-01-21 15:10:31Z whuang $
#$Revision: 28 $
#$HeadURL: file:///Users/whuang/.wei_svn_repository/trunk/jedi-build-tools/bld.py $
#$Date: 2021-01-21 08:10:31 -0700 (Thu, 21 Jan 2021) $
#$Author: whuang $
#########################################################################

import getopt
import os, sys
import subprocess
import time
import datetime

""" Counter """
class Counter:
  """ Constructor """
  def __init__(self, debug=0, workdir=None, caselist=[]):

    """ Initialize class attributes """
    self.debug = debug
    self.workdir = workdir
    self.caselist = caselist

    if(workdir is None):
      print('workdir not defined. Exit.')
      sys.exit(-1)

    if(len(caselist) < 1):
      print('caselist not defined. Exit.')
      sys.exit(-1)

  def get_filename(self, rundir):
    nf = 1
    has_more = True
    while(has_more):
      ftmp = '%s/log.getkf.%d' %(rundir, nf)
      nf += 1

      if(os.path.exists(ftmp)):
        flnm = ftmp
      else:
        has_more = False
    return flnm

  def process(self):
    self.node = 1
    self.core = 36
    self.parstatslist = []
    self.filelist = []
    nc = 0
    for case in self.caselist:
      rundir = '%s/%s/run_80.40t%dn_%dp' %(self.workdir, case,
                self.node, self.core)
     #flnm = '%s/stdoutNerr/stdout.00000000' %(rundir)
      flnm = self.get_filename(rundir)

      nc += 1
      if(os.path.exists(flnm)):
        if(self.debug):
          print('Case ' + str(nc) + ' name: ' + flnm)
       #pstats, gstats = self.stats(flnm)
        par_stats = self.stats(flnm)
        self.filelist.append(flnm)
        self.parstatslist.append(par_stats)
      else:
        print('Filename ' + flnm + ' does not exit. Stop')
        sys.exit(-1)

  def stats(self, flnm):
    if(os.path.exists(flnm)):
      pass
    else:
      print('Filename ' + flnm + ' does not exit. Stop')
      sys.exit(-1)

    print('flnm:', flnm)
    par_stats = {}
    with open(flnm) as fp:
      lines = fp.readlines()
     #line = fp.readline()
      num_lines = len(lines)
      print('Total number of lines: ', num_lines)

      nl = 0
      while(nl < num_lines):
        if(lines[nl].find('Object counts') > 0):
          nl, par_stats = self.object_counts(lines, nl)
          nl += num_lines
        nl += 1

    return par_stats

  def object_counts(self, lines, nl):
   #headleng = 11
    headleng = len('OOPS_STATS ')

    stats = {}
    going = 1
    ns = nl + 3
    while(going):
      line = lines[ns].strip()
      ns += 1
      if(line.find('Object counts') > 0):
        going = 0
        break

     #print('Line ' + str(ns) + ': ' + line)

      item = line.split(': ')
     #print('item=', item)
      namestr = item[0].strip()
     #print('namestr:', namestr)
      if(namestr.find('OOPS_STATS ') >= 0):
        name = namestr[headleng:]
      else:
        name = namestr

      cstr = item[1].strip()
      while(cstr.find('  ') > 0):
        cstr = cstr.replace('  ', ' ')
      clist = cstr.split(' ')
     #print('clist:', clist)

      cinfo = {}
      cinfo['Total'] = int(clist[0])
      cinfo['Simult'] = int(clist[1])

     #print('name: %s, Total: %d, Simult: %d' %(name, cinfo['Total'], cinfo['Simult']))

      stats[name] = cinfo

    return ns, stats

  def output_difference(self):
    txtname = 'object-count.csv'
    OTF = open(txtname, 'w')

    header = 'Object Counts\n'
    OTF.write(header)

    nc = len(self.caselist)
    header = '%-20s' %('Function Name')
    for i in range(nc):
      header = '%s, , %-20s' %(header, self.caselist[i])
    OTF.write(header+'\n')

    header = '%-20s' %('Object Counts')
    for i in range(nc):
      header = '%s, %-10s, %-10s'  %(header, 'Total', 'Simult')
    OTF.write(header+'\n')

   #Select common functions
    common_funcs = []
    namelist = self.parstatslist[0].keys()
    print('namelist =', namelist)
    for func in namelist:
      count = 1
      for n in range(1, nc):
        stats = self.parstatslist[n]
        if(func in stats.keys()):
          count += 1
      if(count == nc):
        common_funcs.append(func)

    print('common_funcs: ', common_funcs)

    stats0 = self.parstatslist[0]
    stats1 = self.parstatslist[1]
    for func in common_funcs:
      if(stats0[func]['Total'] == stats1[func]['Total'] and
         stats0[func]['Simult'] == stats1[func]['Simult']):
        continue
        
      txtinfo = '%-20s' %(func)
      for n in range(nc):
        stats = self.parstatslist[n]
        txtinfo = '%s, %10d, %10d' %(txtinfo, stats[func]['Total'], stats0[func]['Simult'])
      OTF.write(txtinfo+'\n')
    OTF.close()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  workdir = '/work2/noaa/gsienkf/weihuang/jedi/per_core_timing/run'
  caselist = ['anna_halo_aircraft', 'halo_aircraft',
              'anna_roundRobin_aircraft', 'roundRobin_aircraft']

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'workdir=',
                             'caselist='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--workdir'):
      workdir = a
    elif o in ('--corelist'):
      corelist = a
    elif o in ('--nodelist'):
      nodelist = a
    elif o in ('--casename'):
      casename = a
    elif o in ('--output'):
      output = a
    elif o in ('--linear'):
      linear = int(a)
    else:
      assert False, 'unhandled option'

  pr = Counter(debug=debug, workdir=workdir, caselist=caselist)
  pr.process()
  pr.output_difference()

