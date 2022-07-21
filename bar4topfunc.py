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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def cmdout(command):
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
    ostr = result.stdout
    return ostr.strip()

""" Profiler """
class Profiler:
  """ Constructor """
  def __init__(self, debug=0, corelist=[], filelist=[], linear=0, tpn=20, show=0):

    """ Initialize class attributes """
    self.debug = debug
    self.corelist = corelist
    self.filelist = filelist
    self.linear = linear
    self.tpn = tpn
    self.show = show

    self.colorlist = ['red', 'blue', 'green', 'orange', 'royalblue', 'cyan']
    self.function_list = ['total', 'computeHofX', 'GetValues', 'LocalInterpolator', 'ObsSpace', 'State']
    self.fullfunction_list = ['util::Timers::Total', 'oops::LocalEnsembleSolver::computeHofX',
                              'oops::GetValues::GetValues', 'oops::LocalInterpolator::LocalInterpolator',
                              'oops::ObsSpace::ObsSpace', 'oops::State::State']
  def set_core_file(self, corelist=[], filelist=[]):
    self.corelist = corelist
    self.filelist = filelist

  def process(self):
    self.stats_list = {}
   #print('self.corelist = ', self.corelist)
   #print('self.filelist = ', self.filelist)
    for n in range(len(self.corelist)):
      core = self.corelist[n]
      flnm = self.filelist[n]

      if(os.path.exists(flnm)):
        print('Processing file: ', flnm)
        self.stats_list[core] = self.stats(flnm)
      else:
        print('File <' + flnm + '> does not exist.')

  def stats(self, flnm):
    avgtime = {}

    with open(flnm) as fp:
      lines = fp.readlines()
     #line = fp.readline()
      num_lines = len(lines)
     #print('Total number of lines: ', num_lines)

      nl = 0
      while(nl < num_lines):
        if(lines[nl].find('Parallel Timing Statistics') > 0):
         #if(self.debug):
         #  print('Start Parallel Timing Statistics')
          nl, avgtime = self.parallel_time_stats(lines, nl)
          nl += num_lines
        nl += 1

    return avgtime

  def parallel_time_stats(self, lines, nl):
    avgtime = {}
    going = 1
    ns = nl + 3
    while(going):
      line = lines[ns].strip()
      ns += 1
      if(line.find('Parallel Timing Statistics') > 0):
        going = 0
        break

     #print('Line ' + str(ns) + ': ' + line)

      item = line.split(' : ')
     #print(item)
      nlist = item[0].strip().split(' ')
      name = nlist[1]

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      nlist = tstr.split(' ')

      for i in range(len(self.fullfunction_list)):
        funcname = self.fullfunction_list[i]
        if(name == funcname):
          avgtime[i] = float(nlist[2])
         #print('      ' + name + ':' + nlist[2])

    return ns, avgtime

  def plot(self):
    title = 'Top function timing of %d Cores per Node' %(self.tpn)

    x = np.array(self.function_list)
    ts = {}
    line = []

    ymin = 1.0e+36
    ymax = 0.0

    nc = len(self.corelist)
    nf = len(self.function_list)
   #print('nc = ', nc)
   #print('nf = ', nf)

    for n in range(nc):
      core = self.corelist[n]
      ts[n] = []
     #print('core = ', core)
     #print('self.stats_list[core] = ', self.stats_list[core])
      for i in range(nf):
        t = 0.001*self.stats_list[core][i]
        ts[n].append(t)
        if(t > ymax):
          ymax = t
        if(t < ymin):
          ymin = t

   #print('x = ', x)
   #print('ymin = ', ymin)
   #print('ymax = ', ymax)

    if(self.linear):
      tmin = 10.0*int(ymin/10.0) - 10.0
      tmax = 10.0 + 10.0*int(ymax/10.0)
    else:
      tmin = 8192.0
      while(tmin > ymin):
        tmin /= 2.0
      tmax = 1.0
      while(tmax < ymax):
        tmax *= 2.0

    plt.xlim([x[0], x[-1]])
    plt.ylim([tmin, tmax])

    nf = len(self.function_list)
    nc = len(self.corelist)

    labellist = []
    for core in self.corelist:
      label = str(core)
      labellist.append(label)

    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

   #the actual graph:
    fig, ax = plt.subplots(figsize = (10,4))

    idx = np.asarray([i for i in range(nf)])

    width = 1.0/(nf+2)

    for n in range(nc):
      ax.bar(idx+n*width, ts[n], width, color=self.colorlist[n])

    ax.set_xticks(idx)
    ax.set_xticklabels(self.function_list, rotation=65)
    ax.legend(labellist)
    ax.set_xlabel('Functions')
    ax.set_ylabel('Time (second)')
    ax.set_title(title)

    fig.tight_layout()

    if(self.linear):
      imgname = 'lin_bar_%dtpn.png' %(self.tpn)
    else:
      imgname = 'log_bar_%dtpn.png' %(self.tpn)
      plt.yscale("log")
     #plt.yscale("log", base=2)
     #plt.yscale("log", base=10)

    if(self.show):
      plt.show()
    else:
      plt.savefig(imgname)

    plt.cla()
    plt.clf()

  def set_tpn(self, tpn=20):
    self.tpn = tpn

  def set_linear(self, linear=1):
    self.linear = linear

  def set_show(self, show=1):
    self.show = show

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  linear = 1
  tpn = 20
  show = 1

  rundir = '/work2/noaa/gsienkf/weihuang/ufs/soca/new-soca-solver'
  tasklist = [20, 24, 30, 36, 40]
 #tasklist = [20]
  nodelist = [2, 4, 6, 8, 10, 12]
  filelist = []
  corelist = []

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'corelist=', 'filelist='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--corelist'):
      corelist = a
    elif o in ('--filelist'):
      filelist = a
    elif o in ('--linear'):
      linear = int(a)
    elif o in ('--tpn'):
      tpn = int(a)
    elif o in ('--show'):
      show = int(a)
    else:
      assert False, 'unhandled option'

  pr = Profiler(debug=debug, corelist=corelist, filelist=filelist, tpn=tpn, show=show)

  for task in tasklist:
    filelist = []
    corelist = []
    for node in nodelist:
      totalcpus = task*node
      workdir = '%s/soca_solver.%dt%dn_%dp' %(rundir, task, node, totalcpus)
      workdir = '%s/stdoutNerr/stdout.00000000' %(workdir)
      filelist.append(workdir)
      corelist.append(totalcpus)

    pr.set_core_file(corelist=corelist, filelist=filelist)
    pr.process()
    pr.set_tpn(tpn=task)

    for show in [1, 0]:
   #for show in [0]:
      pr.set_show(show=show)
      for linear in [1, 0]:
        pr.set_linear(linear=linear)
        pr.plot()

