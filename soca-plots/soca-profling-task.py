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
  def __init__(self, debug=0, taskpernode=30, output=0,
               nodelist=[], workdir=None, linear=0):

    """ Initialize class attributes """
    self.debug = debug
    self.output = output
    self.workdir = workdir
    self.taskpernode = taskpernode
    self.nodelist = nodelist
    self.linear = linear

    if(workdir is None):
      print('workdir not defined. Exit.')
      sys.exit(-1)

    if(taskpernode < 1):
      print('taskpernode not defined. Exit.')
      sys.exit(-1)

    if(len(nodelist) < 1):
      print('nodelist not defined. Exit.')
      sys.exit(-1)

    self.colorlist = ['red', 'blue', 'green', 'orange', 'royalblue', 'cyan', 'magenta']
    self.function_list = ['total',
                          'computeHofX',
                          'GetValues',
                          'LocalInterpolator',
                          'measurementUpdate',
                          'ObsSpace.read()',
                          'State.read()']
    self.fullfunction_list = ['util::Timers::Total',
                              'oops::LocalEnsembleSolver::computeHofX',
                              'oops::GetValues::GetValues',
                              'oops::LocalInterpolator::LocalInterpolator',
                              'oops::LETKFSolver::measurementUpdate',
                              'oops::ObsSpace::ObsSpace',
                              'oops::State::State']
    self.pmin = 999999999999.0
    self.gmin = 999999999999.0
    self.pmax = 0.0
    self.gmax = 0.0

  def set_output(self, output=0):
    self.output = output

  def set_linear(self, linear=1):
    self.linear = linear

  def process(self):
    self.pstatslist = []
    self.gstatslist = []

    self.filelist = []
    self.corelist = []
    for node in self.nodelist:
      totalcpus = self.taskpernode*node
      rundir = '%s/soca_solver.%dt%dn_%dp' %(self.workdir, self.taskpernode, node, totalcpus)
      flnm = '%s/stdoutNerr/stdout.00000000' %(rundir)

      if(os.path.exists(flnm)):
       #if(self.debug):
       #  print('Case ' + str(nc) + ' name: ' + flnm)
        self.filelist.append(flnm)
        self.corelist.append(totalcpus)
        if(self.debug):
          print('Processing task: %d, node: %d, as file: %s' %(self.taskpernode, node, flnm))
       #pstats, gstats = self.stats(flnm)
       #self.pstatslist.append(pstats)
       #self.gstatslist.append(gstats)
        ptime = self.stats(flnm)
        self.pstatslist.append(ptime)
      else:
        print('Filename ' + flnm + ' does not exit. Stop')
        sys.exit(-1)

  def stats(self, flnm):
    if(os.path.exists(flnm)):
      pass
    else:
      print('Filename ' + flnm + ' does not exit. Stop')
      sys.exit(-1)

    prof = {}

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
          nl, par_stats = self.parallel_time_stats(lines, nl)
          nl += num_lines
       #elif(lines[nl].find('Timing Statistics') > 0):
       #  if(self.debug):
       #    print('Start Timing Statistics')
       #  nl, gen_stats = self.time_stats(lines, nl)
        nl += 1

   #return par_stats, gen_stats

    avgtime = []

    for n in range(len(self.fullfunction_list)):
      name = self.fullfunction_list[n]
      idx = self.get_index(par_stats, name)
      if(idx < 0):
        avgt = 0.0
      else:
        avgt = par_stats[idx]['avg']

      avgtime.append(avgt)

    return avgtime

  def get_index(self, stats, varname):
    idx = -1
    for n in range(len(stats)):
      if(varname == stats[n]['name']):
        idx = n
        break
    return idx

  def time_stats(self, lines, nl):
    stats = []
    going = 1
    ns = nl + 2
    while(going):
      line = lines[ns].strip()

      ns += 1
      
      if(line.find('Timing Statistics') > 0):
        going = 0
        break

     #print('Line ' + str(ns) + ': ' + line)

      item = line.split(': ')
     #print(item)
      nlist = item[0].strip().split(' ')
      name = nlist[1]

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      nlist = tstr.split(' ')
      ft = float(nlist[0])
     #if(ft < 1.0):
     #  continue

      tinfo = {}
      tinfo['name'] = name
      tinfo['time'] = ft
      tinfo['call'] = int(nlist[2])

      stats.append(tinfo)

      if(self.gmin > ft):
        self.gmin = ft
      if(self.gmax < ft):
        self.gmax = ft

    return ns, stats

  def parallel_time_stats(self, lines, nl):
    stats = []
    going = 1
    ns = nl + 3
    while(going):
      line = lines[ns].strip()
      ns += 1
      if(line.find('Parallel Timing Statistics') > 0):
        going = 0
        break

     #print('Line ' + str(ns) + ': ' + line)

      item = line.split(': ')
     #print(item)
      nlist = item[0].strip().split(' ')
      name = nlist[1]

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      nlist = tstr.split(' ')
      ft = float(nlist[0])
     #if(ft < 1.0):
     #  continue

      tinfo = {}
      tinfo['name'] = name
      tinfo['min'] = ft
      tinfo['max'] = float(nlist[1])
      tinfo['avg'] = float(nlist[2])

     #total = 100. * avg / total
     #tinfo['total'] = float(nlist[3])
   
     #imbalance = 100. * (max - min) / avg
     #tinfo['imbalance'] = float(nlist[4])

      stats.append(tinfo)

      if(self.pmin > ft):
        self.pmin = ft
      if(self.pmax < ft):
        self.pmax = ft

    return ns, stats

  def plot(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = 'Top Functions Timing of %d Cores per Node' %(self.taskpernode)

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    labels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      labels.append(lbl)

   #print('x = ', x)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.xscale('linear')
    else:
     #plt.xscale('linear')
      plt.xscale('log', base=2)
      plt.yscale('log', base=10)

      plt.xticks(x, labels)
     #plt.xticks(x, labels, rotation ='vertical')

    pmin = 1.0e20
    pmax = 0.0

    for i in range(len(self.fullfunction_list)):
      for k in range(nl):
        y[k] = 0.001*self.pstatslist[k][i]
        if(pmin > y[k]):
          pmin = y[k]
        if(pmax < y[k]):
          pmax = y[k]
     #print('y = ', y)
      ax.plot(x, y, color=self.colorlist[i], linewidth=2, alpha=0.9, label=self.function_list[i])

    if(self.linear == 0):
      for i in range(len(self.fullfunction_list)):
        for k in range(nl):
          fact = 1.0/np.log2(self.nodelist[k])
          z[k] = 0.001*self.pstatslist[0][i]*fact
       #https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
        ax.plot(x, z, color='black', linewidth=1, alpha=0.5, linestyle='dotted')

    plt.grid()

   #Same limits for everybody!
    print('pmin: %f, pmax: %f' %(pmin, pmax))
    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
   #plt.suptitle(title, fontsize=13, fontweight=0, color='black', style='italic', y=1.02)
    plt.suptitle(title, fontsize=15, fontweight=0, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.2, right=0.70, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (second)', labelpad=10)

   #Create the legend
    fig.legend(ax, labels=self.function_list,
           loc='center right',   # Position of legend
           title='Funcs Legend', # Title for the legend
           fontsize=8,
           borderpad=1.2,
           labelspacing=1.2,
           handlelength=1.5
           )

    if(self.linear):
      imgname = 'lin_%dtasks.png' %(self.taskpernode)
    else:
      imgname = 'log_%dtasks.png' %(self.taskpernode)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
  workdir = '/work2/noaa/gsienkf/weihuang/ufs/soca/new-soca-solver'
  taskpernode = 30
  nodelist = [2, 4, 6, 8, 10, 12]
  linear = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'workdir=',
                             'taskpernode=', 'nodelist='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--workdir'):
      workdir = a
    elif o in ('--taskpernode'):
      taskpernode = int(a)
    elif o in ('--linear'):
      linear = int(a)
    elif o in ('--output'):
      output = int(a)
    else:
      assert False, 'unhandled option'

  pr = Profiler(debug=debug, taskpernode=taskpernode, nodelist=nodelist,
                output=output, workdir=workdir, linear=linear)
  pr.process()
  for linear in [0, 1]:
    pr.set_linear(linear=linear)
    for output in [0, 1]:
      pr.set_output(output=output)
      pr.plot()

