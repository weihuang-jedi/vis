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
  def __init__(self, debug=0, tasklist=[], output=0,
               nodelist=[], workdir=None, linear=0):

    """ Initialize class attributes """
    self.debug = debug
    self.output = output
    self.workdir = workdir
    self.tasklist = tasklist
    self.nodelist = nodelist
    self.linear = linear

    if(workdir is None):
      print('workdir not defined. Exit.')
      sys.exit(-1)

    if(len(tasklist) < 1):
      print('tasklist not defined. Exit.')
      sys.exit(-1)

    if(len(nodelist) < 1):
      print('nodelist not defined. Exit.')
      sys.exit(-1)

    self.colorlist = ['red', 'blue', 'green', 'orange', 'royalblue', 'cyan']
    self.function_list = ['total', 'computeHofX', 'GetValues', 'LocalInterpolator', 'ObsSpace', 'State']
    self.fullfunction_list = ['util::Timers::Total', 'oops::LocalEnsembleSolver::computeHofX',
                              'oops::GetValues::GetValues', 'oops::LocalInterpolator::LocalInterpolator',
                              'oops::ObsSpace::ObsSpace', 'oops::State::State']
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
    for task in self.tasklist:
      flist = []
      clist = []
      plist = []
     #glist = []
      for node in self.nodelist:
        totalcpus = task*node
        rundir = '%s/soca_solver.%dt%dn_%dp' %(self.workdir, task, node, totalcpus)
       #flnm = '%s/stdoutNerr/stdout.00000000' %(rundir)
        flnm = '%s/log.soca_solver' %(rundir)

        if(os.path.exists(flnm)):
         #if(self.debug):
         #  print('Case ' + str(nc) + ' name: ' + flnm)
          flist.append(flnm)
          clist.append(totalcpus)
          if(self.debug):
            print('Processing task: %d, node: %d, as file: %s' %(task, node, flnm))
         #pstats, gstats = self.stats(flnm)
         #plist.append(pstats)
         #glist.append(gstats)
          ptime = self.stats(flnm)
          plist.append(ptime)
        else:
          print('Filename ' + flnm + ' does not exit. Stop')
          sys.exit(-1)

      self.filelist.append(flist)
      self.corelist.append(clist)

      self.pstatslist.append(plist)
     #self.gstatslist.append(glist)

  def print_gen_time(self, nc):
    print('\nTiming Statistics: Time')
    nl = len(self.stats_list[0][1])

    hinfo = '+----+'
    for i in range(nc):
      hinfo = '%s%s+' %(hinfo, 22*'-')
    print(hinfo)

    pinfo = '+-NP-+'
    for i in range(nc):
      pinfo = '%s%s%4d%s+' %(pinfo, 9*' ', self.caselist[i], 9*' ')
    print(pinfo)

    print(hinfo)

    n = 0
    for n in range(nl-10):
      ni = n + 1

      pinfo = '|%3d |' %(ni)
      for i in range(nc):
       #print('self.stats_list[i][1] = ', self.stats_list[i][1])
        idx = self.stats_list[i][1][n]
        tim = self.stats_list[i][0][idx]['time']
        pinfo = '%s %s %10.2f %s |' %(pinfo, 4*' ', tim, 4*' ')
      print(pinfo)
    print(hinfo)

  def print_par_time(self, nc):
    print('\nParallel Timing Statistics: Time')
    nl = len(self.stats_list[0][3])

    hinfo = '+----+'
    for i in range(nc):
      hinfo = '%s%s+' %(hinfo, 22*'-')
    print(hinfo)

    pinfo = '+-NP-+'
    for i in range(nc):
      pinfo = '%s%s%4d%s+' %(pinfo, 9*' ', self.caselist[i], 9*' ')
    print(pinfo)

    print(hinfo)

    n = 0
    for n in range(nl-10):
      ni = n + 1

      pinfo = '|%3d |' %(ni)
      for i in range(nc):
       #print('self.stats_list[i][3] = ', self.stats_list[i][3])
        idx = self.stats_list[i][3][n]
        tim = self.stats_list[i][2][idx]['avg']
        pinfo = '%s %s %10.2f %s |' %(pinfo, 4*' ', tim, 4*' ')
      print(pinfo)
    print(hinfo)

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

  def plot_panel_by_task(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

   #title = 'Timing by Task'
    title = 'Timing by Task (separate reinit obs)'

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

   #print('x = ', x)

    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    ax = axes.flatten()

    pmin = 1.0e20
    pmax = 0.0

    for n in range(len(self.tasklist)):
      for i in range(len(self.fullfunction_list)):
        for k in range(nl):
          y[k] = 0.001*self.pstatslist[n][k][i]
          if(pmin > y[k]):
            pmin = y[k]
          if(pmax < y[k]):
            pmax = y[k]
       #print('y = ', y)
        ax[n].plot(x, y, color=self.colorlist[i], linewidth=2, alpha=0.9)

      label = '%d tasks per node' %(self.tasklist[n])

     #label physical distance to the left and up:
      ax[n].set_title(label, fontsize=8)

   #Adjust pmin, pmax
    pvmin = 1.0
    while(pvmin > pmin):
      pvmin /= 2.0
    pmin = pvmin
    pvmax = 1.0
    while(pvmax < pmax):
      pvmax *= 2.0
    pmax = pvmax

    ylabels = []
    ypoints = []
    pv = pvmin
    while(pv <= pvmax):
      ypoints.append(pv)
      lbl = '%6.1f' %(pv)
      ylabels.append(lbl)
      pv *= 2.0

    if(self.linear):
      plt.xscale('linear')
    else:
     #plt.xscale('linear')
     #plt.yscale('log', base=10)
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)

      plt.xticks(x, xlabels)
     #plt.xticks(x, xlabels, rotation ='vertical')
      plt.yticks(ypoints, ylabels)

   #Same limits for everybody!
    print('pmin: %f, pmax: %f' %(pmin, pmax))
    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)

    plt.grid()
 
   #general title
   #plt.suptitle(title, fontsize=13, fontweight=0, color='black', style='italic', y=1.02)
    plt.suptitle(title, fontsize=15, fontweight=0, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.1, right=0.75, top=0.9)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (second)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=self.function_list,
           loc='center right',   # Position of legend
           title='Funcs Legend', # Title for the legend
           fontsize=8
           )
#          borderpad=1.2,
#          labelspacing=1.2,
#          handlelength=1.5

    if(self.linear):
     #imgname = 'lin_panel_by_task.png'
      imgname = 'separate_reinit_obs_lin_panel_by_task.png'
    else:
     #imgname = 'log_panel_by_task.png'
      imgname = 'separate_reinit_obs_log_panel_by_task.png'

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

  def plot_panel_by_node(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

   #title = 'Timing by Node'
    title = 'Timing by Node (separate reinit obs)'

    nl = len(self.tasklist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    for k in range(nl):
      x[k] = self.tasklist[k]

   #print('x = ', x)

    fig, axes = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=True)
    ax = axes.flatten()

    if(self.linear):
      plt.xscale('linear')
    else:
      plt.xscale('linear')
     #plt.xscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)

    pmin = 1.0e20
    pmax = 0.0

    for n in range(len(self.nodelist)):
      for i in range(len(self.fullfunction_list)):
        for k in range(nl):
          y[k] = 0.001*self.pstatslist[k][n][i]
          if(pmin > y[k]):
            pmin = y[k]
          if(pmax < y[k]):
            pmax = y[k]
       #print('y = ', y)
        ax[n].plot(x, y, color=self.colorlist[i], linewidth=2, alpha=0.9)

      label = '%d nodes' %(self.nodelist[n])

     #label physical distance to the left and up:
      ax[n].set_title(label, fontsize=8)

   #Same limits for everybody!
    print('pmin: %f, pmax: %f' %(pmin, pmax))
    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #Add title
   #fig.title(title, loc='Center', fontsize=12, fontweight=0)
   #plt.legend()
 
   #general title
   #plt.suptitle(title, fontsize=13, fontweight=0, color='black', style='italic', y=1.02)
    plt.suptitle(title, fontsize=15, fontweight=0, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.subplots_adjust(bottom=0.1, right=0.75, top=0.9)

    bs.set_xlabel('Task', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (second)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=self.function_list,
           loc='center right',   # Position of legend
           title='Funcs Legend', # Title for the legend
           fontsize=8,
           borderpad=1.2,
           labelspacing=1.2,
           handlelength=1.5
           )

   #Adjust the scaling factor to fit your legend text completely outside the plot
   #(smaller value results in more space being made for the legend)

    if(self.linear):
     #imgname = 'lin_panel_by_node.png'
      imgname = 'separate_reinit_obs_lin_panel_by_node.png'
    else:
     #imgname = 'log_panel_by_node.png'
      imgname = 'separate_reinit_obs_log_panel_by_node.png'
    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  output = 0
 #workdir = '/work2/noaa/gsienkf/weihuang/ufs/soca/new-soca-solver'
  workdir = '/work2/noaa/gsienkf/weihuang/jedi/run.soca'
 #tasklist = [20, 24, 30, 36, 40]
 #tasklist = [20, 30, 36, 40]
 #tasklist = [20, 32, 36, 40]
  tasklist = [32, 36, 40]
  nodelist = [2, 4, 6, 8, 10, 12]
  linear = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'workdir=', 'output=',
                             'tasklist=', 'nodelist='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--workdir'):
      workdir = a
    elif o in ('--tasklist'):
      tasklist = a
    elif o in ('--linear'):
      linear = int(a)
    elif o in ('--output'):
      output = int(a)
    else:
      assert False, 'unhandled option'

  pr = Profiler(debug=debug, tasklist=tasklist, nodelist=nodelist,
                output=output, workdir=workdir, linear=linear)
  pr.process()
  for linear in [0, 1]:
    pr.set_linear(linear=linear)
    for output in [0, 1]:
      pr.set_output(output=output)
      pr.plot_panel_by_task()
      pr.plot_panel_by_node()

