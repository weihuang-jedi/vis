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
  def __init__(self, debug=0, casename=None, output=0,
               corepernode=40, nodelist=[], workdir=None, linear=0):

    """ Initialize class attributes """
    self.debug = debug
    self.workdir = workdir
    self.nodelist = nodelist
    self.casename = casename
    self.output = output
    self.linear = linear
    self.corepernode = corepernode

    self.pvmax = 256.0
    self.pvmin = 0.125

    if(workdir is None):
      print('workdir not defined. Exit.')
      sys.exit(-1)

    if(corepernode < 1):
      print('corepernode not defined. Exit.')
      sys.exit(-1)

    if(len(nodelist) < 1):
      print('nodelist not defined. Exit.')
      sys.exit(-1)

    if(not os.path.exists(casename)):
     #mode
     #mode = 0o722
     #os.makedirs(casename, mode)
      os.makedirs(casename)

    self.colorlist = ['red', 'blue', 'green', 'orange', 'royalblue', 'cyan',
                      'magenta', 'lime', 'yellowgreen',
                      'violet', 'navy', 'teal']

   #self.get_top_functions()

  def set_minmax(self, vmin=0.125, vmax=256.0):
    self.pvmin = vmin
    self.pvmax = vmax

  def set_linear(self, linear=1):
    self.linear = linear

  def set_output(self, output=1):
    self.output = output

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

  def get_top_functions(self, maxfuncs=6):
    self.max_selected_functions = maxfuncs
    self.num_selected_functions = 0
    self.selected_function_list = []
    self.selected_functions_time = []
    self.selected_functions_name = []
    par_stats = self.parstatslist[0]

    for name in par_stats.keys():
      avgt = par_stats[name]['avg']
      self.add2selected_functions(name, avgt)

    for n in range(self.num_selected_functions):
      pinfo = 'No. %3.3d name: %-40s' %(n, self.selected_functions_name[n])
      pinfo = '%s, time: %8.2f' %(pinfo, self.selected_functions_time[n])
      print(pinfo)

  def add2selected_functions(self, name, avgt):
    n = self.num_selected_functions - 1
    if(self.num_selected_functions < self.max_selected_functions):
      self.selected_functions_time.append(avgt)
      self.selected_functions_name.append(name)
      self.num_selected_functions += 1
    else:
      if(avgt < self.selected_functions_time[n]):
        return
      self.selected_functions_time[n] = avgt
      self.selected_functions_name[n] = name

    while(n > 0):
      if(self.selected_functions_time[n] > self.selected_functions_time[n-1]):
        otime = self.selected_functions_time[n-1]
        oname = self.selected_functions_name[n-1]
        self.selected_functions_time[n-1] = self.selected_functions_time[n]
        self.selected_functions_name[n-1] = self.selected_functions_name[n]
        self.selected_functions_time[n] = otime
        self.selected_functions_name[n] = oname
      n -= 1

  def process(self):
    self.parstatslist = []

    self.filelist = []
    for n in range(len(self.nodelist)):
      totalcores = self.corepernode*self.nodelist[n]
      rundir = '%s/soca_solver.%dt%dn_%dp' %(self.workdir, self.corepernode,
                self.nodelist[n], totalcores)
      flnm = '%s/log.soca_solver' %(rundir)
     #flnm = self.get_filename(rundir)

      if(os.path.exists(flnm)):
       #if(self.debug):
       #  print('Case ' + str(nc) + ' name: ' + flnm)
        if(self.debug):
          print('Processing node: %d, as file: %s' %(self.nodelist[n], flnm))
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

    par_stats = {}
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
        nl += 1

    return par_stats

  def parallel_time_stats(self, lines, nl):
   #headleng = 11
    headleng = len('OOPS_STATS ')

    stats = {}
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
      namestr = item[0].strip()
     #print(namestr)
      if(namestr.find('OOPS_STATS ') >= 0):
        name = namestr[headleng:]
      else:
        name = namestr

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      tlist = tstr.split(' ')

      tinfo = {}
      tinfo['min'] = float(tlist[0])
      tinfo['max'] = float(tlist[1])
      tinfo['avg'] = float(tlist[2])
     #tinfo['percent'] = float(tlist[3])
     #tinfo['imbalance'] = float(tlist[4])

     #print('name: %s, avg: %f' %(name, tinfo['avg']))

      stats[name] = tinfo

    return ns, stats

  def get_minmax(self, statstime):
    pmin = statstime[0][0]
    pmax = pmin

    il = len(self.funclist)
    kl = len(self.nodelist)
    for k in range(kl):
      for i in range(il):
        if(pmin > statstime[i][k]):
          pmin = statstime[i][k]
        if(pmax < statstime[i][k]):
          pmax = statstime[i][k]

    nm = 0
    pvmin = 1.0
    while((pvmin > pmin) and (nm < 5)):
      pvmin *= 0.5
      nm +=1

    pvmax = 1.0
    while(pvmax < pmax):
      pvmax *= 2.0

   #pvmin = 0.125
   #pvmax = 256.0

    return pmin, pmax, pvmin, pvmax

  def get_main_statstime(self, funclabels, funclist):
    self.funclabels = funclabels
    self.funclist = funclist
 
    il = len(self.funclist)
    kl = len(self.nodelist)
    statstime = np.zeros((il, kl))
    for k in range(kl):
      stats = self.parstatslist[k]
      for i in range(il):
        name = self.funclist[i]
        if(name in stats.keys()):
          statstime[i][k] = stats[name]['avg']*0.001/60.0

    self.pmin, self.pmax, self.pvmin, self.pvmax = self.get_minmax(statstime)

    return statstime

  def get_sum_statstime(self, funclabels, funclist):
    self.funclabels = funclabels
    self.funclist = funclist

   #print('in get_sum_statstime')

    il = len(self.funclist)
    kl = len(self.nodelist)
    statstime = np.zeros((il, kl))
    for k in range(kl):
      stats = self.parstatslist[k]
      for i in range(il):
        name = self.funclist[i]
       #print('Node %d Func %d Name %s' %(k, i, name))
        for key in stats.keys():
          if(key.find(name) == 0):
            statstime[i][k] += stats[key]['avg']*0.001/60.0
       #print('statstime[%d][%d] = %f' %(i, k, statstime[i][k]))

    self.pmin, self.pmax, self.pvmin, self.pvmax = self.get_minmax(statstime)

    return statstime

  def get_sum_components(self, sumname):
    self.funclabels = ['sum']
    self.funclist = ['sum']

    print('Working on sumname: %s' %(sumname))

    ns = len(sumname)
    il = 0
    kl = len(self.nodelist)
    stats = self.parstatslist[0]
    for key in stats.keys():
      if(key.find(sumname) == 0):
        il += 1
        self.funclist.append(key)
        compname = key[ns+2:]
        print('sumname: %s No. %d name: %s' %(sumname, il, compname))
        self.funclabels.append(compname)

    il = len(self.funclist)
    kl = len(self.nodelist)
    statstime = np.zeros((il, kl))
    for k in range(kl):
      stats = self.parstatslist[k]
      keys = stats.keys()
      for i in range(1, il):
        name = self.funclist[i]
        if(name in keys):
          statstime[i][k] += stats[name]['avg']*0.001/60.0
          statstime[0][k] += stats[name]['avg']*0.001/60.0

    self.pmin, self.pmax, self.pvmin, self.pvmax = self.get_minmax(statstime)

    return statstime

  def plot(self, statstime, statsname):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

    fig = plt.figure()
    ax = plt.subplot()

    txtname = '%s/timing_%s.csv' %(self.casename, statsname)
    OPF = open(txtname, 'w')
    header = '%s Avg Time (Minutes)\n' %(statsname)
    OPF.write(header)

    print('text file: %s' %(txtname))

    header = '%-40s' %('Function Name')
    for i in range(nl):
      header = '%s, %8d' %(header, i)
    OPF.write(header+'\n')

    for i in range(len(self.funclist)):
      txtinfo = '%-40s' %(self.funclist[i])
      for k in range(nl):
        y[k] = statstime[i][k]
        txtinfo = '%s, %8.2f' %(txtinfo, y[k])
     #print('y = ', y)
      ax.plot(x, y, color=self.colorlist[i], linewidth=2, alpha=0.9)
      OPF.write(txtinfo+'\n')
    OPF.close()

    ylp = []
    ylabels = []
    pv = self.pvmin
    while(pv <= self.pvmax):
      ylp.append(pv)
      lbl = '%6.2f' %(pv)
      ylabels.append(lbl)
      pv *= 2.0

    if(self.linear):
      plt.xscale('linear')
    else:
      plt.xscale('log', base=2)
      plt.yscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xticks(x, xlabels)
     #plt.xticks(x, xlabels, rotation ='vertical')
      plt.yticks(ylp, ylabels)

    if(self.linear == 0):
      for i in range(len(self.funclist)):
        for k in range(nl):
          fact = 1.0/np.log2(2*self.nodelist[k])
          z[k] = statstime[i][0]*fact
       #https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
        ax.plot(x, z, color='black', linewidth=1, alpha=0.5, linestyle='dotted')

    plt.grid()

   #Same limits for everybody!
    plt.xlim(x[0], x[-1])
   #print('pmin: %f, pmax: %f' %(pmin, pmax))
   #plt.ylim(pmin, pmax)
    print('pmin: %f, pmax: %f' %(self.pmin, self.pmax))
    print('pvmin: %f, pvmax: %f' %(self.pvmin, self.pvmax))
    plt.ylim(self.pvmin, self.pvmax)
 
   #general title
   #title = '%s Timing (in Minutes), min: %8.2f, max: %8.2f' %(self.casename, pmin, pmax)
   #title = '%s Timing (in Minutes)' %(self.casename)
   #title = '%s Timing (in Minutes) of %s' %(statsname, self.casename)
    title = '%s Timing (in Minutes)' %(statsname)
    print('plot title: %s' %(title))
   #plt.suptitle(title, fontsize=13, fontweight=0, color='black', style='italic', y=1.02)
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
    plt.subplots_adjust(bottom=0.2, right=0.70, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (Minutes)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=self.funclabels,
           loc="center right",   # Position of legend
           fontsize=8,
           borderpad=1.2,
           labelspacing=1.2,
           handlelength=1.5
           )

   #Adjust the scaling factor to fit your legend text completely outside the plot
   #(smaller value results in more space being made for the legend)

    if(self.linear):
      imgname = '%s/lin_%s_timing.png' %(self.casename, statsname)
    else:
      imgname = '%s/log_%s_timing.png' %(self.casename, statsname)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  casename = 'dev'
  workdir = '/work2/noaa/gsienkf/weihuang/jedi/run.soca'
 #casename = 'anna'
 #workdir = '/work2/noaa/gsienkf/weihuang/ufs/soca/new-soca-solver'
  nodelist = [2, 4, 6, 8, 10, 12]
  corepernode = 36
  output = 0
  linear = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'workdir=', 'output=',
                             'corepernode=', 'nodelist=', 'casename='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--workdir'):
      workdir = a
    elif o in ('--corepernode'):
      corepernode = int(a)
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

  pr = Profiler(debug=debug, corepernode=corepernode, nodelist=nodelist, output=output,
                workdir=workdir, casename=casename, linear=linear)
  pr.process()
  pr.set_linear(linear=linear)
  pr.set_output(output=output)

  main_funclabels = ['total',
                     'computeHofX',
                     'ObsSpace',
                     'LocalInterpolator',
                     'GetValues',
                     'State']
  main_funclist = ['util::Timers::Total',
                   'oops::LocalEnsembleSolver::computeHofX',
                   'oops::ObsSpace::ObsSpace',
                   'oops::LocalInterpolator::apply',
                   'oops::GetValues::process',
                   'oops::State::write']
  statsname = '%s_main' %(casename)
  statstime = pr.get_main_statstime(main_funclabels, main_funclist)
  pr.plot(statstime, statsname)

  sum_funclabels = ['sum(oops::GetValues)',
                    'sum(oops::ObsError)',
                    'sum(oops::ObsFilter)',
                    'sum(oops::ObsSpace)',
                    'sum(oops::ObsVector)',
                    'sum(oops::ObsOperator)',
                    'sum(oops::VariableChange)']
  sum_funclist = ['oops::GetValues',
                  'oops::ObsError',
                  'oops::ObsFilter',
                  'oops::ObsSpace',
                  'oops::ObsVector',
                  'oops::ObsOperator',
                  'oops::VariableChange']
  statsname = '%s_sum' %(casename)
  statstime = pr.get_sum_statstime(sum_funclabels, sum_funclist)
  print('Ready to plot sum functions')
  pr.plot(statstime, statsname)

  for sumname in sum_funclist:
    item = sumname.split('::')
    name = item[1]
    statsname = '%s_%s' %(casename, name)
    statstime = pr.get_sum_components(sumname)
    pr.plot(statstime, statsname)

