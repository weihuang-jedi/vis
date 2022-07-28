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
  def __init__(self, debug=0, corelist=[], casename=None, output=0,
               nodelist=[], workdir=None, linear=0):

    """ Initialize class attributes """
    self.debug = debug
    self.workdir = workdir
    self.corelist = corelist
    self.nodelist = nodelist
    self.casename = casename
    self.output = output
    self.linear = linear

    self.shortest_time = 0.1

    if(workdir is None):
      print('workdir not defined. Exit.')
      sys.exit(-1)

    if(len(corelist) < 1):
      print('corelist not defined. Exit.')
      sys.exit(-1)

    if(len(nodelist) < 1):
      print('nodelist not defined. Exit.')
      sys.exit(-1)

    self.colorlist = ['red', 'green', 'cyan', 'blue', 'magenta',
                      'firebrick', 'lime']
    self.linestyle = ['solid', 'solid', 'solid', 'solid', 'solid',
                      'dashed', 'dashed']
    self.function_list = ['sum(oops::GetValues)',
                          'sum(oops::ObsError)',
                          'sum(oops::ObsFilter)',
                          'sum(oops::ObsSpace)',
                          'sum(oops::ObsVector)',
                          'sum(oops::ObsOperator)',
                          'sum(oops::VariableChange)']
    self.fullfunction_list = ['oops::GetValues',
                              'oops::ObsError',
                              'oops::ObsFilter',
                              'oops::ObsSpace',
                              'oops::ObsVector',
                              'oops::ObsOperator',
                              'oops::VariableChange']

    self.statstime = []
    for i in range(len(self.fullfunction_list)):
      tmax = {}
      for n in range(len(self.nodelist)):
        tmax[n] = 0.0
      self.statstime.append(tmax)
      
  def set_linear(self, linear=1):
    self.linear = linear

  def set_output(self, output=1):
    self.output = output

  def process(self):
    for n in range(len(self.nodelist)):
      rundir = '%s/%s/run_80.40t%dn_%dp' %(self.workdir, self.casename,
                self.nodelist[n], self.corelist[n])
      flnm = '%s/stdoutNerr/stdout.00000000' %(rundir)
      if(self.debug):
        print('Processing node: %d, as file: %s' %(self.nodelist[n], flnm))
      if(os.path.exists(flnm)):
        pass
      else:
        print('Filename ' + flnm + ' does not exit. Stop')
        sys.exit(-1)

      with open(flnm) as fp:
        lines = fp.readlines()
       #line = fp.readline()
        num_lines = len(lines)
       #print('Total number of lines: ', num_lines)

        nl = 0
        while(nl < num_lines):
          if(lines[nl].find('Parallel Timing Statistics') > 0):
           #if(self.debug):
           #  print('lines[%d]: %s' %(nl, lines[nl]))
            self.parallel_time_stats(n, lines, nl)
            break
          nl += 1
  
  def parallel_time_stats(self, n, lines, nl):
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
      tmin = float(nlist[0])
      tmax = float(nlist[1])
      tavg = float(nlist[2])
      
      self.updatestats(n, name, tmax)

  def updatestats(self, n, name, tmax):
    for i in range(len(self.fullfunction_list)):
      if(name.find(self.fullfunction_list[i]) >= 0):
        self.statstime[i][n] = tmax
        return

  def plot(self):
    try:
      plt.close('all')
      plt.clf()
      plt.cla()
    except Exception:
      pass

    title = '%s Timing' %(self.casename)

    nl = len(self.nodelist)
    x = np.zeros((nl), dtype=float)
    y = np.zeros((nl), dtype=float)
    z = np.zeros((nl), dtype=float)
    xlabels = []
    for k in range(nl):
      x[k] = self.nodelist[k]
      lbl = '%d' %(self.nodelist[k])
      xlabels.append(lbl)

    pmin = 1.0/128.0
    pmax = 256.0
    ylabels = []
    yp = []
    pcur = pmin/2.0
    while(pcur < pmax):
      pcur *= 4.0
      lbl = '%6.2f' %(pcur)
      ylabels.append(lbl)
      yp.append(pcur)
      print('yp = ', yp)

    fig = plt.figure()
    ax = plt.subplot()

    if(self.linear):
      plt.xscale('linear')
    else:
     #plt.xscale('log', base=2)
     #plt.yscale('log', base=10)
      plt.xscale('log', basex=2)
      plt.yscale('log', basey=2)
      plt.xticks(x, xlabels)
     #plt.xticks(x, xlabels, rotation ='vertical')
      plt.yticks(yp, ylabels)

    pmin = 1.0e20
    pmax = 0.0

    txtname = 'obs_timing_%s.csv' %(self.casename)
    OPF = open(txtname, 'w')
    header = '%40s' %('Function Name')
    for k in range(nl):
      header = '%s, %12d' %(header, self.nodelist[k])
    OPF.write(header+'\n')

    for i in range(len(self.function_list)):
      print('self.function_list[%d] = %s' %(i, self.function_list[i]))
      print('self.statstime[i] = ', self.statstime[i])
      txtinfo = '%40s' %(self.function_list[i])
      npnts = 0
      for k in range(nl):
        y[k] = 0.001*self.statstime[i][k]/60.0
        if(pmin > y[k]):
          pmin = y[k]
        if(pmax < y[k]):
          pmax = y[k]
        txtinfo = '%s, %12.2f' %(txtinfo, y[k])
        if(y[k] > self.shortest_time):
          npnts += 1
      OPF.write(txtinfo+'\n')
      print('x = ', x[0:npnts])
      print('y = ', y[0:npnts])
      print('self.colorlist[%d] = %s' %(i, self.colorlist[i]))
      print('self.linestyle[%d] = %s' %(i, self.linestyle[i]))
      if(npnts > 1):
        ax.plot(x[0:npnts], y[0:npnts], color=self.colorlist[i], linewidth=2,
                linestyle=self.linestyle[i], alpha=0.9)
    OPF.close()

    if(self.linear == 0):
      for i in range(len(self.fullfunction_list)):
        npnts = 0
        for k in range(nl):
          fact = 1.0/np.log2(2*self.nodelist[k])
          z[k] = 0.001*self.statstime[i][0]*fact/60.0
          if(z[k] > self.shortest_time):
            npnts += 1
       #https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
        if(npnts > 1):
          ax.plot(x[0:npnts], z[0:npnts], color='black', linewidth=1, alpha=0.5, linestyle='dotted')

    plt.grid()

   #Same limits for everybody!
    print('pmin: %f, pmax: %f' %(pmin, pmax))

   #if(pmin < self.shortest_time):
   #  pmin = self.shortest_time
   #if(pmax < 10000.0):
   #  pmax = 10000.0

    pmin = 1.0/128.0
    pmax = 256.0
   
    plt.xlim(x[0], x[-1])
    plt.ylim(pmin, pmax)
 
   #general title
   #title = '%s Timing (in minutes), min: %8.2f, max: %8.2f' %(self.casename, pmin, pmax)
    title = '%s Timing (in minutes)' %(self.casename)
   #plt.suptitle(title, fontsize=13, fontweight=0, color='black', style='italic', y=1.02)
    plt.suptitle(title, fontsize=16, fontweight=1, color='black')

   #Create a big subplot
    bs = fig.add_subplot(111, frameon=False)
   #plt.subplots_adjust(bottom=0.2, right=0.70, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.675, top=0.8)
    plt.subplots_adjust(bottom=0.2, right=0.65, top=0.8)
   #plt.subplots_adjust(bottom=0.2, right=0.5, top=0.8)

   #hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    bs.set_xlabel('Node', labelpad=10) # Use argument `labelpad` to move label downwards.
    bs.set_ylabel('Time (minutes)', labelpad=20)

   #Create the legend
    fig.legend(ax, labels=self.function_list,
           loc="center right",   # Position of legend
           fontsize=6,
           borderpad=1.2,
           handlelength=1.5
           )

#          borderpad=1.2,
#          labelspacing=1.2,
#          handlelength=1.5

   #Adjust the scaling factor to fit your legend text completely outside the plot
   #(smaller value results in more space being made for the legend)

    if(self.linear):
      imgname = 'lin_%s_obs_timing.png' %(self.casename)
    else:
      imgname = 'log_%s_obs_timing.png' %(self.casename)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  casename = 'sondes'
 #workdir = '/work2/noaa/gsienkf/weihuang/jedi/case_study'
  workdir = '/scratch2/BMC/gsienkf/Wei.Huang/tools/vis_tools'
  corelist = [36, 78, 156, 312]
 #corelist = [36, 72, 144, 288]
  nodelist = [1, 2, 4, 8]
  output = 0
  linear = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'workdir=',
                             'corelist=', 'nodelist=', 'casename='])

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

  pr = Profiler(debug=debug, corelist=corelist, nodelist=nodelist, output=output,
                workdir=workdir, casename=casename, linear=linear)
  pr.process()
  for linear in [1, 0]:
    pr.set_linear(linear=linear)
   #for output in [0, 1]:
    for output in [1]:
      pr.set_output(output=output)
      pr.plot()

