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

def cmdout(command):
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
    ostr = result.stdout
    return ostr.strip()

""" Profiler """
class Profiler:
  """ Constructor """
  def __init__(self, debug=0, linear=0, output=0, corelist=None, filelist=None):

    """ Initialize class attributes """
    self.debug = debug
    self.linear = linear
    self.output = output
    self.corelist = corelist
    self.filelist = filelist

    if(corelist is None):
      print('corelist not defined. Exit.')
      sys.exit(-1)

    if(filelist is None):
      print('filelist not defined. Exit.')
      sys.exit(-1)

    nc = 0
    for flnm in self.filelist:
      if(os.path.exists(flnm)):
        pass
      else:
        print('Filename ' + flnm + ' does not exit. Stop')
        sys.exit(-1)

    self.functionlist = ['Total', 'ObsSpace', 'write', 'measurementUpdate', 'Stats',
                         'computeHofX', 'packEigen', 'mask', 'ObsVector',
                         'packEigenSize', 'computeLocalization', 'ones',
			 'ObsVector_print', 'read', 'State_print',
                         'ObsLocalization', 'Geometry']
    self.colorlist = ['red', 'blue', 'cyan', 'magenta', 'orange',
                      'royalblue', 'coral', 'palegreen', 'deepskyblue', 'hotpink',
                      'dimgrey', 'brown', 'yellow', 'teal', 'violet', 'skyblue',
		      'tomato', 'lawngreen', 'steelblue', 'crimson']

  def process(self):
    self.stats_list = []
    nc = 0
    for flnm in self.filelist:
      nc += 1
      if(self.debug):
        print('Processing case ' + str(nc) + ': ' + flnm)
      res = self.stats(flnm)
      self.stats_list.append(res)

     #self.print_gen(nc)
     #self.print_par(nc)

  def print_gen(self, nc):
    print('\nTiming Statistics: Name')
    nl = len(self.stats_list[0][1])

    hinfo = '+----+'
    for i in range(nc):
      hinfo = '%s%s+' %(hinfo, 22*'-')
    print(hinfo)

    pinfo = '+-NP-+'
    for i in range(nc):
      pinfo = '%s%s%4d%s+' %(pinfo, 9*' ', self.corelist[i], 9*' ')
    print(pinfo)

    print(hinfo)

    n = 0
    for n in range(nl-10):
      ni = n + 1

      pinfo = '|%3d |' %(ni)
      tinfo = '|    |'
      for i in range(nc):
       #print('self.stats_list[i][1] = ', self.stats_list[i][1])
        idx = self.stats_list[i][1][n]
        name = self.stats_list[i][0][idx]['name']
       #pinfo = '%s %20s |' %(pinfo, self.stats_list[i][0][idx]['name'])
        nlst = name.split('::')
        pinfo = '%s %20s |' %(pinfo, nlst[-1])

        tim = self.stats_list[i][0][idx]['time']
        tinfo = '%s %s %10.2f %s |' %(tinfo, 4*' ', tim, 4*' ')
      print(pinfo)
      print(tinfo)
      print(hinfo)

  def print_par(self, nc):
    print('\nParallel Timing Statistics: Name')
    nl = len(self.stats_list[0][3])

    hinfo = '+----+'
    for i in range(nc):
      hinfo = '%s%s+' %(hinfo, 22*'-')
    print(hinfo)

    pinfo = '+-NP-+'
    for i in range(nc):
      pinfo = '%s%s%4d%s+' %(pinfo, 9*' ', self.corelist[i], 9*' ')
    print(pinfo)

    print(hinfo)

    n = 0
    for n in range(nl-10):
      ni = n + 1

      pinfo = '|%3d |' %(ni)
      tinfo = '|    |'
      for i in range(nc):
       #print('self.stats_list[i][1] = ', self.stats_list[i][1])
        idx = self.stats_list[i][3][n]
        name = self.stats_list[i][2][idx]['name']
       #pinfo = '%s %20s |' %(pinfo, self.stats_list[i][2][idx]['name'])
        nlst = name.split('::')
        pinfo = '%s %20s |' %(pinfo, nlst[-1])

        tim = self.stats_list[i][2][idx]['avg']
        tinfo = '%s %s %10.2f %s |' %(tinfo, 4*' ', tim, 4*' ')
      print(pinfo)
      print(tinfo)
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
      num_lines = len(lines)
      print('Total number of lines: ', num_lines)

      nl = 0
      while(nl < num_lines):
        if(lines[nl].find('Parallel Timing Statistics') > 0):
          if(self.debug):
            print('Start Parallel Timing Statistics')
          nl, par_stats, index = self.parallel_time_stats(lines, nl)
          prof[2] = par_stats
          prof[3] = index
          nl += num_lines
        elif(lines[nl].find('Timing Statistics') > 0):
          if(self.debug):
            print('Start Timing Statistics')
          nl, gen_stats, index = self.time_stats(lines, nl)
          prof[0] = gen_stats
          prof[1] = index
        nl += 1

    return prof

  def time_stats(self, lines, nl):
    stats = {}
    going = 1
    idx = 0
    ns = nl + 2
    while(going):
      line = lines[ns].strip()

      ns += 1
      
      if(line.find('Timing Statistics') > 0):
        going = 0
        break

      if(line.find('OOPS_STATS Name') >= 0):
        continue

     #print('Line ' + str(ns) + ': ' + line)

      item = line.split(' : ')
     #print(item)
      nlist = item[0].strip().split('::')

      name = nlist[-1]
      if(name == 'print'):
        if(nlist[-2] == 'ObsVector'):
          name = 'ObsVector_print'
        elif(nlist[-2] == 'State'):
          name = 'State_print'

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      tlist = tstr.split(' ')
     #print('tlist = ', tlist)
      ft = float(tlist[0])

      if(ft < 1.0):
        continue

      tinfo = {}
      tinfo['name'] = name
      tinfo['time'] = ft/1000.0
      tinfo['call'] = int(tlist[1])/1000.0

      stats[idx] = tinfo

     #pinfo = '\t%50s%10.2f%8d' %(stats[idx]['name'], stats[idx]['time'], stats[idx]['call'])
     #print(pinfo)
      idx += 1

    index = self.get_index(stats, 'time', 1)
   #if(self.debug):
   #  for idx in index:
   #    pinfo = '\t%50s%10.2f%8d' %(stats[idx]['name'], stats[idx]['time'], stats[idx]['call'])
   #    print(pinfo)

    return ns, stats, index

  def parallel_time_stats(self, lines, nl):
    stats = {}
    going = 1
    ns = nl + 3
    idx = 0
    while(going):
      line = lines[ns].strip()
      ns += 1
      if(line.find('Parallel Timing Statistics') > 0):
        going = 0
        break

     #print('Line ' + str(ns) + ': ' + line)

      item = line.split(' : ')
     #print(item)
      nlist = item[0].strip().split('::')

      name = nlist[-1]
      if(name == 'print'):
        if(nlist[-2] == 'ObsVector'):
          name = 'ObsVector_print'
        elif(nlist[-2] == 'State'):
          name = 'State_print'

      tstr = item[1].strip()
      while(tstr.find('  ') > 0):
        tstr = tstr.replace('  ', ' ')
      nlist = tstr.split(' ')
      ft = float(nlist[0])
     #if(ft < 1.0):
     #  continue

      tinfo = {}
      tinfo['name'] = name
      tinfo['min'] = ft/1000.0
      tinfo['max'] = float(nlist[1])/1000.0
      tinfo['avg'] = float(nlist[2])/1000.0

      tinfo['total'] = float(nlist[3])
   
     #imbalance = 100. * (max - min) / avg
      tinfo['imbalance'] = float(nlist[4])

      stats[idx] = tinfo

     #pinfo = '\t%50s%10.2f' %(stats[idx]['name'], stats[idx]['avg'])
     #print(pinfo)
     #print('\tName: ' + name + ' avg ' + str(stats[name]['avg']))
      idx += 1

    index = self.get_index(stats, 'avg', 2)
   #if(self.debug):
   #  for idx in index:
   #    pinfo = '\t%50s%10.2f' %(stats[idx]['name'], stats[idx]['avg'])
   #    print(pinfo)
    return ns, stats, index

  def get_index(self, stats, varname, nleft):
    s = []
    keys = stats.keys()
    for n in range(len(keys) - nleft):
      t = stats[n][varname]
      s.append(t)
    index = sorted(range(len(s)), key=lambda k: s[k])
    return index[::-1]

  def plot_total(self):
    title = 'Elapes, Measured, and Total time'
    x = np.array(self.corelist)

    y1 = []
    y2 = []

    nc = len(self.corelist)

    nl = len(self.stats_list[0][2])
    name1 = self.stats_list[0][2][nl-1]['name']
    name2 = self.stats_list[0][2][nl-2]['name']
   #print('self.stats_list[0][2][nl-1]=', self.stats_list[0][2][nl-1])
   #print('self.stats_list[0][2][nl-2]=', self.stats_list[0][2][nl-2])

    for i in range(nc):
      nl = len(self.stats_list[i][2])
      t1 = self.stats_list[i][2][nl-1]['avg']
      t2 = self.stats_list[i][2][nl-2]['avg']
      y1.append(t1)
      y2.append(t2)

    nlst1 = name1.split('::')
    nlst2 = name2.split('::')
    label1 = nlst1[-1]
    label2 = nlst2[-1]
    title = label1 + '-' + label2

   #multiple line plot
    lin1 = plt.plot(x, y1, marker='+', markerfacecolor='blue', markersize=12, color='skyblue', linewidth=4, label=label1)
    lin2 = plt.plot(x, y2, marker='x', markerfacecolor='magenta', markersize=12, color='red', linewidth=4, label=label2)

   #Add title and axis names
    plt.title(title)
    plt.xlabel('Procs')
    plt.ylabel('Time(sec)')

    if(self.linear):
      pass
    else:
      plt.xscale("log", base=2)
      plt.yscale("log", base=2)
     #plt.yscale("log", base=10)

    plt.legend()
   #plt.legend(handles=[lin1, lin2])

    if(self.linear):
      imgname = 'lin_totalNmeasured.png'
    else:
      imgname = 'log_totalNmeasured.png'

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

    plt.cla()
    plt.clf()

  def find_gen_name_idx(self, cn, idx):
    name = self.stats_list[0][0][idx]['name']
    if(cn):
      nl = len(self.stats_list[cn][0])
      indx = -1
      for n in range(nl):
        if(name == self.stats_list[cn][0][n]['name']):
         #print('name = ' + name + ', new name = ' + self.stats_list[cn][0][n]['name'])
          indx = n
          break
    else:
      indx = idx
    return indx

  def plot_multiple(self, ns, nl):
    title = 'Function %d-%d' %(ns, ns+nl)

    nc = len(self.corelist)
    x = np.array(self.corelist)
    ts = {}
    namelist = []
    line = []

    ymin = 1.0e+36
    ymax = 0.0

    for n in range(nl):
      nt = ns + n 
      idx = self.stats_list[0][1][nt]
      name = self.stats_list[0][0][idx]['name']
      nlst = name.split('::')
      namelist.append(nlst[-1])
      pl = None
      line.append(pl)
      ts[n] = []
      for i in range(nc):
        ni = self.find_gen_name_idx(i, idx)
        t = self.stats_list[i][0][ni]['time']
       #print('n = ' + str(n) + ', i = ' + str(i) + ', ni = ' + str(ni) + ', t = ' + str(t))
        ts[n].append(t)
        if(t > ymax):
          ymax = t
        if(t < ymin):
          ymin = t

   #t0 = ymax
   #print('t0=', t0)
   #print('x=', x)
   #yt = x
   #for i in range(nc):
   #  t = t0*x[0]/x[i]
   #  yt[i] = t
 
   #multiple line plot
    for num in range(nl):
      line[n] = plt.plot(x, ts[num], marker='', color=self.colorlist[num],
                         linewidth=2, alpha=0.9, label=namelist[num])

   #Add title and axis names
   #plt.title(title)
    plt.xlabel('Procs')
    plt.ylabel('Time(ms)')

    if(self.linear):
      pass
    else:
      plt.xscale("log", base=2)
      plt.yscale("log", base=2)
     #plt.yscale("log", base=10)

    plt.legend()
 
    plt.title(title)

    if(self.linear):
      imgname = 'lin_plot_%d-%d.png' %(ns, ns+nl)
    else:
      imgname = 'log_plot_%d-%d.png' %(ns, ns+nl)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

    plt.cla()
    plt.clf()

  def find_par_name_idx(self, cn, idx):
    name = self.stats_list[0][2][idx]['name']
    if(cn):
      nl = len(self.stats_list[cn][2])
      indx = -1
      for n in range(nl):
        if(name == self.stats_list[cn][2][n]['name']):
         #print('name = ' + name + ', new name = ' + self.stats_list[cn][0][n]['name'])
          indx = n
          break
    else:
      indx = idx
    return indx

  def plot_par_bar(self, pv):
    nv = len(self.stats_list[0][3])
    if(nv < pv):
      print('Cannot plot ' + str(pv) + ' vars, as we only have ' + str(nv) + ' in file.')
      return

    title = '%d Parallel Functions' %(pv)

    nc = len(self.corelist)
    x = np.array(self.corelist)
    ts = {}
    idxlist = []
    namelist = []

    for i in range(nc):
      nv = len(self.stats_list[i][3])
      print('Case ' + str(i) + ' has ' + str(nv) + ' vars')

    for n in range(pv):
      idx = self.stats_list[0][3][n]
      idxlist.append(idx)
      name = self.stats_list[0][2][idx]['name']
      nlst = name.split('::')
      namelist.append(nlst[-1])

    ymin = 1.0e+36
    ymax = 0.0

    for i in range(nc):
      ts[i] = []
      for n in range(pv):
        idx = idxlist[n]
        ni = self.find_par_name_idx(i, idx)
        t = self.stats_list[i][2][ni]['avg']
        ts[i].append(t)
        if(t > ymax):
          ymax = t
        if(t < ymin):
          ymin = t

   #the actual graph:
    fig, ax = plt.subplots(figsize = (10,4))

    idx = np.asarray([i for i in range(pv)])

    width = (1.0-0.1)/len(corelist);

    lengenlist = []

    for n in range(len(corelist)):
      ax.bar(idx+n*width, ts[n], width, color=self.colorlist[n])
      lengenlist.append(str(self.corelist[n]))

    ax.set_xticks(idx)
    ax.set_xticklabels(namelist, rotation=65)
    ax.legend(lengenlist)
    ax.set_xlabel('Functions')
    ax.set_ylabel('Time (second)')

    plt.title(title)

    fig.tight_layout()

    imgname = 'top_bar_%d.png' %(pv)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

    plt.cla()
    plt.clf()

  def find_func_name_idx(self, cn, name):
    nl = len(self.stats_list[cn][2])
    indx = 0
    for n in range(nl):
      if(name == self.stats_list[cn][2][n]['name']):
       #print('name = ' + name + ', new name = ' + self.stats_list[cn][0][n]['name'])
        indx = n
        break
    return indx

  def plot_select_funcs_bar(self, varname='avg'):
    pv = len(self.functionlist)
    nv = len(self.stats_list[0][3])
    if(nv < pv):
      print('Cannot plot ' + str(pv) + ' vars, as we only have ' + str(nv) + ' in file.')
      return

    title = '%s time Bar Plot for selected Functions' %(varname)

    nc = len(self.corelist)
    x = np.array(self.corelist)
    ts = {}
    idxlist = []
    namelist = []

   #for i in range(nc):
   #  nv = len(self.stats_list[i][3])
   #  print('Case ' + str(i) + ' has ' + str(nv) + ' vars')

    ymin = 1.0e+36
    ymax = 0.0

    for i in range(nc):
      ts[i] = []
      for n in range(pv):
        name = self.functionlist[n]
        ni = self.find_func_name_idx(i, name)
        t = self.stats_list[i][2][ni][varname]
        ts[i].append(t)
        if(t > ymax):
          ymax = t
        if(t < ymin):
          ymin = t

   #the actual graph:
    fig, ax = plt.subplots(figsize = (10,4))

    idx = np.asarray([i for i in range(pv)])

    width = (1.0-0.1)/len(corelist);

    lengenlist = []

    print('self.colorlist: ', self.colorlist)

    for n in range(len(corelist)):
      print('No. %d color: %s' %(n, self.colorlist[n]))
      ax.bar(idx+n*width, ts[n], width, color=self.colorlist[n])
      lengenlist.append(str(self.corelist[n]))

    ax.set_xticks(idx)
    ax.set_xticklabels(self.functionlist, rotation=65)
    ax.legend(lengenlist)
    ax.set_xlabel('Functions')
    ax.set_ylabel('Time (second)')

    plt.title(title)

    fig.tight_layout()

    imgname = 'select_functions_%s_CPUtime_bar_plot.png' %(varname)

    if(self.output):
      plt.savefig(imgname)
    else:
      plt.show()

    plt.cla()
    plt.clf()

#--------------------------------------------------------------------------------
if __name__== '__main__':
  debug = 1
  linear = 0
  output = 0
  workdir = '.'
  mtpn = 6
  nodelist = [28, 30, 32, 35, 36, 40, 45, 50, 60]
  corelist = []
  for n in range(len(nodelist)):
    procs = mtpn*nodelist[n]
    corelist.append(procs)

  print("corelist:", corelist)
  print("nodelist:", nodelist)

 #rundir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/surf'
  rundir = '/scratch2/BMC/gsienkf/Wei.Huang/ufs/soca-solver'
  filelist = []
  for n in range(len(corelist)):
    flnm = '%s/soca_solver.%dt%dn_%dp/log.soca_solver' %(rundir, mtpn, nodelist[n], corelist[n])
    filelist.append(flnm)

  print('filelist:', filelist)

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'linear=', 'output='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--linear'):
      linear = int(a)
    elif o in ('--output'):
      output = int(a)
    else:
      assert False, 'unhandled option'

  pr = Profiler(debug=debug, linear=linear, output=output, corelist=corelist, filelist=filelist)
  pr.process()
 #pr.plot_total()
 #pr.plot_par_bar(10)
  pr.plot_select_funcs_bar(varname='min')
  pr.plot_select_funcs_bar(varname='max')
  pr.plot_select_funcs_bar(varname='avg')

