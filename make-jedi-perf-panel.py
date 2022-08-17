import os, sys
import matplotlib.pyplot as plt
from matplotlib import image
import numpy as np
import getopt

import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

#=========================================================================
class MakePanelPlot():
  def __init__(self, debug=0, output=0, caselist=[], funclist=[]):
    self.debug = debug
    self.output = output
    self.caselist = caselist
    self.funclist = funclist

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.caselist = ', self.caselist)
      print('self.funclist = ', self.funclist)

    for funcname in self.funclist:
      self.genit(funcname)

  def genit(self, funcname):
   #initialise grid of axes
   #fig, axes = plt.subplots(2, 2, constrained_layout=True)
    fig, axes = plt.subplots(2, 2)
    axes = axes.ravel()

   #iterate over axes
    for i, ax in enumerate(axes):
      imgname = '%s/log_%s_%s_timing.png' %(self.caselist[i],
                self.caselist[i], funcname)
      im = image.imread(imgname)
      ax.axis('off')
      ax.imshow(im)

   #plt.subplots_adjust(left=0.05, bottom=0.05,
   #                    right=0.95, top=0.95,
   #                    wspace=0.05, hspace=0.05)
   #The parameter meanings (and suggested defaults) are:
   #left  = 0.125  # the left side of the subplots of the figure
   #right = 0.9    # the right side of the subplots of the figure
   #bottom = 0.1   # the bottom of the subplots of the figure
   #top = 0.9      # the top of the subplots of the figure
   #wspace = 0.2   # the amount of width reserved for blank space between subplots
   #hspace = 0.2   # the amount of height reserved for white space between subplots

    plt.tight_layout()

    if(self.output):
      imgname = 'panel-%s' %(funcname)
     #fig.savefig(imgname, dpi=300, bbox_inches="tight")
      fig.savefig(imgname, bbox_inches="tight", quality=100)
    else:
      plt.show()

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  caselist = ['aircraft', 'iasi', 'scatwind', 'sondes']
  funclist = ['main', 'sum', 'GetValues', 'ObsError',
              'ObsSpace', 'ObsVector',
              'ObsOperator', 'VariableChange', 'ObsFilter']

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                             'caselist=', 'funclist='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--caselist'):
      caselist = a
    elif o in ('--funclist'):
      funclist = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

  mpp = MakePanelPlot(debug=debug, output=output,
                      caselist=caselist, funclist=funclist)

