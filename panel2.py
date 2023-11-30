import os, sys
import matplotlib.pyplot as plt
from matplotlib import image
import numpy as np
import getopt

import matplotlib as mpl
mpl.rcParams['figure.dpi']= 300

#=========================================================================
class MakePanelPlot():
  def __init__(self, debug=0, output=0, img1='sfcship_ps',
               img2='surface_pressure'):
    self.debug = debug
    self.output = output
    self.img1 = img1
    self.img2 = img2

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.img1 = ', self.img1)
      print('self.img2 = ', self.img2)

    self.genit()

  def genit(self):
   #initialise grid of axes
   #fig, axes = plt.subplots(2, 2, constrained_layout=True)
    fig, axes = plt.subplots(1, 2)
    axes = axes.ravel()

   #create fake data
    img = [
        'FAKENAME_30.png',
        'FAKENAME_40.png',
    ]

   #iterate over axes
    for i, ax in enumerate(axes):
     #imageprefix = '%s_%s_increment_lev' %(img1, img2)
      imageprefix = '%s_lev' %(img1)
      imgname = img[i].replace('FAKENAME', imageprefix)
      im = image.imread(imgname)
      ax.axis('off')
      ax.imshow(im)

    plt.tight_layout()

    if(self.output):
      imgname = 'panelplot'
     #fig.savefig(imgname, dpi=300, bbox_inches="tight")
      fig.savefig(imgname, bbox_inches="tight", quality=100)
    else:
      plt.show()

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  img1 = 'sfcship_ps'
  img2 = 'surface_pressure'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                             'img1=', 'img2='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--img1'):
      img1 = a
    elif o in ('--img2'):
      img2 = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

  mpp = MakePanelPlot(debug=debug, output=output,
                      img1=img1, img2=img2)

