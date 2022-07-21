#!/bin/bash

 export PROJ_LIB=/work/noaa/gsienkf/weihuang/anaconda3/share/proj
 export PYTHONPATH=/work/noaa/gsienkf/weihuang/jedi/vis_tools/xESMF/build/lib:/work/noaa/gsienkf/weihuang/anaconda3/lib
 export LD_LIBRARY_PATH=/work/noaa/gsienkf/weihuang/anaconda3/lib:${LD_LIBRARY_PATH}
 export PATH=/work/noaa/gsienkf/weihuang/anaconda3/bin:${PATH}

 python interp2latlon.py 

