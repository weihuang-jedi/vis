#!/usr/bin/env python3
import argparse
import os
import sys
import glob
import time
from pathlib import Path

IODA_CONV_PATH = Path(__file__).parent/"../lib/pyiodaconv"
if not IODA_CONV_PATH.is_dir():
    IODA_CONV_PATH = Path(__file__).parent/'..'/'lib-python'
sys.path.append(str(IODA_CONV_PATH.resolve()))

import gsi_ncdiag as gsid


def run_conv_obs(convfile, outdir, platforms):
    print("Processing:"+str(convfile))
    startt = time.time()
    Diag = gsid.Conv(convfile)
    Diag.read()
    Diag.toIODAobs(outdir, platforms=platforms)
    Diag.close()
    print("Time (OBS) %s[%s]: %.3g sec" % (convfile, ",".join(platforms), time.time() - startt))
    return 0


def run_conv_geo(convfile, outdir):
    print("Processing:"+str(convfile))
    startt = time.time()
    Diag = gsid.Conv(convfile)
    Diag.read()
    Diag.toGeovals(outdir)
    print("Time (GEO) %s: %.3g sec" % (convfile, time.time() - startt))
    return 0


ScriptName = os.path.basename(sys.argv[0])

# Parse command line
ap = argparse.ArgumentParser()
ap.add_argument("input_dir", help="Path to concatenated GSI diag files")
ap.add_argument("-o", "--obs_dir",
                help="Path to directory to output observations")
ap.add_argument("-g", "--geovals_dir",
                help="Path to directory to output observations")
ap.add_argument("-d", "--obsdiag_dir",
                help="Path to directory to output observations")
ap.add_argument("-b", "--add_obsbias", default=False,
                help="Add ObsBias group to output observations")
ap.add_argument("-q", "--add_qcvars", default=False,
                help="Add QC variables to output observations")
ap.add_argument("-r", "--add_testrefs", default=False,
                help="Add TestReference group to output observations")

MyArgs = ap.parse_args()

DiagDir = MyArgs.input_dir

# process obs files
if MyArgs.obs_dir:
    ObsDir = MyArgs.obs_dir
    if not Path(ObsDir).is_dir():
        raise Exception("Obs dir: '%s' does not exist." % ObsDir)
    ObsBias = MyArgs.add_obsbias
    QCVars = MyArgs.add_qcvars
    TestRefs = MyArgs.add_testrefs
    # conventional obs first
    # get list of conv diag files
    convfiles = glob.glob(DiagDir+'/*conv*')
    for convfile in convfiles:
        splitfname = convfile.split('/')[-1].split('_')
        if 'conv' in splitfname:
            i = splitfname.index('conv')
            c = "_".join(splitfname[i:i + 2])
        try:
            for p in gsid.conv_platforms[c]:
                run_conv_obs(convfile, ObsDir, [p])
        except (KeyError, IndexError):
            pass
    # radiances next

# process geovals files
if MyArgs.geovals_dir:
    GeoDir = MyArgs.geovals_dir
    # conventional obs first
    # get list of conv diag files
    convfiles = glob.glob(DiagDir+'/*conv*')
    for convfile in convfiles:
        try:
            run_conv_geo(convfile, GeoDir)
        except (KeyError, IndexError):
            pass

