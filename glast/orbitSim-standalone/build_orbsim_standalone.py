#!/usr/bin/env python
#08/27/2008 vernaleo Do something similar to what we do for standalone
#                    gtbin (using some routines from that script).

import os
import shutil
import time
import sys
sys.path.append('../evtbin-standalone')
import build_standalone

SLAC_version="v9r7p1"
GSSC_version="orb-g1_3"
#For an actual release set as true to get rid of the date added to the names
#and to use something closer to HEASARC naming for platforms.
release=0

def copy_files(outname):
    #Put everything in the proper directory structure for the standalone.
    standalonepath="../glast/orbitSim-standalone/"
    fullpath="../glast/orbitSim/"
    shutil.copy2(standalonepath+"README",standalonepath+outname)
    shutil.copy2(standalonepath+"orbitSim-init.csh",standalonepath+outname)
    shutil.copy2(standalonepath+"orbitSim-init.sh",standalonepath+outname)
    dir=standalonepath+outname+"/bin"
    if not os.path.isdir(dir):
        os.mkdir(dir)
    #for cygwin
    if os.path.isfile(fullpath+"src/orbSim/gtorbsim.exe"):
        shutil.copy2(fullpath+"src/orbSim/gtorbsim.exe",dir)
    else:
        shutil.copy2(fullpath+"src/orbSim/gtorbsim",dir)
    dir=standalonepath+outname+"/data"
    if not os.path.isdir(dir):
        os.mkdir(dir)
    shutil.copy2(fullpath+"data/ft2.tpl",dir)
    dir=standalonepath+outname+"/help"
    if not os.path.isdir(dir):
        os.mkdir(dir)
    shutil.copy2(fullpath+"src/orbSim/gtorbsim.txt",dir)
    dir=standalonepath+outname+"/pfiles"
    if not os.path.isdir(dir):
        os.mkdir(dir)
    shutil.copy2(fullpath+"pfiles/gtorbsim.par",dir)

###########MAIN################

if __name__ == "__main__":
    start=time.time()

    (outname,tarname)=build_standalone.prep_output("orbitSim",GSSC_version,release)
    startpath=os.getcwd()
    build_standalone.overlay_src()
    os.chdir("../../BUILD_DIR")

    cmd="./configure --enable-static --disable-shared --enable-readline"
    os.system(cmd)

    if not (os.uname()[0] == 'Darwin'):
        build_standalone.fix_hmakerc()

    os.system("./hmake")
    
    copy_files(outname)

    if not (os.uname()[0] == 'Darwin'):
        build_standalone.unfix_hmakerc()

    #I considered doing this, but ultimately probably not worth it.
    #os.system("./hmake clean")
    os.chdir(startpath)
    build_standalone.tar_output(outname,tarname)
    build_standalone.remove_overlay()

    end=time.time()
    print (end-start)/60.0,"minutes used to build and package."
