# Last update: 2017.07.10
# Author: Ridlo W. Wibowo
# Not GENARAL
# Be careful
#
# Script to run scriptForPI.py without getting some errors
# list of error:
#       - error of CASA VERSION
#       - error of tgz
#
# Run: 
#   python getms.py

import os
import sys
import glob


def targz_to_tgz():
    """scriptForPI.py crash when restoring older calibrations in CASA versions 4.7.0 and above"""
    # Need to change the extensions of some file, from .tar.gz to .tgz
    print("List all *.tar.gz files in this directory: "+os.getcwd())
    listtargz = []
    for ifile in glob.glob("*.tar.gz"):
        print(ifile)
        listtargz.append(ifile)

    for ifile in listtargz:
        print "Moving " + ifile + " to " + ifile[:-7]+".tgz"
        os.system("mv "+ifile+" "+ifile[:-7]+".tgz")


def commentout_version():
    """comment out 'CASA version' checking"""

    print("List all *.py files in this directory: "+os.getcwd())
    listscript = []
    for ifile in glob.glob("*.py"):
        print(ifile)
        listscript.append(ifile)


    for ifile in listscript:
        if ifile != 'getms.py': # except this file
            print "\nSearching on: " + ifile
            
            f = open(ifile, "r")
            searchlines = f.readlines()
            f.close()
            
            for i, line in enumerate(searchlines):
                # comment out 2 lines 
                if ("if re.search(" in line) or ("if casadef.casa_version" in line):
                    searchlines[i] = '# ' + searchlines[i]
                    searchlines[i+1] = '# ' + searchlines[i+1]
                    print searchlines[i], searchlines[i+1]
                    print "Comment out."

            of = open(ifile, "w")
            for line in searchlines:
                of.write(line)
            of.close()

            print "Writing down again"


def run():
    # GETMS in "script" directory
    if os.path.exists("../calibrated"):
        print "Remove previous ./calibrated dir.."
        os.system("rm -rf ../calibrated")

    pipeline = False
    for ifile in glob.glob("*.xml"):
        if "PPR" in ifile:
            pipeline = True
            print "Found PPR file. Use pipeline.."

    if pipeline:
        os.chdir("../calibration")
        targz_to_tgz()
        os.chdir("../script")
        commentout_version()
        os.system('casa --pipeline -c scriptForPI.py')
    else:
        commentout_version()
        os.system('casa -c scriptForPI.py')


def run_project(location="/home/projectname"):
    # list all directory until certain depth
    # run "run" from script dir
    # get all MS for a project (all member inside)
    os.chdir(location)

    print "List all member, in directory: ", os.getcwd()
    filesDepth = glob.glob('*/*/*/*')   # projname/science/group/member
    dirsDepth = filter(lambda f: os.path.isdir(f), filesDepth)
    print dirsDepth

    for dirs in dirsDepth:
        os.chdir(dirs) # enter
        if os.path.exists("script"):
            os.chdir("script")
            print "Run 'getms' in the directory: ", os.getcwd()
            run() # run
        else:
            print "No script directory in ", os.getcwd()

        os.chdir(location) # back



def untar(dirname):
    os.chdir(dirname)

    listfile = glob.glob("*tar") 
    print listfile

    for fl in listfile:
        cmd = "tar xvf %s"%(fl)
        os.system(cmd) 



if __name__ == '__main__':
    # run()
    dirname = "/mnt/sciops/data/rwibowo/projects/2013.1.00532.S"
    untar(dirname)
    run_project(dirname)


