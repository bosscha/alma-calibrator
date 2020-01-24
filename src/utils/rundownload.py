# list all downloadscript file 
# FORM: download*_PROJNAME.sh
# make a directory based on the name (project name)
# mv the script to that dir
# chmod +x
# run the download script


# RUN:
#   python rundownload.py


import os
import sys
import glob
import time

def rundown(dfile):
    homepath = os.getcwd()
    print("Download file: ", dfile)
    projectname = dfile.split("_")[1][:-3] # get the project name, remove ".sh"
    print("Project name: ", projectname)

    if not os.path.exists(projectname):
        print("Create new dir: ", projectname)
        os.makedirs(projectname)

        print("Move the downloadfile to new dir.")
        os.rename(homepath+"/"+dfile, homepath+"/"+projectname+"/"+dfile)

        print("Make executable")
        os.chdir(homepath+"/"+projectname)
        os.system("chmod +x "+dfile)

        print("Run the download file")
        os.system("./"+dfile)

        print("Back to previous path.")
        os.chdir(homepath)


    else:
        print("Directory exist: ", projectname)


if __name__ == '__main__':
    delaytime = 0
    if len(sys.argv) > 1:
        delaytime = float(sys.argv[1]) * 3600.0

    time.sleep(delaytime)

    for dfile in glob.glob("downloadRequest*.sh"):
        start_time = time.time()
        rundown(dfile)
        print("Download time: {0} minutes".format(str((time.time() - start_time)/60.0)))

    print("Finish.")