# To split MS Calibrator
# generate list of .ms.split.cal
# check if it is all 12m array, is fo split using calibratorTools.py
# Run:
#   casa -c extractCal.py PROJECTNAME


import os
import glob
import calibratorTools as cT


def generateMSlist(location):
    # location = project directory
    pwd = os.getcwd()
    
    os.chdir(location)

    # find all MEMBER directory
    filesDepth = glob.glob('*/*/*/*')   # projname/science/group/member
    dirsDepth = filter(lambda f: os.path.isdir(f), filesDepth)
    
    AllMSlist = []
    MSlist_tobe_splitted = []
    No_calibrated = []
    No_mssplitcal = [] # single dish
    array7 = []

    print(dirsDepth)
    print("---------------------------")

    for dirs in dirsDepth:
        os.chdir(dirs) # enter Member dir
        if os.path.exists("calibrated"):
            os.chdir("calibrated")
            newloc = os.getcwd()
            print("Search *ms.split.cal in directory: ", os.getcwd())
            print("-------------------")
            MSlist = []
            for ifile in glob.glob("*.ms.split.cal"):
                print(ifile)
                MSlist.append(newloc+"/"+ifile)
                AllMSlist.append(newloc+"/"+ifile)

            print("-------------------")
            
            if len(MSlist) == 0:
                No_mssplitcal.append(newloc)
            else:
                for ifile in MSlist:
                    print("Checking MS: ", ifile)
                    # check if contain 7m array
                    arr7 = False
                    tb.open(ifile + "/ANTENNA")
                    dish_diameter = tb.getcol("DISH_DIAMETER") 
                    seven = dish_diameter[dish_diameter < 10.]
                    if len(seven) > 0: arr7 = True
                    tb.close()

                    if not arr7:
                        MSlist_tobe_splitted.append(ifile)
                    else:
                        array7.append(ifile)

        else:
            print("No /calibrated directory in this member: ", os.getcwd())
            No_calibrated.append(os.getcwd())

        os.chdir(location)

    # print "Calibrated dir not Found: ", No_calibrated
    # print "All MS :", AllMSlist, len(AllMSlist)
    # print "To be splitted :", MSlist_tobe_splitted, len(MSlist_tobe_splitted)
    # print "Contain 7m-array: ", array7
    # print "No *.ms.split.cal: ", No_mssplitcal

    os.chdir(pwd)

    return MSlist_tobe_splitted, AllMSlist, No_calibrated, No_mssplitcal, array7



def runextract(structFile, MSlist):
    cS = cT.calStructure(structFile, MSlist)
    cS.run()



def printlist(arr):
    for i in arr:
        print(i)


if __name__ == '__main__':
    dirhome = os.getcwd()
    dirname = dirhome + "/" + sys.argv[3]

    MSlist_tobe_splitted, AllMSlist, No_calibrated, No_mssplitcal, array7 = generateMSlist(dirname)
    runextract("structureFile.par", MSlist_tobe_splitted)

    # For logging
    print("Calibrated dir not Found: ", len(No_calibrated))
    printlist(No_calibrated)

    print("All MS: ", len(AllMSlist))
    printlist(AllMSlist)
    
    print("Splitted: ", len(MSlist_tobe_splitted))
    printlist(MSlist_tobe_splitted)

    print("Contain 7m-array: ", len(array7))
    printlist(array7)

    print("No *.ms.split.cal: ", len(No_mssplitcal))
    printlist(No_mssplitcal)