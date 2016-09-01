#!/usr/bin/python

"""

Script to manipulate the MS for the calibrators.


HISTORY:
    2015.09.23:
        -first shot to extract from a calibrated MS the split files of the calibrators
        
    2015.09.24:
        - fixing the splitting
     
    2015.09.25:
        - fixing the splitting
        
    2015.10.01:
        - updating script
        
        
    2015.11.17:
        - keep up the calibrator information in a DB and re-organize the MS of the split-file using the class calStructure.
        - return the list of MS and calibrator in calibrator
        
        
    2015.12.28:
        - adding a class imagingCalibrator to provide some imaging function (continuum or line)
        
        
    2015.12.29:
        - go through the MS list and move the split MS to the corresponding directory
        - adding the band structure
        - !! bug, it seems that the size double if a second run....!!!!!!!!!
        
    2016.01.07:
        - fix of  the previous bug to skip the mv if a ms file exists already.
        - store the information of the dataset in the DB
        
    2016.08.01:
        - add an exception with a problem from splitting
        
    2016.08.11:
        - bug fixes
        
    2016.08.12:
        - add check for split
        
    2016.08.16:
        - add the reference position of the field
        - add an id for the calibrator to identify them uniquely (DB)
        - add reference position refer, m0, m1 for the calibrator (DB)
        
    2016.08.18:
        - add TOLDISTCAL tolerance distance keyword for calibrator ID
        - check if calibrator id exists in DB otherwise create a new one.
        
    
    2016.08.30:
        - fixes for the position check for the calibrator
        
    2016.09.01:
        - check if the calibrator exists to move in the same directory but with the MS name.
        - add a check for solar system object which have no fix coordinates but assuming a current name
        

RUN:

$> casa -c scriptToOrganize.py


"""


__author__="S. Leon @ ALMA"
__version__="0.4.0@2016.09.01"



import sys
sys.path.insert(0,'/home/stephane/git/signalanalysis/SignalAnalysis/Wavelet/')
sys.path.insert(0,'/home/stephane/workspace/AIV/science/analysis_scripts/')

import os
import os.path
import shutil

import numpy as np
import pylab as pl
# from scipy import signal
# from scipy.optimize import curve_fit
# import math 
# import wavelet as wav
import analysisUtils as aU
from casa import casac
from casa import split

import sqlite3 as sql

ARCSECTORAD = 4.84813681109536e-06

SOLARSYSTEMOBJ = ["Titan", "Saturn", "Neptune", "Callisto","Mars","Mercury","Venus","Uranus"]

tb = casac.table()
msmd = casac.msmetadata()

class calibrator:
    
    def __init__(self):
        
            self.listIntent = ["CALIBRATE_BANDPASS","CALIBRATE_FLUX","CALIBRATE_AMPLI","CALIBRATE_PHASE"]
            
            
    def splitCalibrator(self, msName):
        "Split the MS calibrators"
                
        listCalibratorMS = []
        
        es = aU.stuffForScienceDataReduction()
        
        try :
            tb.open(msName)
            dataColNames = tb.colnames()
            tb.close()
            
        except:
            print("###")
            print("### Error to open MS : %s \n"%(msName))
            return(-1)

        ## Band
        
        msmd.open(msName)
        chan_freq = msmd.chanfreqs(1) 
            
        band =  aU.freqToBand(chan_freq[0])
                
        print("## band %s"%(band))
        
        



        if 'CORRECTED_DATA' in dataColNames:
            dataCol = 'corrected'
        else:
            dataCol = 'data'

        calNameList = []
        
        intentSources = es.getIntentsAndSourceNames(msName)
        
        for intentkey in self.listIntent:
                   
            calIds   = intentSources[intentkey]['sourceid']
            calName  = intentSources[intentkey]['name']
            # calIds   = sorted(dict.fromkeys(calIds).keys())
            # calName  = sorted(dict.fromkeys(calName).keys())
            
            # print calName
            # print calIds
            
            iterCalid = iter(calIds)

            for name  in calName:
                
                sourceid =  iterCalid.next()
                
                if name != '':
                    found = False
                else :
                    found = True                        ## avoid empty name
                for checkname in calNameList:
                    if checkname['name'] == name :
                        found = True
                
                if not found :
                    calNameList.append({'name' : name, 'intent': intentkey, 'sourceId': sourceid})
                    
        
        print calNameList
        
        for cal in calNameList:
                           
        ## split the calibrators
            nameMSCal = "%s-%s-%s.ms"%(msName, cal['intent'], cal['name'])
            
            print("## splitting %s \n"%(nameMSCal))
            if os.path.exists(nameMSCal):
                shutil.rmtree(nameMSCal)
            

            success = split(vis=msName, outputvis= nameMSCal,datacolumn= dataCol,field= cal['name'],spw="",width=1,antenna="",
                      timebin="0s",timerange="",scan="",intent="",array="",uvrange="",correlation="",observation="",combine="",
                      keepflags=True,keepmms=False)
            
            print success
            
            if success :
                   refdir = msmd.refdir(cal['sourceId'])
                   listCalibratorMS.append([nameMSCal,cal['name'], band, refdir])
                   print("## direction:")
                   print(refdir)
                   print("#")
                   
                
            else:
                print("### Splitting : Error with %s"%(nameMSCal))
                
        msmd.close()
                   
        return(listCalibratorMS)
    

class dbstore:
    
    def __init__(self,dbname):
        "Connect to dbname"
        
        self.db = dbname      
        self.checkifDB()
    
    def checkifDB(self):
        "Check if the DB exists, if not, it is created with the tables"
        
        if not os.path.isfile(self.db):
            self.createDBTables() 

            
    def createDBTables(self):
        "Create the Tables (in sqlite3) if DB does not exist"
        
        conn = sql.connect(self.db)
        c = conn.cursor()
        
        ## Create the tables
        c.execute('''CREATE TABLE dataset
             (dataid INTEGER PRIMARY KEY , calid int, msfile text, calibrator text, band int, refer text, m0 float, m1 float)''')
                      

        conn.commit()
        conn.close()          
    
    
    def storeMScalibrator(self,item):
        "Store the dataset in the DB"
        
        msName     = item[0]
        calibrator = item[1]
        band       = item[2]
        refdir     = item[3]
        
        m0 = refdir['m0']['value']
        m1 = refdir['m1']['value']
        
        toldistcalrad = item[4] *  ARCSECTORAD
        
        
        #############      
        conn = sql.connect(self.db)
        c = conn.cursor()
        
        ### search for the id calibrator using the position. If not found add a new calibrator id.
        ###
    
        cmd1= "SELECT calid  FROM dataset  WHERE abs(m0 - %f) < %f  AND abs(m1 - %f) < %f "%(m0 , toldistcalrad, m1 , toldistcalrad)
 
        ## check for solar system objects
        if calibrator in SOLARSYSTEMOBJ:
            cmd1= "SELECT calid  FROM dataset  WHERE calibrator = '%s' "%(calibrator)

        print cmd1
        
        c.execute(cmd1)
        rows = c.fetchall()
        
        
        if len(rows) == 0. :
            cmd2 = "SELECT max(calid) FROM dataset "
            c.execute(cmd2)
            maxId = c.fetchone()[0]
            if maxId is None :
                    maxId = 0
            idCal = maxId + 1
            print("## Insert new Calibrator in the DB, ID = %d"%(idCal))
        else :
            idCal = rows[0][0]
            
        print idCal
                
        
        
        c.execute("INSERT INTO dataset(msfile,calid , calibrator,band, refer, m0, m1 ) VALUES('%s',%d, '%s',%d, '%s', %f, %f)"%(msName, \
                 idCal, calibrator, band, refdir['refer'], m0 , m1))
        

        conn.commit()
        conn.close()
        return(0)
    
    
    def checkForCalibrator(self,refdir, toleranceDist):
        "Check if the calibrator exists in the DB and return the current name"
        
        toldistcalrad = toleranceDist *  ARCSECTORAD
        
        m0 = refdir['m0']['value']
        m1 = refdir['m1']['value']

        #############      
        conn = sql.connect(self.db)
        c = conn.cursor()
        
        ### search for the id calibrator using the position. If not found add a new calibrator id.
        ###
    
        cmd1= "SELECT calid , calibrator  FROM dataset  WHERE abs(m0 - %f) < %f  AND abs(m1 - %f) < %f "%(m0 , toldistcalrad, m1 , toldistcalrad)
        print cmd1
        
        c.execute(cmd1)
        rows = c.fetchall()
        
        currentName  = ""
        foundCal = False
        print rows
        
        if len(rows) != 0:
            foundCal = True
            currentName = rows[0][1]
            print("### currentName %s"%(currentName))
         
        conn.commit()
        conn.close()
            
        return(foundCal, currentName)
    
    
 
class calStructure:
    """
    Class to organize the calibrator MS from the split files. Use a sqlite DB to keep the information.
    INPUT:
        stuctureFile: parametre files with the different direcoties and the DB
        listMSFile : list of MS to be split and organized
    """
    
    def __init__(self, structureFile, listMSFile ):
        
        self.structureFile = structureFile
        self.listMSfile    = listMSFile
        
        self.ROOTDIR  = "./"
        self.DBNAME   = "test.db"
        self.MSRM     = False
        self.TOLDISTCAL = 10.0         # maximum distance in arcsec to match one calibrator ID (default)
        self.__readStructureFile()
        
    
    def __readStructureFile(self):
        
        f = open(self.structureFile,"r")
        
        for data in f:
            dataSpl = data.split()
            
            if len(dataSpl) > 1:
                if dataSpl[0] == 'ROOTDIR':
                    self.ROOTDIR = dataSpl[2]
                    if self.ROOTDIR[-1] != "/" :
                        self.ROOTDIR += "/"
                    
                if dataSpl[0] == 'DBNAME':
                    self.DBNAME = dataSpl[2]
                    
                if  dataSpl[0] == 'MSRM':
                    if  dataSpl[2] == "False":
                        self.MSRM     = False  
                    elif dataSpl[2] == "True":
                        self.MSRM     = True
                        
                if dataSpl[0] == 'TOLDISTCAL':
                    self.TOLDISTCAL = float(dataSpl[2])        
                    
                        
        f.close()
        
    def moveMSdirectory(self,listCalMS):
        "Move the split  MS with the Cal name. If the directory does not exist it will create it"
        
        print("moving..............")
        
        ## open the DB
        db = dbstore(self.DBNAME)
        
        
        
        for calms in listCalMS:
            
            print calms
            
            band           = calms[2][0]
            calibratorName = calms[1]
            refdir         = calms[3]
            destdir = self.ROOTDIR + calms[1] + '/' + "Band%d"%(band) + '/'

            
            ## check if it is a new calibrator or not to keep the same name...
            ## 
            print("## check if calibrator exists ..")
            
            foundCalibrator , currentName = db.checkForCalibrator(refdir, self.TOLDISTCAL)
            
            if foundCalibrator:
                destdir = self.ROOTDIR + currentName + '/' + "Band%d"%(band) + '/'

            
            print("destdir")
            print destdir
            
            if not os.path.exists(destdir):
                print("### Create a new directory %s"%(destdir))
                os.makedirs(destdir)
                
            msdest = destdir + calms[0].split('/')[-1]
            print"msdest"
            print msdest
            
            if os.path.exists(msdest):
                print("### %s exists already. We skip it."%(msdest))
                shutil.rmtree(calms[0])
                
                
            else:
                shutil.move(calms[0],msdest)
                
                itemDB = [msdest, calibratorName, band, refdir, self.TOLDISTCAL]
                db.storeMScalibrator(itemDB)
                
 
           
    def run(self):
        "run the extraction and tree-organization"
        
        cal = calibrator()
        
        ## check and create the DB if needed
        db = dbstore(self.DBNAME)
        db.checkifDB()
        
        ## go thru the list of MS..
        

        
        for msfile in self.listMSfile:

            
            #try :
            ## extract the calibrator MS from each MS
        
            print("## splitting calibrators in %s"%(msfile))
            listCalMS = cal.splitCalibrator(msfile)
                
                
            ## move the calibrator MS to the corresponding place

            if listCalMS != -1:
                self.moveMSdirectory(listCalMS)
        
        
        
            
            ## remove the original MS if it proceeds.
            
            
            #except:
            #    print("##")
            #    print("## Error to move split file from  MS: %s"%(msfile))
            #    print('## \n\n')
        
        

        

class imagingCalibrator:
    """
    class providing methods to image calibrators.
    """
    
    def __init__(self):
        
        pass
    
    
    def continuum(self):
        "imaging the continuum"

        pass
    
        
    def line(self):
        "imaging a line after continuum subtraction"
        
        pass

    
##############################################
############### Main program #################
if __name__ == "__main__":
    
    ## Testing...
    
    pass

            



        
