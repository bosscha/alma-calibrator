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

RUN:

$> casa -c scriptToOrganize.py

"""


__author__="S. Leon @ ALMA"
__version__="0.1.4@2016.01.07"



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


tb = casac.table()

class calibrator:
    
    def __init__(self):
        
            self.listIntent = ["CALIBRATE_BANDPASS","CALIBRATE_FLUX","CALIBRATE_AMPLI","CALIBRATE_PHASE"]
            
            
    def splitCalibrator(self, msName):
        "Split the MS calibrators"
                
        listCalibratorMS = []
        
        es = aU.stuffForScienceDataReduction()
        
        tb.open(msName)
        dataColNames = tb.colnames()
        tb.close()

        if 'CORRECTED_DATA' in dataColNames:
            dataCol = 'corrected'
        else:
            dataCol = 'data'

        calNameList = []
        
        intentSources = es.getIntentsAndSourceNames(msName)
        
        for intentkey in self.listIntent:
                   
            calIds   = intentSources[intentkey]['sourceid']
            calName  = intentSources[intentkey]['name']
            calIds   = sorted(dict.fromkeys(calIds).keys())
            calName  = sorted(dict.fromkeys(calName).keys())

            for name in calName:
                
                if name != '':
                    found = False
                else :
                    found = True                        ## avoid empty name
                for checkname in calNameList:
                    if checkname['name'] == name :
                        found = True
                
                if not found :
                    calNameList.append({'name' : name, 'intent': intentkey})
                    
        
        print calNameList
        
        for cal in calNameList:
                           
        ## split the calibrators
            nameMSCal = "%s-%s-%s.ms"%(msName, cal['intent'], cal['name'])
            
            print("## splitting %s \n"%(nameMSCal))
            if os.path.exists(nameMSCal):
                shutil.rmtree(nameMSCal)
            
            split(vis=msName, outputvis= nameMSCal,datacolumn= dataCol,field= cal['name'],spw="",width=1,antenna="",
                timebin="0s",timerange="",scan="",intent="",array="",uvrange="",correlation="",observation="",combine="",
                keepflags=True,keepmms=False)
            
          
            
            spwInfo = es.getSpwInfo(nameMSCal)
            spwIds = sorted(spwInfo.keys())
                       
            ff =  aU.getFrequencies(nameMSCal,spwIds[0]) 
            band = aU.freqToBand(ff[0][0])
            
            listCalibratorMS.append([nameMSCal,cal['name'], band])  
                   
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
             (dataid INTEGER PRIMARY KEY , msfile text, calibrator text, band int)''')
                      

        conn.commit()
        conn.close()          
    
    
    def storeMScalibrator(self,item):
        "Store the dataset in the DB"
        
        msName     = item[0]
        calibrator = item[1]
        band       = item[2]
        
        
        #############      
        conn = sql.connect(self.db)
        c = conn.cursor()
        
        c.execute("INSERT INTO dataset(msfile,calibrator,band) VALUES('%s','%s',%d)"%(msName, calibrator, band))
        
        conn.commit()
        conn.close()
        return(0)
        
 
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
        self.__readStructureFile()
        
    
    def __readStructureFile(self):
        
        f = open(self.structureFile,"r")
        
        for data in f:
            dataSpl = data.split()
            
            if len(dataSpl) > 1:
                if dataSpl[0] == 'ROOTDIR':
                    self.ROOTDIR = dataSpl[2]
                    
                if dataSpl[0] == 'DBNAME':
                    self.DBNAME = dataSpl[2]
                    
                if  dataSpl[0] == 'MSRM':
                    if  dataSpl[2] == "False":
                        self.MSRM     = False  
                    elif dataSpl[2] == "True":
                        self.MSRM     = True
                        
        f.close()
        
    def moveMSdirectory(self,listCalMS):
        "Move the split  MS with the Cal name. If the directory does not exist it will create it"
        
        print("moving..............")
        
        ## open the DB
        db = dbstore(self.DBNAME)
        
        
        
        for calms in listCalMS:
            band           = calms[2][0]
            calibratorName = calms[1]
            destdir = self.ROOTDIR + calms[1] + '/' + "Band%d"%(band) + '/'
            print destdir
            
            if not os.path.exists(destdir):
                print("### Create a new directory %s"%(destdir))
                os.makedirs(destdir)
                
            msdest = destdir + calms[0].split('/')[-1]
            print msdest
            
            if os.path.exists(msdest):
                print("### %s exists already. We skip it."%(msdest))
                shutil.rmtree(calms[0])
                
            else:
                shutil.move(calms[0],msdest)
                
                itemDB = [msdest, calibratorName, band]
                db.storeMScalibrator(itemDB)
                
 
           
    def run(self):
        "run the extraction and tree-organization"
        
        cal = calibrator()
        
        ## check and create the DB if needed
        db = dbstore(self.DBNAME)
        db.checkifDB()
        
        ## go thru the list of MS..
        

        
        for msfile in self.listMSfile:
            
            try:
                tb.open(msfile)
                dataColNames = tb.colnames()
                tb.close()
            
            ## extract the calibrator MS from each MS
        
                print("## splitting calibrators in %s"%(msfile))
                listCalMS = cal.splitCalibrator(msfile)
                
                print listCalMS
                
                self.moveMSdirectory(listCalMS)
        
        
            ## move the calibrator MS to the corresponding place
                print listCalMS
        
            
            ## remove the original MS if it proceeds.
            
            
            except:
                print("##")
                print("## Error to move split file from  MS: %s"%(msfile))
                print('## \n\n')
        
        

        

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

            



        