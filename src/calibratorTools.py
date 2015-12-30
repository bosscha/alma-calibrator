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
        

RUN:

"""


__author__="S. Leon @ ALMA"
__version__="0.1.1@2015.12.28"



import sys
sys.path.insert(0,'/home/stephane/git/signalanalysis/SignalAnalysis/Wavelet/')
sys.path.insert(0,'/home/stephane/workspace/AIV/science/analysis_scripts/')

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
            split(vis=msName, outputvis= nameMSCal,datacolumn= dataCol,field= cal['name'],spw="",width=1,antenna="",
                timebin="0s",timerange="",scan="",intent="",array="",uvrange="",correlation="",observation="",combine="",
                keepflags=True,keepmms=False)
            
            listCalibratorMS.append([msName,cal['name']])           
        
        
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
             (dataid INTEGER PRIMARY KEY , msfile text, filedata text, calibrator text, lines int)''')
        
    
        c.execute('''CREATE TABLE calibrator
             (calid INTEGER PRIMARY KEY,dataset_id INTEGER, source text ,FOREIGN KEY(dataset_id) REFERENCES dataset(dataid) )''')
             

        conn.commit()
        conn.close()          
    
 
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
        self.__readStuctureFile()
        
    
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
        
        
    def run(self):
        "run the extraction and tree-organization"
        
        ## check and create the DB if needed
        db = dbstore(self.dbname)
        
        ## go thru the list of MS..
        
        
        ## extract the calibrator MS from each MS
        
        
        ## move the calibrator MS to the corresponding place
        
        
        ## remove the original MS if it proceeds.
        

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

            



        