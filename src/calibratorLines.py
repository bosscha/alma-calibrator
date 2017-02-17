"""
This program aims to detect lines in ALMA calibratos:
- Extract the spw from the calibrators
- analyze the spectra (wavelet filtering, etc)
- give a scoring from the presence of lines
- produce plots of the spw

HISTORY:
    2014.06.05:
        - first shot
        
    2014.06.06:
        - create report
        - create analysisSpw
        - test the workflow w/o plotms
        - create the searchLines
        - detect lines, channels and amp
        
    2014.06.07:
        - create the line report
        - add S/N
        - plot optional
        - work with CASA if astropy is not imported
        - fix some CASA issues
        
    2014.06.09:
        - extract CASA data
        
        
    2014.06.10:
        - working directory
        
    2014.06.11:
        - update plotting
        
    2014.06.17:
        - create an artificial spectra with lines  and export to  a file.
        
        
    2014.06.19:
        - check if the number of channels is okay for the filtering. Otherwise, no filtering
        - change the scale of the plots.
        - plot the lines if detected
        
    
    2015.01.23:
        - add the Gaussian fitting for a line (lineFitting)
        
    2015.01.25:
        - Adding the channel and noise and noise filtered in the result
        
        
    2015.10.02:
        - Revising calibratorLine to make it more robust vs. the data
        - add information about the spectra in the report
        - add the option DIRPLOT in the extractline.par file
        
    2015.10.05:
        - add the option DIRDATA in the extractline.par file. Note that now the data file will be redone if it is again analyzed.
        
    2015.10.21:
        - update frequency resolution
        
    2015.10.30:
        - store the results in a database. The API is in the class dbstore. the first implementation with sqlite
        - add the full ms name in listSpwSuccess to uniquely identify the spw with spwName+msName
        
    2015.10.31:
        - add DBSTORE keyword to store in a db if present.
        - create the Table dataset and populate it with real data
        
    2015.11.01:
        - add the contrast for the lines (Amplitude / continuum)
        
    2015.11.02:
        - remove the duplicated entries for dataset/lines
    
    2015.11.05:
        -fix the lines table. Note that with plotms the minimum frequency resolution is 1MHZ
        -add contrast
        
    2015.11.09:
        -fix the database ingestion
        
    2015.11.10:
        - many fixes to make the extraction and analysis more robust.
        - add source to lines DB
        
        
    2015.11.16:
        - add column in lines about the molecular specie  and status of the line (likely real ?, edge, etc..). Default is Null
        - add the maxChannel column to eliminate the edges.
        
    2015.11.21:
        - check if the calibrator name has two or more names separated by ; and take the first one for the source
        
    2015.11.23:
        - change some column of lines table (flag instead of status)
        - exception if line fitting fails.
     
    2015.12.27:
        - adding spw column in dataset table and updating the ingestion.
        
    
    2017.02.02 :
        - modifying analysisSpw to add filepar for the parameter file (default extractLine.par)
        - Add the field coordinates the dataset table (for easy search later). The extract SPW was accordingly modify
        
    2017.02.09:
        - change the method to get the spwIds
        
    2017.02.10:
        - comment the sys.path for wavelet to use the local one
        
    2017.02.17:
        - add an exception in extractSpw

RUN:
"""

from os.path import curdir

__author__="S. Leon @ ALMA"
__version__="0.7.0@2017.02.17"


import sys
#sys.path.insert(0,'/home/stephane/git/signalanalysis/SignalAnalysis/Wavelet/')
sys.path.insert(0,'/home/stephane/workspace/AIV/science/analysis_scripts/')
sys.path.insert(0,'/users/sleon/AIV/science/analysis_scripts')

import os
import os.path
import shutil

import numpy as np
import pylab as pl
from scipy import signal
from scipy.optimize import curve_fit
import math 
import wavelet as wav

import analysisUtils as aU
from casa import casac
from casa import plotms

import sqlite3 as sql

tb = casac.table()


RAD2DEGREES = 180. / math.pi 

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
             (dataid INTEGER PRIMARY KEY , msfile text, filedata text, calibrator text, coordSky1 real, coordSky2 real, spw int, lines int)''')
        
    
        c.execute('''CREATE TABLE lines
             (lineid INTEGER PRIMARY KEY,dataset_id INTEGER, source text, channelnumber int, \
            freq1  REAL, freq2  REAL, amplitude real, sn real, maxChannel int, chan1 int, chan2 int, noise real,noiseFiltered real,\
            contrast real, A_fit real, mu_fit real, sigma_fit real , molecule text, flag text, comment text, FOREIGN KEY(dataset_id) REFERENCES dataset(dataid) )''')
             

        conn.commit()
        conn.close()
        
           
    def storeMSCalSpwLines(self,arg):
        "Store the MS+spw+Cal and lines"
        
        conn = sql.connect(self.db)
        c = conn.cursor()
        
        msName     = arg[0]
        calibrator = arg[1][0]
        coord1     = arg[1][1]
        coord2     = arg[1][2]
        fileData   = arg[2]
        spwindow   = arg[3]
        lines      = arg[4]
        maxChannel = arg[5]
        
        ############
        ## check if the dataset exists, if yes remove the old ones...
        ############
        
        t = (msName, fileData)
        c.execute('SELECT * FROM dataset WHERE msfile=? AND filedata=?', t)
        duplication = c.fetchall()
        
        if len(duplication)>0:
            # remove it and eventually remove the lines
            for data in duplication:
                id = data[0]
                c.execute('SELECT * FROM lines WHERE dataset_id=?', (id,))
                lineduplicated = c.fetchall()
                print lineduplicated
                
                if len(lineduplicated)>0:
                    for li in lineduplicated:
                        line_id = li[0]
                        c.execute('DELETE FROM lines WHERE lineid=?', (line_id,))
                        
                c.execute('DELETE FROM dataset WHERE dataid=?', (id,))
                
        #############      
        
        c.execute("INSERT INTO dataset(msfile,filedata,calibrator, coordSky1, coordSky2,spw,lines) VALUES('%s','%s','%s',%f, %f, %d,%d)"%(msName , fileData, calibrator, coord1, coord2, spwindow, len(lines)))
        
        datasetid = c.lastrowid
        
        if len(lines) > 0:
            for li in lines:
                nc = li[0]
                f1 = li[1]
                f2 = li[2]
                am = li[3]
                sn = li[4]
                line0 = li[5]
                line1 = li[6]
                noise = li[7]
                noiseFilt = li[8]
                contrast  = li[9]
                a_fit     = 0.
                mu_fit    = 0.
                sigma_fit = 0.
                
                
                print"li"
                print len(li)
                
                if len(li) == 12:
                    a_fit  = li[10][0]
                    mu_fit = li[10][1]
                    sigma_fit = li[10][2]
                
                c.execute("INSERT INTO lines(dataset_id, source, channelnumber, freq1, freq2, amplitude, sn, maxChannel, chan1, chan2, noise,noiseFiltered,\
                    contrast, A_fit, mu_fit, sigma_fit) VALUES(%d, '%s' ,%d,%f, %f,%f,%f, %d,  %d, %d, %f, %f, %f, %f, %f, %f)"%(datasetid, \
                    calibrator, nc,f1,f2,am,sn, maxChannel, line0, line1, noise, noiseFilt, contrast, a_fit, mu_fit, sigma_fit))
                
        
        print("datasetid %d"%(c.lastrowid))

        conn.commit()
        conn.close()
        return(0)
    

        
      
        
        
        
class extractSpwField:
    
    def __init__(self,listFileMS, workingdir ='./'):
        
        self.listFileMS = listFileMS
        self.report = ""
        self.listIntent = ["CALIBRATE_BANDPASS","CALIBRATE_FLUX","CALIBRATE_AMPLI","CALIBRATE_PHASE"]
        self.listSpwSuccess = []
        self.workingDirectory = workingdir
        
        
    def extractSpw(self,msName,intent = 'CALIBRATE_BANDPASS', listCal = []):
        "Extract the spectra from a spw of a given intent"
        

        es = aU.stuffForScienceDataReduction()
        # vm = aU.ValueMapping(msName)

        try:        
            tb.open(msName)
            dataColNames = tb.colnames()
            tb.close()
            
        except:
            print("##")
            print("## Error to open : %s"%(msName))
            print('## \n\n')
            return(False,[], [],[])
            

        if 'CORRECTED_DATA' in dataColNames:
            dataCol = 'corrected'
        else:
            dataCol = 'data'

        ### Extraction information MS
        
        try:
            tb.open(msName+'/FIELD')
            sourceIds = tb.getcol('SOURCE_ID')
        
            ##### new coordinates
            srccoords = tb.getcol('REFERENCE_DIR')

            ######
        
            tb.close()
                    
            intentSources = es.getIntentsAndSourceNames(msName)
            
            tb.open(msName+'/DATA_DESCRIPTION')
            spwIds = tb.getcol('SPECTRAL_WINDOW_ID').tolist()
            tb.close()
            
        except:
            print("##")
            print("## Error (FIELD/SPW) to open : %s"%(msName))
            print('## \n\n')
            return(False,[], [],[])        


     
        
        if  intentSources[intent]['name'][0] == '':
            return(False, [], [], [])

        
        calIds   = intentSources[intent]['sourceid']
        calName  = intentSources[intent]['name']
              
        print calIds
        print calName
         
           
        
        print spwIds
        
        listFile    = []
        listCalFile = []
        listSuccess = []
        listSpw     = []
        
        indexOnCalId = 0
        
        for calNameJ in calName:
            
            ## coordinates in degrees
            print calNameJ
            print calIds[indexOnCalId]
            sourceId = calIds[indexOnCalId]
            coordSky1 = ( srccoords[0][0][sourceId] % (2 * math.pi) ) * RAD2DEGREES
            coordSky2 = srccoords[1][0][sourceId]  * RAD2DEGREES
            print coordSky1
            print coordSky2
            
            indexOnCalId += 1
                   
            
            found = False
            if calNameJ in listCal:
                found = True

            
            if not found:
                listCal.append(calNameJ)
            
                for spwj in spwIds:
                    # spwChanFreqs = vm.spwInfo[j]['chanFreqs']  
                    # print spwChanFreqs
            
                    # print("Test coordinates")
                    
 
                    spwId = "%s"%(spwj)
                
                    msNameWithoutSlash = msName.split('/')[-1]
                    msNameDat = self.workingDirectory+msNameWithoutSlash
                    calNameJwithoutmultiplename = calNameJ.split(";")[0]
                    
                    fileName = "%s-%s-%s-spw%s.dat"%(msNameDat,intent,calNameJwithoutmultiplename,spwj)
                    print fileName
            
                    successSpw = plotms(vis = msName ,xaxis="frequency",xdatacolumn = "",yaxis="amp",ydatacolumn = dataCol ,selectdata=True,field = calNameJ,
                                spw = spwId ,timerange="",uvrange="",antenna="",scan="",correlation="",array="",observation="",msselect="",averagedata=True,
                                avgchannel="1",avgtime="1e9",avgscan=True,avgfield=True, avgbaseline=True,avgantenna= False,avgspw=False,scalar=False,
                                transform=False,freqframe="",restfreq="",veldef="RADIO",shift=[0.0, 0.0],extendflag=False,extcorr=False,extchannel=False,
                                iteraxis="",plotfile = fileName ,expformat="txt",
                                highres=False, overwrite=False,showgui=False)
            
                    listFile.append(fileName)
                    listCalFile.append([calNameJwithoutmultiplename, coordSky1, coordSky2] )
                    listSpw.append(spwj)
                    listSuccess.append(successSpw)
                    
        
                    if not successSpw :
                        print("## Problem with the extraction of %s ..."%(fileName))
        
        
        print("vvvvvvvvvvvvvvvvvvvvvvvv")
        print listFile
        print("vvvvvvvvvvvvvvvvvvvvvvv")
        
        return(listSuccess, listFile, listCalFile, listSpw)
    
        
    def extractListMS(self):
        "Read in the file the list of MS"
        
        
        f = open(self.listFileMS)
        
        self.listMS = []
        
        for line in f:
            if line[0] != "#":
                dat = line.split()
                self.listMS.append(dat[0])
                
                
            
            
    def extractSpwIntent(self):
        "Extract the different Intent"
        
        
        ## loop on the MS
        for ms in self.listMS:
            
            
            for intent in self.listIntent:
                
                # print("# Intent: %s \n"%(intent))     
                listCal = []
                successList = []
                fileList = []
                listSpw  = []
                
                
                successList , fileList , listCal , listSpw = self.extractSpw(ms, intent, listCal)
                
                print("dddddddddddddddddddddddddd")
                print fileList
                print("dddddddddddddddddddddddddddddddddddddddddddddddddddd")

                
                print("##start ..")
                
                
                if len(fileList) > 0:
                    print("# %s - %s : ok"%(ms, intent))
                    print("# %s created "%(fileList))
                    
                    index = 0
                    for filename in fileList:
                        calibrator = listCal[index]
                        spw        = listSpw[index]
                        
                        
                        if successList[index]:
                            self.listSpwSuccess.append([filename,ms, calibrator, spw])


                        index += 1
                    
                else :
                    print("# %s - %s : failed"%(ms, intent))
                    
            
     
    def createReport(self):
        "Create the report on the SPW extraction"
        
        pass
           
        
    def run(self):
        "Run the extraction of the spw  / field"    
        
        self.extractListMS()
        self.extractSpwIntent()
        self.createReport()
        

        return(self.listSpwSuccess)

#######################

class analysisSpw:
    
    def __init__(self,file, filepar = "extractLine.par"):
        
        self.DataFile = file
        self.parameterFile = filepar
        
        self.DIRPLOT = ''
        self.DIRDATA = ''
        self.DBNAME  = ''
        self.readParameterSearch()
    
    def readParameterSearch(self):
        "Extract the parameter for the Search"
        
        f = open(self.parameterFile,"r")
        
        for data in f:
            dataSpl = data.split()
            
            if len(dataSpl) > 1:
                if dataSpl[0] == 'WAVSN':
                    self.WAVSN = float(dataSpl[2])
                
                if dataSpl[0] == 'WAVSCALE':
                    self.WAVSCALE = int(dataSpl[2])
                    
                if dataSpl[0] == 'THRESHOLDLINE':
                    self.THRESHOLDLINE = float(dataSpl[2])
                    
                if dataSpl[0] == 'NUMCHANDETECTMIN':
                    self.NUMCHANDETECTMIN = int(dataSpl[2])
                
                if dataSpl[0] == 'SHOWPLOT':
                    self.SHOWPLOT = dataSpl[2]
                    
                if dataSpl[0] == 'FITLINE':
                    self.FITLINE =  dataSpl[2]
                    
                if dataSpl[0] == 'DIRPLOT':
                    self.DIRPLOT = dataSpl[2]
                    
                if dataSpl[0] == 'DIRDATA':
                    self.DIRDATA = dataSpl[2]
                    
                if dataSpl[0] == 'DBNAME':
                    self.DBNAME = dataSpl[2]               
                
                
        
        f.close
        

    def  searchLines(self, freq, amp):
        """
        Main script to search for line with a scoring at the end.
        The parameters for the search are in self.parameterFile
        Return the lines parameters detected. (channels, amplitude)
        """
        
        curdir = os.getcwd()
        
        if self.DIRPLOT != '':
            if not os.path.isdir(self.DIRPLOT):
                os.mkdir(self.DIRPLOT)
                os.chdir(self.DIRPLOT)
            else :
                os.chdir(self.DIRPLOT)
                
        
        if self.SHOWPLOT == 'Y':
                      
            
            yminA = min(amp)
            ymaxA = max(amp)
            ymax = ymaxA + (ymaxA-yminA) *1.1
            ymin = ymaxA - (ymaxA-yminA) *2.5
            
            fig = pl.figure(1, figsize=(14,9))
            pl.clf()
            figfile = self.DataFile+".raw.png"
            pl.title(figfile.split('/')[-1])
            pl.xlabel("Frequency (GHz)")
            pl.ylabel("Amplitude")
            pl.xlim([freq[0],freq[-1]])
            pl.ylim([ymin,ymax])
            pl.plot(freq, amp)
            fig.savefig(figfile)
            pl.draw()
        
       
        ## filtering with the wavelet
           
        ##
        
        noise = np.std(amp)
        mean  = np.mean(amp)
        print mean
        
        nChannels = len(amp)
        nChanP2 = math.log(nChannels) / math.log(2.0)
        
        
        wt = wav.wt()        
        
        
        if nChanP2 > self.WAVSCALE + 2 :   
            wspw = wt.atrous1d(amp, self.WAVSCALE)
            wspwFilt = wt.filtering1d(wspw, self.WAVSN, waveletNoise = True, spectralNoise = noise)
            spwRestore = wt.restore1d(wspwFilt, 0, self.WAVSCALE)
            
        else :
            print("# No filtering because of the number of channels is too low for the WT")
            spwRestore = amp
            
            
            
        ## new noise
        noiseFiltered = np.std(spwRestore)
        print("Noise : %f \n"%(noise))
        print("Noise filtered : %f \n"%(noiseFiltered))
        
        ampZeroed = np.abs(spwRestore - mean)
        indAboveThreshold = np.where(ampZeroed > self.THRESHOLDLINE * noiseFiltered )
    
        ## number of contiguous channel
        ##
        numOfChannel = 0
        
        for i in range(len(indAboveThreshold[0])-1):
            dChan = abs(indAboveThreshold[0][i+1] - indAboveThreshold[0][i])
            if dChan == 1:
                numOfChannel += 1
                
        lines = []
        
        if numOfChannel > self.NUMCHANDETECTMIN:
            lines = self.detectLineChannel(indAboveThreshold[0], ampZeroed , self.NUMCHANDETECTMIN)   
            print("# Lines Detected: %d"%(len(lines)))
            print("# Num. Total of Channels: %d"%(numOfChannel))
            print("# Line Channels:")
            print lines

        finalLines = []
        if len(lines) == 0:
            
            if self.SHOWPLOT == 'Y':
            
                yminA = min(spwRestore)
                ymaxA = max(spwRestore)
                ymax = ymaxA + (ymaxA-yminA) *1.1
                ymin = ymaxA - (ymaxA-yminA) *2.5
            
            
                fig2 = pl.figure(2)
                pl.clf()
                figfile = self.DataFile+".filt.png"
                pl.title(figfile.split('/')[-1])
                pl.xlabel("Frequency (GHz)")
                pl.ylabel("Amplitude")
                pl.xlim([freq[0],freq[-1]])
                pl.ylim([ymin,ymax])
                pl.plot(freq, spwRestore)
                fig2.savefig(figfile)
                pl.draw()    
            
            os.chdir(curdir)
            return(finalLines)
        
        
        else :
            lineTemp = []
            
            ampMean = spwRestore - mean
            
            for line in lines:
                if abs(min(ampMean[line[0]:line[1]])) > max(ampMean[line[0]:line[1]]) and min(ampMean[line[0]:line[1]]) < 0. :
                    am = min(ampMean[line[0]:line[1]])
                else :
                    am = max(ampMean[line[0]:line[1]])
                             
                
                nC = line[1]-line[0]+1
                f1 = freq[line[0]]
                f2 = freq[line[1]]
                # am = line[2]
                sn = abs (am / noiseFiltered)
                contrast = am / mean
                
                if self.FITLINE == 'Y':
                    fit, coeff = self.lineFitting(freq , amp - mean , f1 , f2 , am)
                    finalLines.append([nC, f1,f2 ,am , sn , line[0], line[1], noise, noiseFiltered, contrast, fit , coeff])
                else :
                    finalLines.append([nC,f1,f2 ,am , sn , line[0], line[1], noise, noiseFiltered, contrast ])
                
        
        
            if self.SHOWPLOT == 'Y':
            
                yminA = min(spwRestore)
                ymaxA = max(spwRestore)
                ymax = ymaxA + (ymaxA-yminA) *1.1
                ymin = ymaxA - (ymaxA-yminA) *2.5
            
                print "draw line"
                fig2 = pl.figure(2)
                pl.clf()
                figfile = self.DataFile+".filt.png"
                pl.title(figfile.split('/')[-1])
                pl.xlabel("Frequency (GHz)")
                pl.ylabel("Amplitude")
                pl.xlim([freq[0],freq[-1]])
                pl.ylim([ymin,ymax])
                pl.plot(freq, spwRestore)
                # plot lthe lines
                for res in finalLines:
                    xx = [res[1],res[2]]
                    yline = ymin + (ymaxA-yminA) * 0.15
                    yy = [yline,yline]
                    pl.plot(xx, yy,"r-", linewidth = 2.5)
                
                    fig2.savefig(figfile)
                    pl.draw()    
            
            
        os.chdir(curdir)      
        return(finalLines)
                


            
    def   detectLineChannel(self,indexLines, amp, minChan):
        "Detect the range of the line and the amplitude normalized (+)  with minChan"
        
        lines =[]
        
        chan1 = indexLines[0]
        chan2 = indexLines[0]

        ampmax = amp[chan1]
        
        lineEnd = False
        
        for index in indexLines: 
            if abs(index - chan2) <= 1:
                chan2 = index
                if amp[chan2] > ampmax :
                    ampmax = amp[chan2]
                    
            else :
                lineEnd = True
                
            
            if abs(chan2 - chan1) >= minChan and lineEnd :
                lines.append([chan1,chan2,ampmax])
            
            if lineEnd :
                lineEnd = False
                chan1 = index
                chan2 = index
                ampmax = amp[chan1]
            
        
        ## evaluate the last values        
        if abs(chan2 - chan1) >= minChan:
                lines.append([chan1,chan2, ampmax])       
        
        return(lines)
                    
    def gauss(self, x, *p):
        "Gaussian model for fitting"
        
        A, mu, sigma = p
        return A * np.exp(-(x-mu)**2/(2.*sigma**2))
    
    
    def lineFitting(self, freq, amp, freq1, freq2, ampmax):
        """
        Gaussian fitting on the detected line
        freq1 , freq2 are the extrema frequencies for the detection
        ampmax is the maximum amplitude 
        """
        
        A     = ampmax 
        mu    = (freq1 + freq2) / 2.
        sigma = abs((freq2 - freq1)  / 2.) 
        p0 = [A, mu, sigma]
        
        try:
            fitParams, fitCovariances = curve_fit(self.gauss , freq , amp , p0 = p0)
        
        except:
            print("### Fit Line error.")
            fitParams = [0.,0.,0.]
            fitCovariances = [0., 0., 0.]
            
        return(fitParams , fitCovariances)
        
    
    def extractData(self,dataFile,casa = False):
        "Extract the data from the SPW file, special case for CASA file"
        
   
        if not casa:
            f = open(dataFile)
        
            nData = 0
            freqtemp = []
            amptemp = []
        
            for line in f:
                data = line.split()
                freqtemp.append(float(data[0]))
                amptemp.append(float(data[1]))
                nData += 1
            
            freq = np.zeros(nData)
            amp = np.zeros(nData)
        
            for i in range(nData):
                freq[i] = freqtemp[i]
                amp[i]  = amptemp[i]
            
            f.close()
            
        if casa:
            
            f = open(dataFile)
            
            nData = 0
            freqtemp = []
            amptemp = []
        
            for line in f:
                data = line.split()
                if line[0] != "#" :
                    freqtemp.append(float(data[0]))
                    amptemp.append(float(data[1]))
                    nData += 1
                
            f.close()
        
            if nData == 0:
                return([],[])
            
            ftemp = np.zeros(nData)
            atemp = np.zeros(nData)
            nDat  = np.zeros(nData)
            
            ## assumes that the frequency are contiguous
            ## Note that the frequency resolution is 1MHz with plotms
            
            fcurrent = freqtemp[0]
            idat = 0
            nSpwData = 0
            
            for i in  range(nData):
                if freqtemp[i] == fcurrent:
                    ftemp[idat] = fcurrent
                    atemp[idat] += amptemp[i]
                    nDat[idat] += 1
                else :
                    nSpwData += 1
                    idat += 1
                    fcurrent = freqtemp[i]
                    ftemp[idat] = fcurrent
                    atemp[idat] += amptemp[i]
                    nDat[idat] += 1
                
            freq = np.zeros(nSpwData)
            amp = np.zeros(nSpwData)
            
    
            print nSpwData
        
            for i in range(nSpwData) :
                freq[i] = ftemp[i]
                amp[i]  = atemp[i] / nDat[i]  
            
            
            ## Sorting the frequency
            
            ## moving the dataFile in DIRDATA
            
            if self.DIRPLOT != '':
                if not os.path.isdir(self.DIRDATA):
                    os.mkdir(self.DIRDATA)
                    shutil.move(dataFile,"./%s/%s"%(self.DIRDATA,dataFile))
                    
                else :
                    shutil.move(dataFile,"./%s/%s"%(self.DIRDATA,dataFile))
            
        return(freq, amp)
                           
            
    
    def createReportLines(self, lines, nchannel, freqRes):
        "Create and return  a report on the line detection"
        
        outStr = "#\n# Lines \n"
        outStr += "# Channel number : %d \n"%(nchannel)
        outStr += "# Frequency Resolution: %6.4f MHz \n"%(abs(freqRes)*1e3)
        outStr += "# \n"
        
        if len(lines) == 0:
            outStr += "# Detected : 0 \n"
        else :
            nlines = len(lines)
            outStr += "# Detected : %d \n"%(nlines)
            outStr += "#------------- \n"
            index = 0
            for line in lines:
                outStr += "# Line : %d \n"%(index)
                outStr += "# Number of channels : %d\n"%(line[0])
                outStr += "# Start freq. (channel) : %f (%d)\n"%(min(line[1],line[2]), line[5])
                outStr += "# End   freq. (channel) : %f (%d)\n"%(max(line[1],line[2]), line[6])
                outStr += "# Linewidth (MHz)       : %6.3f\n"%(abs(line[1]-line[2]) * 1e3)
                outStr += "# Amplitude             : %f\n"%(line[3])
                outStr += "# Noise (filtered)      : %f (%f) \n"%(line[7],line[8])
                outStr += "# Contrast              : %f \n"%(line[9])
                outStr += "# S/N                   : %3.1f\n"%(line[4])
                
                if self.FITLINE == 'Y' :
                    outStr += "# Gaussian fit (A , mu, sigma)  : %f , %f , %f\n"%(line[10][0],line[10][1],line[10][2])
                    
                outStr += "#----------\n"
                index += 1
                
        return(outStr)
        
        
        
    def simulSpectra(self,fileSpectra, nBin, continuumLevel, sigma, linePositionChan, lineWidth, lineIntensity, startFreq = 0., resFreq =1.):
        """ Save in fileSpectra of nBin channels with the linePosition, lineWidth and lineIntensity list 
        starFreq and resFreq are optional"""
        
        freq = np.arange(nBin)
        spectra = np.random.normal(continuumLevel, sigma, nBin)
        
        index = 0
        for pos in linePositionChan:
            
            nChan = 4 * int(lineWidth[index] / resFreq)    
            spec = lineIntensity[index] * signal.gaussian(nChan,lineWidth[index])
            
            startPos = pos - nChan / 4
            spectra[pos:pos+nChan] = spectra[pos:pos+nChan] + spec
        
            index += 1
            
            
        f = open(fileSpectra,"w")
        
        index = 0
        for frequency  in freq :
            strOut = "%f   %f \n"%(frequency, spectra[index])
            f.write(strOut)
            index += 1
            
        f.close()
       
        
        
        
    def run(self, casa = False):
        "Run the analysis on each file"
        
        report = ""
        linesDetected = []
        nchannel = 0
                
        print("# Analyzing %s"%(self.DataFile))
        freq, amp = self.extractData(self.DataFile , casa = True)
 

        
        if len(freq) > 1:
            linesDetected = self.searchLines(freq,amp)
            nchannel = len(freq)
            freqRes  = freq[1]-freq[0]
            report = self.createReportLines(linesDetected, nchannel, freqRes)
        
        return(report,linesDetected, nchannel)
        
###################################

class extractLines:
    "Main class to extract the line"
    
    def __init__(self, inputMSfile, fileLinePar = "extractLine.par"):

        self.inputMS = inputMSfile
        self.linepar = fileLinePar
    
        
        
    def run(self):
        "Run the extraction"
        
        reportFileName = self.inputMS+".info"
        finalReport = ""
        
    ## Extract the spw
        print("## Starting the extraction of the SPW...")
        
        eSpw = extractSpwField(self.inputMS)
        listFileData = eSpw.run()
        
        
    ## Analysis of the spectra and storing the result if necessary in the db
    
        print("## Starting the analysis")
        
        print listFileData 
        
        ## check the parameters to connect if necessary to a DB
        storeondb = False
        aSpwparam = analysisSpw('dummy', filepar = self.linepar)
        
        if aSpwparam.DBNAME != '':
            db = dbstore(aSpwparam.DBNAME)
            storeondb = True

        if len(listFileData) != 0:
            for fileDataMS in listFileData:
 
                fileData   = fileDataMS[0]
                msName     = fileDataMS[1]
                calibrator = fileDataMS[2]
                spwindow   = fileDataMS[3]
                
                
                aSpw = analysisSpw(fileData, filepar = self.linepar )
                reportFileData , linesDetected , maxChannel  = aSpw.run(casa = True)
            
                if storeondb:
                    
                    db.storeMSCalSpwLines([msName,calibrator,fileData, spwindow, linesDetected, maxChannel])
                                       
                
                finalReport += "\n################# \n"
                finalReport += "# Data File : %s \n"%(fileData)
                finalReport += reportFileData
            
    
        print finalReport
          
        fout = open(reportFileName,"w")
        fout.write(finalReport)
        fout.close()
            
        print("### Done")          
        
        

        
        
    
        

##############################################
############### Main program #################
if __name__ == "__main__":
    
    ### Extraction MS using CASA
    print("## Use the  class extractLines to extract the absorption lines")


    