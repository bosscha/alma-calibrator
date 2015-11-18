"""
This module is designed to analyze the performances of the extraction of spectra 
through simulations.
It uses the extraction methods coded in calibratorLines.


HISTORY:
    2015.01.23:
        - create two classes
        - needs the file extractLines.par

    2015.01.25:
        - add the analysis of the line
        - inheritance of class extraction from simulSpectra
        - add absorption

RUN:



"""

__author__="S. Leon @ ALMA"
__version__="0.0.2@2015.01.25"


import calibratorLines as cL
import numpy as np
import numpy.random as rd
from scipy import signal
import time


class simulSpectra:
    
    
    def __init__(self):
        
        self.nFreq      = 10000
        self.nLines     = 20
        self.ARange     = [2.,10.]    # absolute but can be in absoprtion
        self.muRange    = [500,9500]      # integer
        self.sigmaRange = [2.0 , 20.0]
        self.noise      = 5.0
        
        self.freq = np.arange(self.nFreq)
        self.amp  = np.zeros(self.nFreq)
        
    def generateNoise(self):
        "Gaussian noise with no continuum"
        
        self.amp = rd.normal(0., self.noise, self.nFreq)
        
        
    def generateLines(self):
        "Geneerate the list of lines in the spectra "
        
        ts = int(time.time()) % 999999
        rd.seed(ts)
        
        self.Lines = []
        
        center = rd.randint(self.muRange[0], self.muRange[1], self.nLines)
        sigma  = rd.uniform(self.sigmaRange[0], self.sigmaRange[1], self.nLines)
        A      = rd.uniform(self.ARange[0], self.ARange[1], self.nLines)
        
        sign   = rd.uniform(-1.0, 1.0, self.nLines)
        
        for i in range(self.nLines):
            
            signA = 1.0
            if sign[i] < 0. :
                signA = -1.
                
            self.Lines.append({"center" : center[i] , "sigma" : sigma[i] , "A" : signA * A[i]})
            

    def generateSpectra(self): 
        " generate the spectra with the line"
        
        for line in self.Lines:
            
            print line
            
            center = line['center']
            sigma  = line['sigma']
            A      = line['A']
            
            nChan = 20 * int(sigma)
            spec = A * signal.gaussian(nChan , sigma)
            
            self.amp[center - nChan / 2  : center + nChan / 2] +=  spec
            
            
        
    def run(self):
        
        self.generateNoise()
        self.generateLines()
        self.generateSpectra()
        
###############################
class extraction(simulSpectra):
    
    def __init__(self):
        
        simulSpectra.__init__(self)
        self.nSample = 10
        self.result  = []
        
    def generateExtraction(self):
        
        
        self.run()
        
        al = cL.analysisSpw("sample")
        linesDetected = al.searchLines(self.freq, self.amp)
    
        self.result.append([self.Lines , linesDetected])
        
        
    def generateSample(self):
        "generate the list of spectra and extraction"
        
        self.result = []
        
        for i in range(self.nSample):
            self.generateExtraction()
            

    
    def linesStatistics(self):
        "Compute the statistics on the line detected"
        
        print("# analyzing the results ...")
        # compute the fraction of true line detected vs. snr
        
        nbin = 100
        trueLineSnr  = np.zeros(nbin) 
        totalLineSnr = np.zeros(nbin)

        
        dsnr   = 0.2
        snrArr = dsnr * np.arange(nbin)
        
        for spec in self.result:
            lines = spec[0]
            linesDetected = spec[1]
            
            for linefound in linesDetected:
                
                found = self.isLineReal(linefound, lines)
                snr = linefound[4]
                ind = int(snr/dsnr)
                
                if found:
                    trueLineSnr[ind]  += 1.
                    totalLineSnr[ind] += 1.
                else :
                    totalLineSnr[ind] += 1.
                    
        
        for i in range(nbin):
            if totalLineSnr[i] > 0. :
                trueLineSnr[i] = trueLineSnr[i] / totalLineSnr[i]
                
         ## analisis of the real line
         
        LineDetectedSnr  = np.zeros(nbin) 
        totalRealLineSnr = np.zeros(nbin)

          
        for spec in self.result:
                
            lines = spec[0]
            linesDetected = spec[1]
                
            for realLines in lines:
                    
                found = self.isLineDetected(realLines, linesDetected)
                    
                snr = realLines['A'] / self.noise
                ind = int(snr/dsnr)
                    
                if found:
                    LineDetectedSnr[ind]  += 1.
                    totalRealLineSnr[ind] += 1.
                else :
                    totalRealLineSnr[ind] += 1.
        
        for i in range(nbin):
            if totalRealLineSnr[i] > 0. :
                LineDetectedSnr[i] = LineDetectedSnr[i] / totalRealLineSnr[i]         
                    
                    
        return([snrArr, trueLineSnr], [snrArr, LineDetectedSnr])            
            

    
    def   isLineReal(self, linefound, lines) :
        "check if the line found is real"
        
        freq1 = linefound[1]
        freq2 = linefound[2]
        freq  = (freq1 + freq2) / 2.
        
        found = False
        TolFreq = 3.      # sigma
        
        for l in lines:
            A = l['A']
            center = l['center']
            sigma  = l['sigma']
            d = abs((freq - center) / sigma)
            
            if d < TolFreq :
                found = True
                
        return(found)
            
    def isLineDetected(sef , line , linesDetected):
        "check if a real line was detected"
        
        A = line['A']
        center = line['center']
        sigma  = line['sigma']
        
        found = False
        TolFreq = 3.   # sigma
        
        for l in linesDetected :          
            freq1 = l[1]
            freq2 = l[2]
            freq = (freq1 + freq2) / 2.
            d = abs((freq - center) / sigma)
            
            if d < TolFreq :
                found = True
                
        return(found)            
            
            
        

    def analysisSample(self):
        "Analyze the results comparing the simulated lines with the extracted ones"
        
        self.generateSample()
        stat1 = self.linesStatistics()
        
        
        return(stat1)
        
    
    

        
        
        
