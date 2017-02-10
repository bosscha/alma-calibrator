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
        
        
    2017.02.09:
        - add the main part to be used standalone
        
    2017.02.10:
        - update the statistics to include the total detected line

RUN:



"""

__author__="S. Leon @ ALMA"
__version__="0.1.0@2017.02.09"


import calibratorLines as cL
import numpy as np
import numpy.random as rd
from scipy import signal
import time


class simulSpectra:
    
    
    def __init__(self):
        
        self.nFreq      = 10000
        self.nLines     = 20
        self.ARange     = [2.,10.]    # absolute but can be in absorption
        self.muRange    = [500,9500]      # integer
        self.sigmaRange = [0.1 , 20.0]
        self.noise      = 1.0
        
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
    
    def __init__(self,  filepar = 'extractLine.par'):
        
        simulSpectra.__init__(self)
        self.nSample = 10
        self.result  = []
        self.filepar = filepar
        
    def generateExtraction(self):
        
        
        self.run()
        
        al = cL.analysisSpw("sample", filepar = self.filepar)
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

        LineDetectedTotal = 0
        FalseLineTotal = 0.0
        
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
                    LineDetectedTotal += 1.0
                else :
                    totalLineSnr[ind] += 1.
                    FalseLineTotal += 1.0 
                    
        
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
                    
         
        totalLine = self.nSample * self.nLines
                   
        fractionDetected  = LineDetectedTotal / totalLine
        fractionFalseLine = FalseLineTotal / totalLine
                    
                    
        
        return([snrArr, trueLineSnr], [snrArr, LineDetectedSnr], fractionDetected, fractionFalseLine)            
            

    
    def   isLineReal(self, linefound, lines) :
        "check if the line found is real"
        
        freq1 = linefound[1]
        freq2 = linefound[2]
        freq  = (freq1 + freq2) / 2.
        
        found = False
        TolFreq = 0.3      # sigma
        
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
        TolFreq = 0.3   # sigma
        
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
        
    
##############################################
############### Main program #################
if __name__ == "__main__":
    
    ## run a simulation
       
    a =  extraction()
    stat = a.analysisSample()
    
    print("Statistics ...")
    print stat
        
        
        
