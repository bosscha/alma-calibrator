#!/usr/bin/python

"""
Script to query the ALMA archive with filter to select the list of Calibrators to be analyzed


HISTORY:
    2015.09.15:
        - first shot
        - parse the list of calibrators and eliminate the duplication 
        - query the region around the calibrators for public data
        
    2015.09.16:
        - add a selectDeepfield to constrain the selection of the target
        - method to write the final report
        
    2015.09.23:
        - add maxFreqRes in selectDeepField (mainly for absorption)
        
        
    2017.03.28:
        - fixing a bung 
        
"""

__author__="S. Leon @ ALMA"
__version__="0.1.0@2015.09.23"


from astroquery.alma import Alma
from astropy import coordinates
from astropy import units as u
import time

class queryCal:
    
    def __init__(self,fileCal, fluxmin):
        
        self.listCal = self.readCal(fileCal, fluxmin)
        

    def readCal(self,file, fluxmin = 0.2):
        "Read a list of calibrators in CSV format from the Source Catalogue web interface"
        
        listcal = []
        
        fcal = open(file)
        for line in fcal:
            if line[0] != "#":
                tok       = line.split(",")
                band      = tok[0].split(" ")[1]
                flux      = float(tok[7])
                name      = tok[13].split("|")[0]
                alpha2000 = float(tok[3])
                delta2000 = float(tok[5])
                
                if flux >= fluxmin:
                    found = False
                    for nameYet in listcal:
                        if nameYet[0] == name:
                            found = True
                    
                    if not found:
                        listcal.append([name, alpha2000, delta2000])
                        
            
        return(listcal)
    
    
    def queryAlma(self,listname, public = True ):
        "Query the data for the list of ALMA name"
        
        result = []
        
        for name in listname:
            
            region = coordinates.SkyCoord(name[1], name[2], unit='deg')
            alma_data = Alma.query_region(region, 0.003*u.deg, science = False, public =  public)
            
            result.append([name,alma_data])
            
        return(result)
            
            
            
    def selectDeepField(self, data, minTimeBand ={3:1e9,6:1e9,7:1e9}, maxFreqRes = 1000.0, verbose = True):
        """
        From the queryAlma result we filter to select the deep field
        minTimeBand  : dictionary to select the minimum time  per band 
        maxFreqRes: maximum of the frequency resolution (in kHz)
        
        Return the names and projects in a report with the date !!
        """
        
        
        finalReport = []
        
        for item in data:
            name      = item[0][0]
            alpha2000 = item[0][1]
            delta2000 = item[0][2]
            
            nuids = len(item[1])
            
            projects = []
 
            
            totalTime = {3:0., 4:0., 5:0., 6:0., 7:0., 8:0., 9:0., 10:0.}
            
            reportSource  = "#### Source : %s #####\n"%(name)
            reportSource += "#### Coord. 2000: %f  %f \n"%(alpha2000, delta2000)
            reportSource += "\n"
            
            
            
        
            for uids in item[1]:
                selectSG = False
                
                
                code        = uids['Project code']
                source      = uids['Source name']
                band        = int(uids['Band'])
                integration = uids['Integration']
                frequency   = uids['Frequency support']
                obsdate     = uids['Observation date']
                res         = uids['Spatial resolution']
                asdm        = uids['Asdm uid']
                freqRes     = uids['Frequency resolution']
                
                if freqRes < maxFreqRes:
                    selectSG = True
                    totalTime[band] += integration
                
                if verbose and selectSG:
                    reportSource += "## %s  %20s Band:%d obsdate:%s FreqRes:%6.1f Res:%4.2f %s \n"%(code, source, band, obsdate, freqRes, res, asdm)
                
                foundProject = False
                for p in projects:
                    if p == code:
                        foundProject = True
                    

                    
                if not foundProject and selectSG:
                    projects.append(code)
                    
                

            
                
            reportSource += "\n"
            for key in totalTime:
                if totalTime[key] > 0.:
                    reportSource += "Time Band %d : %6.0fs (%3.1fh) \n"%(key, totalTime[key], totalTime[key] / 3600.)
                
            reportSource += "\n Project codes: \n"
            for p in projects:
                reportSource += " %s "%(p)
                
            reportSource += "\n\n"
            
            reportSource += "Total uids: %d \n"%(nuids)
            
            reportSource += "\n"
            
            
            selectSource = True
            for k in minTimeBand:
                if totalTime[k] < minTimeBand[k]:
                    selectSource = False
                    
            if selectSource:
                
                finalReport.append([nuids, item[0], reportSource])
                
        ## sorting according to the number of uids
        finalReportSorted = sorted(finalReport, key=lambda data: data[0])
        
        return(finalReportSorted)
        
        
    def writeReport(self, report, file = "deepfieldRG.txt"):
        "output the final repor from selectDeepField"
        
        fout = open(file,"w")
        
        nsource = 0
        for rep in report:
            nsource += 2
            fout.write(rep[2])
            print(rep[2])
            
        endText = "#################################\n"
        endText = "### Total Number of Sources : %d \n"%(nsource)
        
        fout.write(endText)
        print(endText)
        
        fout.close()
        

##############Main program########################
#### example ....#######################
if __name__=="__main__":
    " main program"   
    
    fileCal = "CalSept2015.list"
    q       = queryCal(fileCal, fluxmin = 3.0)
    data    = q.queryAlma(q.listCal, public = True)
    report  = q.selectDeepField(data, minTimeBand = {3:500., 6:100., 7:1000.}, maxFreqRes = 1000.0, verbose = False)
    q.writeReport(report, file = fileCal+".report")
     
            
    
        
        
    