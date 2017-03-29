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
        - fixing a bug
        - adapted to async by RWW
        - add an option
        - still a bug with verbose = True ...

    2017.03.29
        - fixing verbose = True
        - add safety if query is fail (just skip) | if fail, try to change the radius in query
        - remove the temporary xml-file
    
        
        
"""

__author__="S. Leon @ ALMA, RWW @ ITB"
__version__="0.2.0@2017.03.28"




from astroquery.alma import Alma
from astropy import coordinates
from astropy import units as u
import time
import os
import numpy as np
import pandas as pd
from xml.etree import ElementTree as ET
from bs4 import BeautifulSoup


class queryCal:
    
    def __init__(self,fileCal, fluxrange, readFile = True):
        
        if readFile:
            self.listCal = self.readCal(fileCal, fluxrange)
        else:
            self.listCal = []



    def readCal(self, ifile, fluxrange = [0.5, 1.0]):
        "Read a list of calibrators in CSV format from the Source Catalogue web interface"

        listcal = []

        with open(ifile, 'r') as fcal:
            for line in fcal:
                if line[0] != "#":
                    tok       = line.split(",")
                    band      = tok[0].split(" ")[1]
                    flux      = float(tok[7])
                    name      = tok[13].strip().split("|")[0]
                    alpha2000 = tok[3]
                    delta2000 = tok[5]

                    if (flux > fluxrange[0] and flux < fluxrange[1]):
                        found = False
                        for nameYet in listcal:
                            if nameYet[0] == name:
                                found = True

                        if not found:
                            coord = coordinates.SkyCoord(alpha2000 + delta2000, unit=(u.hourangle, u.deg), equinox='J2000') # converted using astropy
                            listcal.append([name, coord.ra.value, coord.dec.value, flux])

        return(listcal)



    def readCalold(self,file, fluxrange = [0.5, 1.0]):
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
                
                if (flux > fluxrange[0] and flux < fluxrange[1]):
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
            alma_data = Alma.query_region_async(region, 0.004*u.deg, science = False, public =  public)
            
            # alma_data is in unicode
            xmldata = BeautifulSoup(alma_data.text, "lxml")

            if len(xmldata.contents) < 2:
                print name, " query FAIL!"
            else:
                with open('xmldata.xml', 'w') as ofile: # write in a file, 
                    ofile.write(str(xmldata.contents[1]))

                tree = ET.parse("xmldata.xml") # read the file again. Find a better way.
                root = tree.getroot()

                projects = []
                for i, proj in enumerate(root[0][0][0][2][36][0]): # many child
                    projects.append([])
                    for data in proj:
                        if data.text is None:
                            projects[i].append('None')
                        else:
                            projects[i].append(data.text)

                columns=['Project code', 'Source name', 'RA', 'Dec', 'Galactic longitude', 'Galactic latitude', \
                         'Band', 'Spatial resolution', 'Frequency resolution', 'Array', 'Mosaic', 'Integration', \
                         'Release date', 'Frequency support', 'Velocity resolution', 'Pol products', \
                         'Observation date', 'PI name', 'SB name', 'Proposal authors', 'Line sensitivity (10 km/s)', \
                         'Continuum sensitivity', 'PWV', 'Group ous id', 'Member ous id', 'Asdm uid', 'Project title', \
                         'Project type', 'Scan intent', 'Field of view', 'Largest angular scale', 'QA2 Status',\
                         'Pub', 'Science keyword', 'Scientific category', 'ASA_PROJECT_CODE']

                #convert python list to Pandas DataFrame
                df = pd.DataFrame(projects, columns=columns)

                result.append([name, df])

        if os.path.exists('xmldata.xml'): # remove the temporary xml file
            os.remove('xmldata.xml')
            
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
            
            projects = []
            
            totalTime = {3:0., 4:0., 5:0., 6:0., 7:0., 8:0., 9:0., 10:0.}
            
            reportSource  = "#### Source : %s #####\n"%(name)
            reportSource += "#### Coord. 2000: %f  %f \n"%(alpha2000, delta2000)
            reportSource += "\n"
            
            uids = item[1]
            nuids = len(item[1])

            code        = uids['Project code']
            source      = uids['Source name']
            band        = uids['Band']
            integration = uids['Integration']
            frequency   = uids['Frequency support']
            obsdate     = uids['Observation date']
            res         = uids['Spatial resolution']
            asdm        = uids['Asdm uid']
            freqRes     = uids['Frequency resolution']
            pol         = uids['Pol products']
            

            # selection
            for i in range(nuids):
                selectSG = False
                if float(freqRes[i]) < maxFreqRes:
                    selectSG = True
                    totalTime[int(band[i])] += float(integration[i])

                if verbose and selectSG:
                    reportSource += "## %s  %20s  Band:%s  obsdate:%s  FreqRes:%s  Res:%s  Pol:%s    asdm:%s\n"%(code[i], source[i], band[i], obsdate[i], freqRes[i], res[i], pol[i], asdm[i])
                
                foundProject = False
                for p in projects:
                    if p == code[i]:
                        foundProject = True
                    
                if not foundProject and selectSG:
                    projects.append(code[i])

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
                
        # sorting according to the number of uids
        finalReportSorted = sorted(finalReport, key=lambda data: data[0])
        
        return(finalReportSorted)

        
        
    def writeReport(self, report, file = "deepfieldRG.txt"):
        "output the final report from selectDeepField"
        
        fout = open(file,"w")
        
        nsource = 0
        for rep in report:
            nsource += 1
            fout.write(rep[2])
            print(rep[2])
            
        endText  = "#################################\n"
        endText += "### Total Number of Sources : %d \n"%(nsource)
        
        fout.write(endText)
        print(endText)
        
        fout.close()
        

##############Main program########################
#### example ....#######################
if __name__=="__main__":
    " main program"   
    
    fileCal = "callist_20170329.list"
    q       = queryCal(fileCal, fluxrange = [0.95, 1.0])
    data    = q.queryAlma(q.listCal, public = True)
    report  = q.selectDeepField(data, minTimeBand = {3:500., 6:100., 7:1000.}, maxFreqRes = 1000.0, verbose = False)
    q.writeReport(report, file = fileCal+".report")