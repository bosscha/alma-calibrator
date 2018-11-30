"""
Class to anlayze the lines DB

2015.11.11:
    - first shot on a Source search
   
   
2015.11.16:
    - return the list of sources
    - add a flagLines method
    
2015.11.21: 
    - update  flagLines with Telluric flagging (to be tested, e.g. if same source with different name ...)
    - set the findSpeciesSource method using the splatalogue DB
    
2015.11.23:
    - update telluric line flagging
    
2015.12.07:
    - add flagging for duplicated lines
    - add method clearFlags to clear all flags
    
    
2015.12.10:
    - add astroquery.splatalogue
    
2017.02.04:
    - modify output for the findSpeciesSource

   
2017.02.06:
    - add a finetuning find source if already some hints of velocity are found


2017.02.14:
    - fix minor bug
    
        
2017.02.16:
    - Add method to get information about a DB
    
2017.02.17:
    - add number of lines
    - add redshift in output
    
    
2017.02.18:
    - add information to findSpeciesSource
    
    
2017.02.19:
    - reshaping the flagging to use the sky position and not only the name
    
    
2017.02.20:
    - select source with different coordinates
    - search NED redshift from coordinates
    
    
2017.02.21:
    - adding a method to scan over a list of lines of a source over a list of transition to match several lines at the same time

2017.02.22:
    - minor updates
    - adding findLineSource
    
2017.02.23:
    - adding vel diff to scanningLineSources
    - return result of infoDB
    
2017.02.24:
    - add galactic coordinates in infoDB
    - add a new method for a specific Galactic scanning.
    - add a class for line plotting.
    

2017.02.27:
    - add plotting for one source.
    
    
2017.02.28:
    - important changes in findLinesSource and findLinesCoord (add dataset_id)
    - add wavelet filtering for the display
    
    
2017.03.01:
    - add fitcont for the line plotting.
    
    
2017.03.03:
    - add a scanning method search over a range of velocity on lineid.
    
2017.03.04:
    - add coordinates to the infoDB
    
    
2017.03.06:
    - new class splat to deal directly with the splatalogue.db w/o astroquery (should be faster and no timeout)
    
    
2017.03.07:
    - query to the local splatalogue DB.
    - check for the atmospheric lines..
    
2017.03.14:
    - minor update of getInfoDB
    - add an option to write the lineid selected with the spectra plotting.
    
2017.03.15:
    - update the line plotting method.
    - add a plotEnv to show environment from NED in the list of sources.
    
2017.03.16:
    - add standard deviation to the plot of the spectra.
    
2017.03.18:
    - testing eclipse on mac ....
    
2017.03.18:
    - testing eclipse on mac ....
    
2017.03.20:
    - adding the option to change the splatdb in analysisLines
    - updating splat db. some local splat db must be still missing ...
    
2017.03.28:
    - add resolved_QNs in splat lines.
    
2017.05.25:
    - add an emission/absorption boolean to search for lines
    
2017.05.26:
    - bug fixing.
    
2017.05.29:
    - completion for redshift/velocity  scanning
    
2018.11.23:
    -adding SQL function to query lines
    
RUN:

"""


__author__="S. Leon @ ALMA"
__version__="0.6.9@2017.05.29"


import sys
sys.path.append("../../src/utils")

import numpy as np
import pylab as pl
from scipy import signal
from scipy.optimize import curve_fit
import math 
import os

from astroquery.splatalogue import Splatalogue as spla
from astropy import units as u
from astropy import constants as const
from astropy import coordinates
from astroquery.ned import Ned
from astropy.cosmology import WMAP9 as cosmo

import galenv

import sqlite3 as sql

import matplotlib.pyplot as plt
import wavelet as wav

SOLARSYSTEM = ['Mars','Jupiter','Callisto','Saturn','Titan','Pallas','Ceres','Neptune','Uranus']

DEFAULTSPLAT = "/home/stephane/Science/RadioGalaxy/ALMA/absorptions/splatalogue/splatalogue.db"

class splat():
    "Class to deal with splatalogue.db for faster access"
    
    def __init__(self, splatdb =  DEFAULTSPLAT):
        
        self.splatdb = splatdb        
        self.connectSplat()        
        self.createViewSplat()
        
        
    def connectSplat(self):
        "Connect  to the splatalogue.db"

        self.conn = sql.connect(self.splatdb)
        self.cursor = self.conn.cursor()
        

        
        
    def closeSplat(self):
        "Close the connection to splatalogue.db"
        
        self.conn.commit()
        self.conn.close()
        
        
    def createViewSplat(self):
        "Create a view for the splatalogue"
        
        cmdview = '''
                CREATE temporary VIEW splatLines AS
                select  main.ll_id, main.intintensity , main.sijmu2, main.aij, main.upper_state_energy_K, main.resolved_QNs, main.species_id, main.orderedfreq, 
                species.chemical_name, 
                species.s_name,
                species.atmos,
                species.extragalactic
                from  main
                left join species   on  species.species_id  =  main.species_id  
                                  
                '''
                    
        self.cursor.execute(cmdview) 
    
    
    def querySplat(self,freqlow, freqhigh, energy_max = 1e9, atmosCheck = False):
        """
        Search for lines with the following filters:
            freqlow, freqhigh: frequency range in GHz.
            emax: upper state energy level in K
            atmosCheck : True, check for  atmospheric lines
        """
    
        f1 = min(freqlow,freqhigh)  * 1e3
        f2 = max(freqlow,freqhigh)  * 1e3
        
        if atmosCheck :
            cmd = "SELECT s_name, chemical_name, resolved_QNs, upper_state_energy_K, orderedfreq ,  species_id FROM splatLines WHERE orderedfreq > %f AND orderedfreq < %f AND upper_state_energy_K <= %f AND atmos = 1"%(f1, f2, energy_max)
            
        else:
            cmd = "SELECT s_name, chemical_name, resolved_QNs, upper_state_energy_K, orderedfreq ,  species_id FROM splatLines WHERE orderedfreq > %f AND orderedfreq < %f AND upper_state_energy_K <= %f AND atmos = 0"%(f1, f2, energy_max)
            
        self.cursor.execute(cmd)   
        lines = self.cursor.fetchall()
        
        return(lines)
    
    
##############################################################################################       
class analysisLines:
    "class to analyze the lines DB"
    
    def __init__(self,dbname, splatdb = DEFAULTSPLAT  ):
        
        self.dbname = dbname
        self.splatdb = splatdb
        
      
    def query(self,cmd):
        "Query for the line.db"
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()        
        c.execute(cmd)
          
        res = c.fetchall()
        
        conn.commit()
        conn.close()   
        
        return(res)     
        
        
    def findSource(self,source, flag = False):
        "Find information on a specific source"
                
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        if not flag:
            c.execute('SELECT * FROM lines WHERE source=?', (source,))
        else:
            c.execute('SELECT * FROM lines WHERE source=? AND flag IS NULL', (source,))
        lines = c.fetchall()
        
        conn.commit()
        conn.close()
        
        return(lines)
    
    
    def findLinesSource(self,source, flag = True):
        """
        Same as findLinesCoord with only source name
        """
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        cmdview = '''
                    CREATE temporary VIEW linesky AS
                    select  lines.lineid, lines.dataset_id, lines.source, lines.sn, lines.freq1, lines.freq2, lines.A_fit, lines.mu_fit, lines.sigma_fit, lines.flag,
                    dataset.coordSky1 , 
                    dataset.coordSky2 
                    from  lines
                    left join  dataset   on  dataset.dataid  =  lines.dataset_id                      
                '''
                    
        c.execute(cmdview)
        
        if flag :
            cmd = "SELECT lineid, dataset_id, sn, A_fit, mu_fit, sigma_fit FROM linesky WHERE source = '%s' AND flag IS NULL"%(source)
            
        else :
            cmd = "SELECT lineid, dataset_id, sn, A_fit, mu_fit, sigma_fit FROM linesky WHERE source = '%s' "%(source)

        c.execute(cmd)
        lines = c.fetchall()
        
        conn.commit()
        conn.close()
        
        return(lines)
               
        
    
    def findLinesCoord(self,coord, flag = True):
        """
        return all lines for a given coordinate.
        flag: use the flag or not
        """
        
        TOLPOS = 5e-4
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        cmdview = '''
                    CREATE temporary VIEW linesky AS
                    select  lines.lineid, lines.dataset_id, lines.source, lines.sn, lines.freq1, lines.freq2, lines.A_fit, lines.mu_fit, lines.sigma_fit, lines.flag,
                    dataset.coordSky1 , 
                    dataset.coordSky2 
                    from  lines
                    left join  dataset   on  dataset.dataid  =  lines.dataset_id                      
                '''
                    
        c.execute(cmdview)
        
        if flag :
            cmd = "SELECT lineid, dataset_id, sn, A_fit, mu_fit, sigma_fit FROM linesky WHERE ABS(%f - coordSky1) < %f AND ABS(%f - coordSky2) < %f AND flag IS NULL"%(coord[0], TOLPOS, coord[1], TOLPOS)
        else :
            cmd = "SELECT lineid, dataset_id, sn, A_fit, mu_fit, sigma_fit FROM linesky WHERE ABS(%f - coordSky1) < %f AND ABS(%f - coordSky2) < %f"%(coord[0], TOLPOS, coord[1], TOLPOS)


        c.execute(cmd)
        lines = c.fetchall()
        
        conn.commit()
        conn.close()
        
        return(lines)
        
    
    def listSources(self):
        "Return the list of sources with different coordinates"
        
        TOLPOS = 5e-4
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        c.execute('SELECT DISTINCT calibrator FROM dataset')
        sources = c.fetchall()
        
        listDiffSource = []
        
        for source in sources:
            c.execute("SELECT coordSky1, coordSky2 FROM dataset where calibrator = '%s'"%(source))           
            coord = c.fetchone()
            
            cmd = "SELECT DISTINCT calibrator FROM dataset WHERE ABS(%f - coordSky1) < %f AND ABS(%f - coordSky2) < %f"%(coord[0], TOLPOS, coord[1], TOLPOS)
            c.execute(cmd)
            sourceTemp = c.fetchall()
            
            found = False
            for s in listDiffSource:
                if sourceTemp[0][0] in s:
                    found = True
                    
            
            if not found:
                templist =[]
                for s in sourceTemp:
                    templist.append(s[0])
                listDiffSource.append(templist)       
        
        
        conn.commit()
        conn.close()
        
        return(listDiffSource)       
        
        
    def flagLines(self,edgeFlag = False, telluricFlag = False, duplicateFlag = False):
        "Flag the line with differnt keyword"
        
        
        EDGE_TOL          = 10              ## Flag the line if within EDGE_TOL from the edge
        TELLURIC_TOL      = 1.0e-3          ## tolerance for the lines in different calibrators to be considered telluric (in GHz)
        TELLURIC_MIN_LINE = 3               ## minimum line with same frequencies to consider a telluric line.
        DUPLICATE_TOL     = 2.0e-3          ## tolerance in frequency to consider a duplication
        
        edgestr      = "EDGE"
        telluricstr  = "TELLURIC"
        duplicatestr = "DUPLICATE"
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        ## Create a temporary view withthe lines and  sky coordinates
        ## not yet executed ...
        
        cmdview = '''
                    CREATE temporary VIEW linesky AS
                    select  lines.lineid, lines.source, lines.sn, lines.freq1, lines.freq2, lines.flag,
                    dataset.coordSky1 , 
                    dataset.coordSky2 
                    from  lines
                    left join  dataset   on  dataset.dataid  =  lines.dataset_id                      
                    '''
                    
        ##  c.execute(cmdview)
        print("## View not yet created...")
        ####

        nFlag = 0
        
        if edgeFlag or telluricFlag:
            c.execute('SELECT lineid, maxChannel, chan1, chan2 FROM lines')
            lines = c.fetchall()
            for li in lines:
                
                id         = li[0]
                maxChannel = li[1]
                chan1      = li[2]
                chan2      = li[3]
                
            
                if chan1 < EDGE_TOL or (maxChannel-chan2) < EDGE_TOL:
                    nFlag += 1
                    print("## EDGE - Flagged line %d"%(id))  
                    cmd = "UPDATE lines SET  flag = '%s' WHERE lineid = %d"%(edgestr, id)
                    c.execute(cmd)
                    
            print("\n")
                

        if  telluricFlag:
            "Check if the same line is in another calibrator"
            
            print("## Connect to local splatalogue for telluric lines ...")
            spl = splat(splatdb = self.splatdb)
            
            
            c.execute('SELECT lineid, source, freq1, freq2, flag FROM lines')
            lines = c.fetchall()
            for li in lines:
                
                id     = li[0]
                sou    = li[1]
                f1     = min(li[2],li[3])
                f2     = max(li[2],li[3])
                flag   = li[4]
                
                if flag != "EDGE":
                    cmd = "SELECT lineid, source, freq1, freq2, flag FROM lines  WHERE lineid != %d AND  ABS(freq1 - %f) < %f AND ABS(freq2 - %f) < %f AND source != '%s' AND (flag IS NULL or flag != '%s') "%(id, f1, TELLURIC_TOL, f2, TELLURIC_TOL, sou, "EDGE") 
                                                                                                
                    c.execute(cmd)
                    
                    ## note that we do not check yet if the same source is repeated...
                    
                    linesTell =  c.fetchall()
                    
                    if len(linesTell) > TELLURIC_MIN_LINE -1:
                        print("## TELLURIC - Flagged line %d"%(id))
                        print("## TELLURIC - Source : %s"%(sou))
                        print("## TELLURIC - Frequencies : %f - %f GHz"%(f1,f2))
                        print("## TELLURIC - other sources (frequencies) :")
                        for item in linesTell:
                            idtell  = item[0]
                            soutell = item[1]
                            f1tell  = item[2]
                            f2tell  = item[3]
                            
                            print("#### %s (%f,%f)"%(soutell, f1tell, f2tell))
                            
                            
                        ## check if telluric is tagged in splatalogue
                        telluric = spl.querySplat(f1, f2, energy_max = 1e9, atmosCheck = True)
                        if len(telluric) > 0:
                            print("## TELLURIC - SPLATALOGUE:")
                            print telluric                               
                                                             
                        print("\n")
                        
                        nFlag += 1
                        cmd = "UPDATE lines SET  flag = '%s' WHERE lineid = %d"%(telluricstr, id)
                        c.execute(cmd)
                    
        spl.closeSplat()
        
        
        if duplicateFlag:
            "Check for duplicated lines"
            
            c.execute('SELECT lineid, source, sn,  freq1, freq2, flag FROM lines WHERE flag IS  NULL')
            lines = c.fetchall()
            for li in lines:
                
                id     = li[0]
                sou    = li[1]
                sn     = li[2]
                f1     = li[3]
                f2     = li[4]
                flag   = li[5]
                
                cmd = "SELECT lineid, source, sn, freq1, freq2, flag FROM lines  WHERE source = '%s' AND ABS(freq1 - %f) < %f AND ABS(freq2 - %f) < %f AND flag IS NULL AND lineid != %d"%(sou,f1,DUPLICATE_TOL, f2,DUPLICATE_TOL,id)                                   
                c.execute(cmd)
        
                linesdupli = c.fetchall()
                
                
                if len(linesdupli) > 0:
                    
                    for item in linesdupli:
                        iddupli = item[0]
                        sndupli = item[2]
                        
                        if sn > sndupli:
                            idflag = iddupli
                        else:
                            idflag = id
                        
                           
                        print("## DUPLICATE - Flagged line %d"%(idflag))
                        print("## DUPLICATE - Source : %s"%(sou))
                        print("\n")
                        nFlag += 1 
                        cmd = "UPDATE lines SET  flag = '%s' WHERE lineid = %d"%(duplicatestr, idflag)
                        c.execute(cmd)
                        
                        
        conn.commit()
        conn.close()
        
        return(nFlag)             
        
        
    def clearFlags(self):
        "Clear all flags"
        
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        c.execute("UPDATE lines SET flag = NULL")
        
        conn.commit()
        conn.close()
        
        return(0)
        
    def findSpeciesSource(self, sourceName, redshift, DV, flag = False, outputline = 50):
        """
        Try to identify the lines for a given source with redshift using the splatalogue DB. 
        If flag = True search only for lines w/o flagging (EDGE, TELURIC, ...)
        DV : uncertainty on the velocity in km/s
        """
        
        lines = self.findSource(sourceName, flag)
        
        if len(lines) == 0:
            print("## Species - No lines found for source %s"%(sourceName))
            return(0)
        
        resLines = []
        for li in lines:
            freq1 = li[4] * (1. + redshift)
            freq2 = li[5] * (1. + redshift)
            
            df = freq1 * u.GHz * DV * 1e3 * u.m / u.s /  const.c
            
            print("------")
            print("Source: %s"%(sourceName))
            print("Redshift: %f"%(redshift))
            print("Frequency redshifted: %f"%(li[4]))
            print("Frequency at rest: %f"%(freq1))
            print("Velocity offset: +/- %f km/s"%(DV))
            
            columns = ('Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)','Linelist')
            
            transitions = spla.query_lines(freq1*u.GHz-df,freq2*u.GHz+df)[columns]
            
            transitions.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
            transitions.rename_column('E_U (K)','EU_K')
            transitions.rename_column('Resolved QNs','QNs')
            transitions.sort('EU_K')
            
        
            transitions.pprint(outputline)
            resLines.append(transitions)
            
        return(transitions)
      
    
    def findSpeciesFineTuningSource(self, sourceName, redshift, DV1, DV2, flag = False, outputline = 50):
        """
        Try to identify the lines for a given source with redshift using the splatalogue DB. 
        If flag = True search only for lines w/o flagging (EDGE, TELURIC, ...)
        DV1,Dv2 :  range in velocity
        """  
        
        lines = self.findSource(sourceName, flag)
        
        if len(lines) == 0:
            print("## Species - No lines found for source %s"%(sourceName))
            return(0)
        
        resLines = []
        for li in lines:
            freq1 = li[4] * (1. + redshift)
            freq2 = li[5] * (1. + redshift)
            
            df1 = freq1 * u.GHz * DV1 * 1e3 * u.m / u.s /  const.c
            df2 = freq2 * u.GHz * DV2 * 1e3 * u.m / u.s /  const.c
            
            print("------")
            print("Source: %s"%(sourceName))
            print("Redshift: %f"%(redshift))
            print("Frequency redshifted: %f"%(li[4]))
            print("Frequency at rest: %f"%(freq1))
            print("Frequency offset:")
            print df1, df2
            
            
            columns = ('Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)','Linelist')
            
            transitions = spla.query_lines(freq1*u.GHz+df1,freq2*u.GHz+df2)[columns]
            
            transitions.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
            transitions.rename_column('E_U (K)','EU_K')
            transitions.rename_column('Resolved QNs','QNs')
            transitions.sort('EU_K')
            
        
            transitions.pprint(outputline)      
            resLines.append(transitions)
            
        return(resLines)
      
    
    def findSpeciesFrequencyVelocity(self, freq, vel, dvel):
        " Find species at an observed frequency and for a given vel, dvel at low redshift ..."
        
        freq1 = freq / (1 + 1e3 * (vel - dvel) / const.c.value)
        freq2 = freq / (1 + 1e3 * (vel + dvel) / const.c.value)
    
        transitions = spla.query_lines(freq1*u.GHz,freq2*u.GHz)[columns]
        transitions.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
        transitions.rename_column('E_U (K)','EU_K')
        transitions.rename_column('Resolved QNs','QNs')
        transitions.sort('EU_K')
            
        
        transitions.pprint(10)  
        
    
    
    def getCoordSource(self,source):
        "Get coordinates for a source "
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        c.execute("SELECT coordSky1, coordSky2 FROM dataset where calibrator = '%s'"%(source))           
        coord = c.fetchone()
        
        conn.commit()
        conn.close()
        
        return(coord)
    
    
    def getMS(self,coord1, coord2):
        "Get MS for a position"
                      
        TOLPOS = 5e-4
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        cmd = "SELECT DISTINCT msfile FROM dataset WHERE ABS(%f - coordSky1) < %f AND ABS(%f - coordSky2) < %f"%(coord1, TOLPOS, coord2, TOLPOS)
        c.execute(cmd)
        data = c.fetchall()

        return(data)
        
      
    def getInfoDB(self, flag = True , verbose = True):
        """
        To get infos about the DB
            - number of sources
            - number of lines
            - galactic coordinate
        """
        
        
        s = sorted(self.listSources())
        
        print("## Sources: %d"%(len(s)))
        print("##")
        totalLines = 0
        resultList = []
        
        for sourceDiff in s:
            sText =  ""
            totalLinesOne = 0
            for source in sourceDiff:
                lines = self.findSource(source, flag)
                sText += "%s  "%(source)
                totalLinesOne += len(lines)
            if verbose : print("## %s : %d lines "%(sText, totalLinesOne))
            
            ## redshift
            coord = self.getCoordSource(sourceDiff[0])
            
            if verbose : print("## Coordinates (J2000): %f %f"%(coord[0],coord[1]))
            ## Get Galactic coordinates:
                
            coordGal = self.getGalacticCoord(coord[0],coord[1])
            if verbose : print("## Galactic Coordinates (l,b): (%f , %f) "%(coordGal.galactic.l.value,coordGal.galactic.b.value))
            
            try:
                redshift = self.findNEDredshift(coord[0],coord[1])
                if verbose:
                    print("## Redshift:")
                    print redshift
                
                
                ## Number of MS
                msfile = self.getMS(coord[0],coord[1])
                
                if verbose:
                    print("##")
                    print("## Number of MS: %d"%(len(msfile)))
                resultList.append([sourceDiff[0], redshift])
            
            except:
                if verbose : print("## No redshift")
                ## put -99 for redshift if unknown
                resultList.append([sourceDiff[0], -99])

            if verbose : print("")
            totalLines += totalLinesOne
        if verbose : print("\n Total Line: %d"%(totalLines))
        
        return(resultList)


   
    def findNEDredshift(self, coord1, coord2):
        """
        Try to find redshift for ccoordinates
        """
        
        RADIUSCONE = 5e-4
        
        co = coordinates.SkyCoord(ra=coord1, dec=coord2, unit=(u.deg, u.deg), frame='fk5')
        result_table = Ned.query_region(co, radius= RADIUSCONE * u.deg, equinox='J2000.0')
        result_table.sort('Distance (arcmin)')
        redshift = result_table['Redshift'].quantity
        
        return(redshift)
        
        
        
    def getGalacticCoord(self, coord1, coord2):
        """
        Transform equatorial coordinates to Galactic coordinates
        
        """
    
        co = coordinates.SkyCoord(ra=coord1, dec=coord2, unit=(u.deg, u.deg), frame='fk5')
        
        cogal = co.galactic
        
        return(cogal)
    
    
    
    def scanningLinesSources(self, source, lines, templateLine, maxRedshift, nummatchmin = 2 ):
        """
        Scan over the lines of a source using the templateLine from redshift 0-dz up to maxRedshift+dz
        
        templateLine format : array of [string transition, frequency (GHz)]
        
        lineid = if len > 0 selec on the lineid only for the scanning
        
        """
    
        DZ        = 1e-5
        FREQTOL     = 5e-4      ## in GHz
        
        minRedshift = -0.05
        
        NSCAN       = int((maxRedshift-minRedshift) / DZ)
        
         
        print("## Scanning lines ...")
        print("## %d scans"%(NSCAN))
        
        lineMatchingTotal = []
    
        for i in range(NSCAN):
            z = minRedshift + i * DZ

            matchline = []
            
            for line in lines:
                diffmin = 1e9
                transmatch = []
                
                for  trans in templateLine:
                    diffFreq = abs( (trans[1] / (1.+z)) - line[3])
                    
                    if diffFreq <  diffmin :
                        diffmin  = diffFreq
                        transmatch = trans
                    
                if diffmin < FREQTOL :
                        veldiff = const.c * 1e-3 * diffmin * (1. +z)*(1.+z) / line[3]                        
                        matchline.append([transmatch, veldiff.value,line[1]  ])

                
            nmatch =   len(matchline) 
            
            if  nmatch >= nummatchmin :
                dist = cosmo.comoving_distance(z)
                vel  = z * const.c
                print("## %d matches for redshift: %f"%(nmatch, z))
                print("## Velocity (km/s) : %f"%(vel.value / 1000.))
                print("## Distance (Mpc)  : %f"%(dist.value))
                print("## Line, vel. diff.  wrt. redshift (km/s), S/N:")
                for match in matchline:
                    print match
                print("## \n")
                
                lineMatchingTotal.append([matchline,z])
                
        return(lineMatchingTotal)   



    def scanningSplatVelocitySourceLineid(self, lineid, velmin, velmax, dv, nrao = True, emax = 50., local = True, absorption = True, emission = True):
        """
            Scan over a list of lineid (for a similar source but not checked..) using splatalogue. It scans over velocity range and should 
            be mainly Galactic.
            Input:
                lineid : array of lineid.
                velmin, velmax: range of velocity
                dv: velocity step (km/s)
                nrao: True if only nrao recommended
                emax: maximum energy upper level
                local: use splat class and local DB
                
        """
        
        ## Get the lines... 
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        lines =[]
        
        for lineidit in lineid:
            if emission and absorption:
                cmd = "SELECT lineid, dataset_id, source, sn, A_fit, mu_fit, sigma_fit FROM lines WHERE lineid = %d"%(lineidit)
            if not absorption and emission:
                cmd = "SELECT lineid, dataset_id, source, sn, A_fit, mu_fit, sigma_fit FROM lines WHERE lineid = %d AND A_fit > 0."%(lineidit)
            if not emission and absorption:
                cmd = "SELECT lineid, dataset_id, source, sn, A_fit, mu_fit, sigma_fit FROM lines WHERE lineid = %d AND A_fit < 0."%(lineidit)
            
            c.execute(cmd)
            line = c.fetchall()
            
            if len(line) > 0:
                lines.append(line)
        
        conn.commit()
        conn.close()
        
        
        print("## Number of lines: %d"%(len(lines)))
        
        
        ### Scan over the velocity range ..
        print("## Galactic scanning ...")
        
        if local:
            print("## Connect to local splatalogue ...")
            spl = splat(splatdb = self.splatdb)
            
        else:
            columns = ('Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)','Linelist')
  
        
        nv = int((velmax -velmin) / dv)       
        vel = np.arange(nv) * dv + velmin
        
        linematches = []
        
        for v in vel:
            lineVel = []
            lineVel.append(v)            
            
            completion = 100. * (v-velmin) / (velmax - velmin)
            print("## Velocity: %5.5f   (%3.1f%% done)"%(v, completion))
            
            nline = 0
            
            for l in lines:
                freq = l[0][5]
                source = l[0][2]
                
                
                freq1 = freq / (1 + 1e3 * (v + dv /2.) / const.c.value)
                freq2 = freq / (1 + 1e3 * (v - dv /2.) / const.c.value)
                
                
                if local :                    
                    transitions = spl.querySplat(freq1, freq2, energy_max = emax)
                    
                else :         
                    transitions = spla.query_lines(freq1*u.GHz,freq2*u.GHz,only_NRAO_recommended = nrao, energy_max= emax, energy_type='eu_k' )[columns]               
                    transitions.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
                    transitions.rename_column('E_U (K)','EU_K')
                    transitions.rename_column('Resolved QNs','QNs')
                    transitions.sort('EU_K')
                
                
                # transitions.pprint(5)

                lineVel.append(transitions) 
                   
                if len(transitions) > 0:
                    nline += 1 
    
            print("### Lines found: %d"%(nline))

 
            
            linematches.append(lineVel)
            
        if local:
            spl.closeSplat()

        
        return(linematches)
        
        

    def scanningSplatRedshiftSourceLineid(self, lineid, zmin, zmax, dz, nrao = True, emax = 50., local = True,  absorption = True, emission = True):
        """
            Scan over a list of lineid (for a similar source but not checked..) using splatalogue. It scans over redshift range and should 
            be mainly Galactic.
            Input:
                lineid : array of lineid.
                zmin, zmax: range of redshift
                dz: redshift step
                nrao: True if only nrao recommended
                emax: maximum energy upper level
                local: use splat class and local DB
                
        """
        
        
        ## Get the lines... 
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        lines =[]
        
        for lineidit in lineid:
            if emission and absorption:
                cmd = "SELECT lineid, dataset_id, source, sn, A_fit, mu_fit, sigma_fit FROM lines WHERE lineid = %d"%(lineidit)
            if not absorption and emission:
                cmd = "SELECT lineid, dataset_id, source, sn, A_fit, mu_fit, sigma_fit FROM lines WHERE lineid = %d AND A_fit > 0."%(lineidit)
            if not emission and absorption:
                cmd = "SELECT lineid, dataset_id, source, sn, A_fit, mu_fit, sigma_fit FROM lines WHERE lineid = %d AND A_fit < 0."%(lineidit)
            
            c.execute(cmd)
            line = c.fetchall()
            
            
            if len(line) > 0:
                lines.append(line)
            
        print("## Number of lines: %d"%(len(lines)))
        
        conn.commit()
        conn.close()
        
        
        ### Scan over the redshift range ..
        print("## Redshift scanning ...")
        
        if local:
            print("## Connect to local splatalogue ...")
            
            spl = splat(splatdb = self.splatdb)
            
        else:
            columns = ('Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)','Linelist')
 
        nz = int((zmax -zmin) / dz)       
        redshift = np.arange(nz) * dz + zmin
        
        linematches = []
        
        for z in redshift:
            nline = 0
            lineZ = []
            lineZ.append(z)
            
            completion = 100. * (z - zmin) / (zmax - zmin)
            print("## Redshift: %5.5f   (%3.1f%% done)"%(z, completion))
            
            for l in lines:
                freq = l[0][5]
                source = l[0][2]
                
                freq1 = freq * (1. + z - dz / 2.)
                freq2 = freq * (1. + z +  dz / 2.)
                
                #print freq1, freq2
                
                if local:
                    transitions = spl.querySplat(freq1, freq2, energy_max = emax)
                else:
                    transitions = spla.query_lines(freq1*u.GHz,freq2*u.GHz, only_NRAO_recommended = nrao, energy_max= emax, energy_type='eu_k')[columns]               
                    transitions.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
                    transitions.rename_column('E_U (K)','EU_K')
                    transitions.rename_column('Resolved QNs','QNs')
                    transitions.sort('EU_K')
                
                    # transitions.pprint(5)

                lineZ.append(transitions) 
                    
                
                if len(transitions) > 0:
                    nline += 1 
    
            print("### Lines found: %d"%(nline))

 
            
            linematches.append(lineZ)
            
        if local:
            spl.closeSplat()
                
        return(linematches)   
         
############################################################################################################       
        

class plotLines:
    "Class to plot the lines"
    
    
    def __init__(self,dbname, dirdata, dirplot):
        
        self.dbname  = dbname        
        self.dirdata = dirdata
        self.dirplot = dirplot
    
    
    
    def extractData(self,dataFile, casa = True):
        "Extract the data from the dataFile, special case for CASA file"
        
   
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
            
    
            # print nSpwData
        
            for i in range(nSpwData) :
                freq[i] = ftemp[i]
                amp[i]  = atemp[i] / nDat[i]  
            
            
            
        return(freq, amp)

             
    
    def getDataSource(self, source , coordCheck = True):
        """ 
        Get all the spectra for a gven source with the dataid.
        If coordCheck is True it uses the coordinates to get the data.
        """
        
        TOLPOS = 5e-4
        
        
        if coordCheck:
            al = analysisLines(self.dbname)
            coord = al.getCoordSource(source)
            
        
        ## loop over the data files
               
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        if not coordCheck:
            cmd = "SELECT dataid, filedata  FROM dataset where calibrator = '%s'"%(source)
            c.execute(cmd)          
            data = coord = c.fetchall()
         
        else:   
            cmd = "SELECT dataid, filedata FROM dataset WHERE ABS(%f - coordSky1) < %f AND ABS(%f - coordSky2) < %f"%(coord[0], TOLPOS, coord[1], TOLPOS)
            c.execute(cmd)
            data = c.fetchall()
  

        if len(data) == 0 :
            print("## No data...")
            return([])
        
        
        dataList = []
        
        ## change to data directory
        ##
        
        curDir = os.getcwd()
        os.chdir(self.dirdata)
        
        for filedata in data:
            # print filedata[0]
            # print filedata[1]
            
            try:
                freq, amp = self.extractData(filedata[1])            
                dataList.append([source, filedata[0], freq, amp])
                
            except:
                print("## Problem to extract data in %s"%(filedata[1]))
            
        # back to the initial directory
        os.chdir(curDir)
        
        print("## Data loaded in memory...")
        
        return(dataList)
    
    
    def plotDataSet(self, data, dataid, freqmin, freqmax, ampmin, ampmax):
        """
        plot for a single dataid 
         - data extracted for the source
         - freqmin, freqmax, ampmin, ampmax : range for the plot
         
         
        """
        
        ## only one result !!
        found = False
        
        for d in data:
            if d[1] == dataid:
                res = d
                found = True
                
                
        if found:
            fig, ax = plt.subplots()
            
            ax.plot(res[2],res[3])
            plt.show()
            
        
    
    
    def   plotLineSource(self, source, data, flag = True, coordCheck = True, wavFiltering = True, wavscale = 4, wavsn = 2.7, fitcont = 0, lineidShow = True, plotWidthFac = 35.):
        """
        plot the lines for a source with the following inputs:
            - source : name of the source
            - data  : data for the source
            - flag : True select only unflagged lines
            - coordCheck : True if check the coord instead of the source
            - wavFiltering : True if perform wavelet filtering
            - if fitcont > 0 a polynom is fitted and line/cont computed
            - if lineidShow is True will display the array of "s" selected
            - plotWidthFac : the range of the plot in term of linewidth.
            
            
        """
        al = analysisLines(self.dbname)
        
        lineidSelected = []
        
        print("## Source : %s"%(source))
        
        ## wavelet parameter.. freq and amp are numpy array
        if wavFiltering:
            wt = wav.wt(verbose = False)
            print("## Wavelet filtering ...")
        
        ## Get coord
        if coordCheck:
            coord = al.getCoordSource(source)
            lines = al.findLinesCoord(coord, flag)
            
        else:
            lines = al.findLinesSource(source)
        

        ## interactive mode on
        plt.ion()
        
        if lineidShow:
            print("##")            
            print("## type ''s'' to select the line.")
            print("##")
        
        for line in lines:
            datasetid = line[1]
            lineid    = line[0]
            
            for d in data:
                if d[1] == datasetid:
                    print("## Data found (plotLineSource)")
                    print line
                    
                    freq = d[2]
                    amp  = d[3]
                    
                    ## if fitcont > 0 fits the continuum and produces line/cont
                    if fitcont > 0:
                        print("## Fitting the continuum ...")
                        z = np.polyfit(freq, amp, fitcont)
                        p = np.poly1d(z)
            
                        ampTemp = amp / p(freq)
                        amp = ampTemp
                    
                    print("## Dataid : %d"%(datasetid))
                    print("## Lineid : %d"%(lineid))
                    ## Filtering with wavelet ...
                    
                    if wavFiltering:
                            noise = np.std(amp)
                            mean  = np.mean(amp) 
                            wspw = wt.atrous1d(amp, wavscale)
                            wspwFilt = wt.filtering1d(wspw, wavsn, waveletNoise = True, spectralNoise = noise)
                            spwRestore = wt.restore1d(wspwFilt, 0, wavscale)
                            stdFiltered = np.std(spwRestore)
                            print("## Standard Deviation filtered: %f"%(stdFiltered))
                            
                    else :
                        spwRestore = amp
                        stdNoFiltered = np.std(spwRestore)
                        print("## Standard Deviation : %f"%(stdNoFiltered))
                    
                    ## plotting 
                    xmin = line[4] - line[5] * plotWidthFac
                    xmax = line[4] + line[5] * plotWidthFac
                    
                    plt.clf()
                    figfile = "toto.png"
                    
                    # plt.title(figfile.split('/')[-1])
                    plt.xlabel("Frequency (GHz)")
                    plt.ylabel("Amplitude")
                    plt.xlim(xmin,xmax)
                    # plt.ylim([ymin,ymax])
                    plt.plot(freq, spwRestore)
                    #fig.savefig(figfile)
                    
                    c = raw_input()
                    
                    if c == 's':
                        lineidSelected.append(lineid)
                        
        
        if lineidShow:
            print("## Selected lineid =")
            print lineidSelected
                  
        return(lineidSelected)
    
    
    def plotEnvSources(self, theta= 0.1 ,  showFig = True):
        """
        Show the environment for the sources in the DB
        theta: cone angle in degrees
        """
        
        
        ge = galenv.Galenv()
        
        al = analysisLines(self.dbname)
        
        s = sorted(al.listSources())
        
        for sourceDiff in s:
            
            print("##")
            print("## Source: %s"%(sourceDiff[0]))
            print("##")
            
            coord = al.getCoordSource(sourceDiff[0])
            co = coordinates.SkyCoord(ra = coord[0], dec = coord[1], unit=(u.deg, u.deg), equinox= "J2000")
        
            ge.cone_bypos(co, theta, show = showFig, savefig=False, imgname="plot.png")    
        