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

RUN:
 
"""


__author__="S. Leon @ ALMA"
__version__="0.2.3@2017.02.17"


import numpy as np
import pylab as pl
from scipy import signal
from scipy.optimize import curve_fit
import math 

from astroquery.splatalogue import Splatalogue as spla
from astropy import units as u
from astropy import constants as const

import sqlite3 as sql

class analysisLines:
    "class to analyze the lines DB"
    
    def __init__(self,dbname):
        
        self.dbname = dbname
        
        
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
    
    
    def listSources(self):
        "Return the list of sources"
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        c.execute('SELECT DISTINCT source FROM lines')
        sources = c.fetchall()
        
        conn.commit()
        conn.close()
        
        return(sources)       
        
        
    def flagLines(self,edgeFlag = False, telluricFlag = False, duplicateFlag = False):
        "Flag the line with differnt keyword"
        
        
        EDGE_TOL          = 10              ## Flag the line if within EDGE_TOL from the edge
        TELLURIC_TOL      = 2.0e-3          ## tolerance for the lines in different calibrators to be considered telluric (in GHz)
        TELLURIC_MIN_LINE = 2               ## minimum line with same frequencies to consider a telluric line.
        DUPLICATE_TOL     = 2.0e-3          ## tolerance in frequency to consider a duplication
        
        edgestr      = "EDGE"
        telluricstr  = "TELLURIC"
        duplicatestr = "DUPLICATE"
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        

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
            "Check if the same line was in another calibrator"
            
            c.execute('SELECT lineid, source, freq1, freq2, flag FROM lines')
            lines = c.fetchall()
            for li in lines:
                
                id     = li[0]
                sou    = li[1]
                f1     = li[2]
                f2     = li[3]
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
                        print("\n")
                        nFlag += 1
                        cmd = "UPDATE lines SET  flag = '%s' WHERE lineid = %d"%(telluricstr, id)
                        c.execute(cmd)
                    
        
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
        
    def findSpeciesSource(self, sourceName, redshift, DV, flag = False):
        """
        Try to identify the lines for a given source with redshift using the splatalogue DB. 
        If flag = True search only for lines w/o flagging (EDGE, TELURIC, ...)
        DV : uncertainty on the velocity in km/s
        """
        
        lines = self.findSource(sourceName, flag)
        
        if len(lines) == 0:
            print("## Species - No lines found for source %s"%(sourceName))
            return(0)
        
        for li in lines:
            freq1 = li[4] * (1. + redshift)
            freq2 = li[5] * (1. + redshift)
            
            df = freq1 * u.GHz * DV * 1e3 * u.m / u.s /  const.c
            
            print("------")
            print("Frequency redshifted: %f"%(li[4]))
            print("Frequency at rest: %f"%(freq1))
            print df
            
            columns = ('Species','Chemical Name','Resolved QNs','Freq-GHz','Meas Freq-GHz','Log<sub>10</sub> (A<sub>ij</sub>)','E_U (K)','Linelist')
            
            transitions = spla.query_lines(freq1*u.GHz-df,freq2*u.GHz+df)[columns]
            
            transitions.rename_column('Log<sub>10</sub> (A<sub>ij</sub>)','log10(Aij)')
            transitions.rename_column('E_U (K)','EU_K')
            transitions.rename_column('Resolved QNs','QNs')
            transitions.sort('EU_K')
            
        
            transitions.pprint(100)
    
    
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
            
        return(transitions)
      
      
      
    def getInfoDB(self, flag = True):
        """
        To get infos about the DB
            - number of sources
            - number of lines
        """
        
        
        s = self.listSources()
        
        print("## Sources: %d"%(len(s)))
        totalLines = 0
        
        for source in s:
            lines = self.findSource(source[0], flag)
            totalLines += len(lines)
            print("%s : %d lines "%(source[0], len(lines)))
            
        print("\n Total Line: %d"%(totalLines))
        
        