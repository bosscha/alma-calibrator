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
    
    
RUN:
 
"""


__author__="S. Leon @ ALMA"
__version__="0.1.3@2015.11.21"


import numpy as np
import pylab as pl
from scipy import signal
from scipy.optimize import curve_fit
import math 


import sqlite3 as sql

class analysisLines:
    "class to analyze the lines DB"
    
    def __init__(self,dbname):
        
        self.dbname = dbname
        
        
    def findSource(self,source):
        "Find information on a specific source"
                
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        
        c.execute('SELECT * FROM lines WHERE source=?', (source,))
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
        
        
    def flagLines(self,edgeFlag = False, telluricFlag = False):
        "Flag the line with differnt keyword"
        
        
        EDGE_TOL = 10           ## Flag the line if within EDGE_TOL from the edge
        TELLURIC_TOL = 2.0      ## tolerance for the lines in different calibrators to be considered telluric (in MHz)
        
        edgestr = "EDGE"
        
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
                    cmd = "UPDATE lines SET  status = '%s' WHERE lineid = %d"%(edgestr, id)
                    c.execute(cmd)
                    
        
        if  telluricFlag:
            "Check if the same line was in another calibrator"
            
            c.execute('SELECT lineid, source, freq1, freq2, status FROM lines')
            lines = c.fetchall()
            for li in lines:
                
                id     = li[0]
                sou    = li[1]
                f1      = li[2]
                f2      = li[3]
                stat     = li[4]
                
                if stat != "EDGE":
                    c.execute('SELECT lineid, source, freq1, freq2, status FROM lines  WHERE lineid != ? AND  \
                        ABS(freq1 - ?) < ? AND ABS(freq1 - ?) < ? AND source != ? AND status != ?', (id, f1, f2, sou,"EDGE" ))
                    
                    linesTell =  c.fetchall()
                    
                    if len(linesTell) > 0:
                        print("## TELLURIC - Flagged line %d"%(id))
                        print("## TELLURIC - other sources (frequencies) :")
                        for item in linesTell:
                            idtell  = item[0]
                            soutell = item[1]
                            f1tell  = item[2]
                            f2tell  = item[3]
                            
                            print("#### %s (%f,%f"%(soutel, f1tell, f2tell))
                                                    
        
        conn.commit()
        conn.close()
        
        return(nFlag)             
        
        
    def findSpeciesSource(self, sourceName, redshift, flag = False):
        """
        Try to identify the lines for a given source with redshift using the splatalogue DB. 
        If flag = True search only for lines w/o flagging (EDGE, TELURIC, ...)
        """
        
        lines = self.findSources(sourceName)
        
        if len(lines) == 0:
            print("## Species - No lines found for source %s"%(sourceName))
            return(0)
        
        for li in lines:
            pass
        