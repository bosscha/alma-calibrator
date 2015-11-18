"""
Class to anlayze the lines DB

2015.11.11:
    - first shot on a Source search
   
   
2015.11.16:
    - return the list of sources
    - add a flagLine method
    
    
RUN:
 
"""


__author__="S. Leon @ ALMA"
__version__="0.1.1@2015.11.16"


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
        
        
        conn = sql.connect(self.dbname)
        c = conn.cursor()
        

        nFlag = 0
        
        if edgeFlag or telluricFlag:
            c.execute('SELECT lineid, maxChannel, chan1, chan2 FROM lines')
            lines = c.fetchall()
            for li in lines:
                
                lineid     = li[0]
                maxChannel = li[1]
                chan1      = li[2]
                chan2      = li[3]
            
                
            
                if chan1 < EDGE_TOL or (maxChannel-chan2) < EDGE_TOL:
                    nFlag += 1
                    print("## Flagged line %d"%(lineid))            
            
        
        conn.commit()
        conn.close()
        
        return(nFlag)             
        
        