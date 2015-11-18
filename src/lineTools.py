"""
Class to anlayze the lines DB

2015.11.11
    - first shot on a Source search
   
   
RUN:
 
"""


__author__="S. Leon @ ALMA"
__version__="0.1.0@2015.11.11"


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