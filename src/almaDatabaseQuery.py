#!/home/ridlo/anaconda2/bin/python

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

    2017.03.30
        - add option selectPol
        - add flux and acceptedProject in report

    2017.04.07
        - add option: save query result to sql database

    2017.04.09
        - selectdeepfield from sql file 

    2018-02.10
        - start to re-code again, data structure changed by ALMA

    2018-03.30
        - fix some bugs
        - rename the class and functions accoding to PEP8 Python style guide


        
        
"""

__author__="S. Leon @ ALMA, RWW @ ITB"
__version__="0.2.0@2017.03.28"


from astroquery.alma import Alma
from astropy import coordinates
from astropy import units as u

import numpy as np
import pandas as pd

import sys 

reload(sys) #! for Unicode conversion problem
sys.setdefaultencoding('utf8')

import sqlite3 as sql


with open('ALMACAL_object.dat') as f:
    almacal_list = f.read().splitlines()
    

class databaseQuery:
    """Collection of function to query and filter the database of 
    ALMA Calibrator/Archive"""

    def __init__(self):
        pass

    def read_calibratorlist(self, ifile, fluxrange=[0.5, 1.0]):
        '''Read a list of calibrators in CSV format 
        from the Source Catalogue web interface'''
        
        listcal = []
        
        with open(ifile, 'r') as fcal:
            for line in fcal:
                if line[0] != "#":
                    tok       = line.split(",")
                    band      = tok[0].split(" ")[1]
                    name      = tok[13].strip().split("|")[0]
                    alpha2000 = tok[3]
                    delta2000 = tok[5]
                    flux      = float(tok[7]) # flux in a particular Band

                    if (flux >= fluxrange[0] and flux <= fluxrange[1]):
                        found = False
                        for nameYet in listcal:
                            # to remove the duplicate
                            if nameYet[0] == name:
                                found = True

                        if not found:
                            coord = coordinates.SkyCoord(alpha2000 + delta2000, unit=(u.hourangle, u.deg), equinox='J2000') 
                            # convert using astropy
                            listcal.append([name, coord.ra.value, coord.dec.value, flux, band])

        return(listcal)
    
    
    
    def query_object(self, obj, theta = 0.005, science = False, public = True, savedb=False, dbname='calibrators.db'):
        """
        query per object
        obj is an array consist of at least [objname, ra, dec, ...]
        ra and dec are in degrees
        """

        region = coordinates.SkyCoord(obj[1], obj[2], unit='deg')

        alma_data = Alma.query_region(region, theta*u.deg, science = science, public = public)
        # it is in astropy Table format

        df = alma_data.to_pandas() # convert to pandas Dataframe
        
        #! change Numpy Masked-Array to str, so it can be converted to SQL
        df['Band'] = df['Band'].astype("str")

        # save the query result in sql database 
        # with Table's name = calibrator's name
        if savedb:
            conn = sql.connect(dbname)
            conn.text_factory = str
            if not df.dropna().empty:
                df.to_sql(obj[0], conn, if_exists='replace') # if the Table exist in the db, just replace it
                

        return df
    
    
    
    def query_list(self, listobj, science = False, public = True, savedb=False, dbname='calibrators.db'):
        '''
        Query the data for the list of ALMA calibrator (bulk)
        name is an array consist of at least [objname, ra, dec]
        '''
        
        # result = []
        for obj in listobj:
            df = self.query_object(obj, science = science, public = public, savedb=savedb, dbname=dbname)
        #    if not df.dropna().empty:
        #        result.append([obj[0], df])
            
        # return(result)
    
    
    
    def select_object_from_sqldb(self, dbname, \
                                 maxFreqRes = 1000.0, \
                                 array='12m', \
                                 excludeCycle0=True, \
                                 selectPol=False, \
                                 minTimeBand = {3:1e9,6:1e9,7:1e9}, \
                                 verbose = True, silent=True):
        """
        From the sql database we can select some calibrators based on:
        - maxFreqRes    : maximum of the frequency resolution (in kHz)
        - array         : only select project from '12m' array obs.
        - excludeCycle0 : if True, ignore project from Cycle 0
        - selectPol     : if True, only select uid with full pol.
        - minTimeBand   : dictionary to select the minimum time per band 
        
        Return the names and projects in a report with the date.
        """
        connection = sql.connect(dbname)
        connection.text_factory = str # return str, not unicode
        cursor = connection.cursor()

        # make a list of calibrator from table name.
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tablefetch = cursor.fetchall()
        tables = [table[0] for table in tablefetch] # remove tuple

        nsource = 0 # number of selected source
        finalReport = []

        for tab in tables: 
            # table name is the name of source
            # subquery = '['+tab+']'
            subquery = "SELECT * FROM [{0}] WHERE (Array LIKE '%{1}%') AND ([Frequency resolution] < {2})".format(tab, array, maxFreqRes)

            if selectPol:
                subquery += " AND (`Pol products` LIKE '%XY%')"

            if excludeCycle0:
                subquery += " AND (`Project code` NOT LIKE '2011.%')"
            
            
            # and here we execute the query
            # total integration time is calculated after selection for other criterias
            totalTime = {3:0., 4:0., 5:0., 6:0., 7:0., 8:0., 9:0., 10:0.}
            for key in totalTime:
                sqlcmd = "SELECT SUM(Integration) FROM ({0}) WHERE Band LIKE '%{1}%'".format(subquery, key)
                cursor.execute(sqlcmd)
                totalTime[key] = cursor.fetchone()[0]

            for key in totalTime: # remove None
                if totalTime[key] == None:
                    totalTime[key] = 0.0

            # last selection is integration time
            selectSource = True
            sum_time = 0
            for key in minTimeBand:
                if totalTime[key] < minTimeBand[key]:
                    selectSource = False

                sum_time += totalTime[key]

            if sum_time == 0.0:
                selectSource = False

            # if tab in almacal_list:
            #     selectSource = False
            
            if selectSource:
                if not silent:
                    print 'Accepted', tab, totalTime
                
                nsource += 1

                if tab in almacal_list:
                    reportSource = "\n######## Source name: {0} (In ALMACAL) ########\n".format(tab)
                else:
                    reportSource  = "\n######## Source name: {0} ########\n".format(tab)

                reportSource += "\n"
                for key in totalTime:
                    if totalTime[key] > 0.:
                        reportSource += "Time Band %d : %6.0fs (%3.1fh) \n"%(key, totalTime[key], totalTime[key] / 3600.)

                reportSource += "\nList of obs:"
                
                total_number_of_uids = 0
                totalprojects = [] # all project for this source
                for band in totalTime:
                    if totalTime[band] > 0.: # report per Band
                        reportSource += "\n\n### Band: {0} ###".format(band)

                        sqlcmd  = "SELECT `Project code`, `Source name`, RA, Dec, "
                        sqlcmd += "Band, Integration, `Observation date`, `Spatial resolution`, "
                        sqlcmd += "`Asdm uid`, `Frequency resolution`, Array, `Pol products` " + subquery[8:]
                        sqlcmd += " AND Band LIKE '%{0}%'".format(band)
                        sqlcmd += " ORDER BY Integration DESC;" 
                        # Select some columns from criteria
                        # Sort it base on integration time 
                        # in each band
                        
                        cursor.execute(sqlcmd)
                        fetch = cursor.fetchall()

                        projects = [] # project for this Band
                        for uid in fetch:
                            foundProject = False
                            for p in projects:
                                if p == uid[0]:
                                    foundProject = True

                            if not foundProject:
                                projects.append(uid[0])

                            if verbose:
                                reportSource += "\n{0} {1} ra:{2} dec:{3} Band:{4} Int:{5} Obsdate:{6} SpatRes:{7} ASDM:{8} FreqRes:{9} Array:{10} Pol:{11}".format(uid[0], uid[1], str(uid[2])[:9], str(uid[3])[:9], \
                                    uid[4], uid[5], uid[6][:9], str(uid[7])[:5], uid[8], str(uid[9])[:7], uid[10], uid[11])

                        number_of_uid = len(fetch)
                        total_number_of_uids += number_of_uid

                        reportSource += "\n\nTotal accepted uid for Band {0} = {1}".format(band, number_of_uid)
                        reportSource += "\n\nProject codes for Band %d: \n"%(band)

                        for p in projects:
                            reportSource += "%s "%(p)

                            # Calculate total number of projects
                            # one project may (likely) have several uids in different band
                            foundProject = False
                            for q in totalprojects:
                                if p == q:
                                    foundProject = True

                            if not foundProject:
                                totalprojects.append(p)



                total_number_of_projects = len(totalprojects)
                endText  = "\n\n###\nTotal accepted uid for this object = {0}".format(total_number_of_uids)
                endText += "\nTotal accepted project for this object = {0}\n".format(total_number_of_projects)
                endText += "###############################################\n\n\n\n\n\n\n"
                reportSource += endText

                finalReport.append([total_number_of_projects, tab, reportSource])

            else:
                if not silent:
                    print 'Not accepted', tab, totalTime
                pass


        connection.close()
        print "Number of accepted source: ", nsource

        # sorting according to the number of uids
        finalReportSorted = sorted(finalReport, key=lambda data: data[0])
        
        return(finalReportSorted)



    def write_report(self, report, file = "deepfieldRG.txt", silent=True):
        "output the final report from selectDeepField"
        
        fout = open(file,"w")
        
        nsource = 0
        for rep in report:
            nsource += 1
            fout.write(rep[2])
            if not silent:
                print(rep[2])
            
        endText  = "#################################\n"
        endText += "### Total Number of Sources : %d \n"%(nsource)
        
        fout.write(endText)
        print(endText)
        
        fout.close()
        


if __name__=="__main__":
    file_listcal = "callist_20170329.list"

