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

    2017.03.30
        - add option selectPol
        - add flux and acceptedProject in report

    2017.04.07
        - add option: save query result to sql database

    2017.04.09
        - selectdeepfield from sql file 
        
        
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

import sqlite3 as sql


columns=['Project_code', 'Source_name', 'RA', 'Dec', 'Galactic_longitude', 'Galactic_latitude', \
         'Band', 'Spatial_resolution', 'Frequency_resolution', 'Array', 'Mosaic', 'Integration', \
         'Release_date', 'Frequency_support', 'Velocity_resolution', 'Pol_products', \
         'Observation_date', 'PI_name', 'SB_name', 'Proposal_authors', 'Line_sensitivity_10km_s', \
         'Continuum_sensitivity', 'PWV', 'Group_ous_id', 'Member_ous_id', 'Asdm_uid', 'Project_title', \
         'Project_type', 'Scan_intent', 'Field_of_view', 'Largest_angular_scale', 'QA2_Status',\
         'Pub', 'Science_keyword', 'Scientific_category', 'ASA_PROJECT_CODE']



class queryCal:
    
    def __init__(self):
        pass
        # if readFile:
        #     self.listCal = self.readCal(fileCal, fluxrange)
        # else:
        #     self.listCal = []



    def readCal(self, ifile, fluxrange = [0.5, 1.0]):
        '''Read a list of calibrators in CSV format from the Source Catalogue web interface'''

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

                    if (flux >= fluxrange[0] and flux <= fluxrange[1]):
                        found = False
                        for nameYet in listcal:
                            if nameYet[0] == name:
                                found = True

                        if not found:
                            coord = coordinates.SkyCoord(alpha2000 + delta2000, unit=(u.hourangle, u.deg), equinox='J2000') # converted using astropy
                            listcal.append([name, coord.ra.value, coord.dec.value, flux])

        return(listcal)
    

    def queryObject(self, obj, theta = 0.005, public = True, savedb=False, dbname='calibrators.db'):
        '''query per object'''
        # obj is an array consist of at least [objname, ra, dec]
        # ra dec in deg
        
        region = coordinates.SkyCoord(obj[1], obj[2], unit='deg')
        try:
            alma_data = Alma.query_region_async(region, theta*u.deg, science = False, public =  public)

            # alma_data is in unicode
            xmldata = BeautifulSoup(alma_data.text, "lxml")

            if len(xmldata.contents) < 2:
                print "WARNING: ", obj[0], " query FAIL! [Empty response]"
                res = 0
            else:
                with open('xmldata.xml', 'w') as ofile: # write in a file, 
                    ofile.write(str(xmldata.contents[1]))

                tree = ET.parse("xmldata.xml") # read the file again. TODO: Find a better way!
                root = tree.getroot()

                projects = []
                for i, proj in enumerate(root[0][0][0][2][36][0]): # many child
                    projects.append([])
                    for data in proj:
                        if data.text is None:
                            projects[i].append('None')
                        else:
                            projects[i].append(data.text)

                #convert python list to Pandas DataFrame
                df = pd.DataFrame(projects, columns=columns)

                if len(df) == 0: # empty data
                    print "WARNING: ", obj[0], " query success, but no data!"
                    res = 0
                else:
                    # change the type of some columns
                    df['Spatial_resolution']    = df['Spatial_resolution'].astype(float)
                    df['Frequency_resolution']  = df['Frequency_resolution'].astype(float)
                    df['Integration']           = df['Integration'].astype(float)
                    df['Velocity_resolution']   = df['Velocity_resolution'].astype(float)
                    df['Field_of_view']         = df['Field_of_view'].astype(float)
                    df['Largest_angular_scale'] = df['Largest_angular_scale'].astype(float)

                    # save the query result in sql database with Table's name = calibrator's name
                    if savedb:
                        conn = sql.connect(dbname)
                        df.to_sql(obj[0], conn, if_exists='replace') # if table exist, just replace
                    
                    res = df

            # remove the temporary xml file is exist
            if os.path.exists('xmldata.xml'): 
                os.remove('xmldata.xml')

            return res

        except:
            print "WARNING: ", obj[0], " query FAIL! [No response from server]"
            return 0



    def queryAlma(self, listname, public = True, savedb=False, dbname='calibrators.db'):
        '''Query the data for the list of ALMA name (bulk)'''
        # name is an array consist of at least [objname, ra, dec]
        
        result = []
        for name in listname:
            df = self.queryObject(name, public = public, savedb=savedb, dbname=dbname)
            if isinstance(df, pd.DataFrame): # if not 0
                result.append([name, df])
            
        return(result)
            
              
    def selectDeepField(self, data, minTimeBand ={3:1e9,6:1e9,7:1e9}, maxFreqRes = 1000.0, selectPol=False, verbose = True):
        """
        From the queryAlma result we filter to select the deep field
        minTimeBand  : dictionary to select the minimum time  per band 
        maxFreqRes: maximum of the frequency resolution (in kHz)
        selectPol: is set True, only select uid with XY pol.
        
        Return the names and projects in a report with the date !!
        """
        
        finalReport = []
        # for every object/source
        for idx, item in enumerate(data):
            name      = item[0][0]
            alpha2000 = item[0][1]
            delta2000 = item[0][2]
            flux      = item[0][3]
            
            projects = []
            
            totalTime = {3:0., 4:0., 5:0., 6:0., 7:0., 8:0., 9:0., 10:0.}
            
            reportSource  = "\n\n#### Source : %s #####\n"%(name)
            reportSource += "#### No : %d \n"%(idx)
            reportSource += "#### Coord. 2000: %f  %f \n"%(alpha2000, delta2000)
            reportSource += "#### Flux : %f\n"%(flux)
            reportSource += "\n"
            
            uids = item[1]
            nuids = len(item[1])

            code        = uids['Project_code']
            source      = uids['Source_name']
            band        = uids['Band']
            integration = uids['Integration']
            frequency   = uids['Frequency_support']
            obsdate     = uids['Observation_date']
            res         = uids['Spatial_resolution']
            asdm        = uids['Asdm_uid']
            freqRes     = uids['Frequency_resolution']
            pol         = uids['Pol_products']
            

            # selection process
            # if selectPol == True: check the Polarization AND freqRes else: only check freqRes
            # if accepted: sum up the integration time
            for i in range(nuids):
                selectSG = False

                try: # for an error in ALMA data, Band ("6 7").. only found 1 data so far
                    bandi = int(band[i])
                except:
                    print "Ew, double band in ", code[i]
                    bandi = [int(x) for x in band[i].split(" ")]
                    

                if float(freqRes[i]) < maxFreqRes:
                        if selectPol:
                            if pol[i] =='XX XY YX YY':
                                selectSG = True
                        else:
                            selectSG = True

                if selectSG:
                    if isinstance(bandi, int):
                        totalTime[bandi] += integration[i]
                    else:
                        for band_i in bandi:
                            totalTime[band_i] += integration[i]


                if verbose and selectSG:
                    reportSource += "## %s  %20s   Band:%s   Obsdate:%s   FreqRes:%f   SpatRes:%f   Pol:%s    Asdm:%s  Integration:%f\n"%(code[i], source[i], band[i], obsdate[i], freqRes[i], res[i], pol[i], asdm[i], integration[i])
                
                foundProject = False
                for p in projects:
                    if p == code[i]:
                        foundProject = True
                    
                if not foundProject and selectSG:
                    projects.append(code[i])

            
            acceptedProject = len(projects) # number of accepted projects

            reportSource += "\n"
            for key in totalTime:
                if totalTime[key] > 0.:
                    reportSource += "Time Band %d : %6.0fs (%3.1fh) \n"%(key, totalTime[key], totalTime[key] / 3600.)
            
            # select source if number of uid > 0 AND totalTime > minTimeBand
            # add report to the finalReport
            selectSource = False
            if acceptedProject > 0:
                selectSource = True
                reportSource += "\n Project codes: \n"
                for p in projects:
                    reportSource += " %s "%(p)
                
                reportSource += "\n\n"
                reportSource += "Total uids: %d \n"%(nuids)
                reportSource += "Total accepted projects: %d"%(acceptedProject)
                reportSource += "\n"

                for k in minTimeBand:
                    if totalTime[k] < minTimeBand[k]:
                        selectSource = False
            
                    
            if selectSource:
                finalReport.append([nuids, item[0], reportSource])
                
        # sorting according to the number of uids
        finalReportSorted = sorted(finalReport, key=lambda data: data[0])
        
        return(finalReportSorted)



    def selectDeepfield_fromsql(self, dbname, maxFreqRes = 1000.0, array='12m', \
        excludeCycle0=True, selectPol=False, minTimeBand = {3:1e9,6:1e9,7:1e9}, \
        verbose = True, silent=True):

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
            subquery = "SELECT * FROM [{0}] WHERE (Array LIKE '%{1}%') AND (Frequency_resolution < {2})".format(tab, array, maxFreqRes)

            if selectPol:
                subquery += " AND (Pol_products LIKE '%XY%')"

            if excludeCycle0:
                subquery += " AND (Project_code NOT LIKE '2011.%')"
            
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


            if selectSource:
                if not silent:
                    print 'Accepted', tab, totalTime
                
                nsource += 1

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

                        sqlcmd  = "SELECT Project_code, Source_name, RA, Dec, "
                        sqlcmd += "Band, Integration, Observation_date, Spatial_resolution, "
                        sqlcmd += "Asdm_uid, Frequency_resolution, Array, Pol_products " + subquery[8:]
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
                                reportSource += "\n{0} {1} ra:{2} dec:{3} Band:{4} Int:{5} Obsdate:{6} SpatRes:{7} ASDM:{8} FreqRes:{9} Array:{10} Pol:{11}".format(uid[0], uid[1], uid[2][:9], uid[3][:9], \
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



    def writeReport(self, report, file = "deepfieldRG.txt", silent=True):
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




    

##############Main program########################
#### example ....#######################
if __name__=="__main__":
    " main program"   
    
    fileCal = "callist_20170329.list"
    q       = queryCal()
    listCal = q.readCal(fileCal, fluxrange = [0.5, 10000000])
    data    = q.queryAlma(listCal, public = True, savedb=True, dbname='calibrators_gt_0.5Jy.db')
    report = q.selectDeepfield_fromsql("calibrators_gt_0.5Jy.db", maxFreqRes=10000000, array='12m', \
        excludeCycle0=True, selectPol=False, minTimeBand={3:3600., 6:3600., 7:3600.}, verbose=True, silent=True)
    q.writeReport(report, "report.txt", silent=True)