#!/home/ridlo/anaconda2/bin/python

# python version of http://ned.ipac.caltech.edu/forms/denv.html

"""
HISTORY:

2017.03.10:
    - first stab

2017.03.14:
    - add function to query from position
    - add function to plot all data even without redshift information (cone only)
    
2017.03.15:
    - put a circle in plot_cone (SL)

BUG:
    + ERROR if 'theta' (opening angle of cone) or the number of data is too large
        error message: Exception: Query failed: HTTPConnectionPool(host='ned.ipac.caltech.edu', port=80): Read timed out.

    + Web service is so much faster, can't do anything

RUN: 

"""

# copyleft
__author__="Ridlo W. Wibowo @ ITB"
__version__="0.1.0@2017.03.10"


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck15 # FlatLambdaCDM with parameters from Planck Collab 2015, Paper XIII
from astropy import coordinates
from astropy.constants import c

from astroquery.ned import Ned


class Galenv:
    def __init__(self):
        self.CosmoModel = Planck15
        self.C = c.value/1000.0 # c in km/s
        self.cone = 0.0
        self.conedv = 0.0


    def queryobject_byname(self, objname):
        obj_table = Ned.query_object(objname)
        
        z = obj_table[0]['Redshift']
        v0 = obj_table[0]['Velocity']
        ra = obj_table[0]['RA(deg)']
        dec = obj_table[0]['DEC(deg)']

        return z, v0, ra, dec


    def theta_from_tangen_dist_and_dist(self, tangen_dist, dist):
        theta = np.degrees(float(tangen_dist)/dist) # Flat Universe
        return theta


    def calc_dA_theta(self, z, tangen_dist):
        angular_diameter_dist = self.CosmoModel.angular_diameter_distance(z)
        theta = self.theta_from_tangen_dist_and_dist(tangen_dist, angular_diameter_dist.value)
        return angular_diameter_dist.value, theta

    
    def z_from_vel(self, v):
        return np.sqrt((1.0 + v/self.C)/(1.0 - v/self.C)) - 1.0


    def searchobject_in_cone(self, coord, theta, withz=True):
        # coord in astropy.coordinates
        result = Ned.query_region(coord, radius=theta*u.deg, equinox=coord.equinox) 
        if withz: # select objects which has redshift information
            res = result[~result['Redshift'].mask] # return False if data is masked (no redshift)
        else:
            res = result # all objects

        self.cone = res
        print "Number of objects in cone (only if z exist = "+str(withz)+"): ", len(res)
        return res


    def searchobject_in_cone_dv(self, coord, theta, v0, dv, withz=True):
        res = self.searchobject_in_cone(coord, theta, withz=withz)
        vmax, vmin = v0+dv, v0-dv
        cut1 = res[res['Velocity'] > vmin]
        cut2 = cut1[cut1['Velocity'] < vmax]

        self.conedv = cut2
        print "Number of objects in cone that has data of velocity v +- dv: ", len(cut2)
        return cut2


    def plot_env(self, coord, theta, res, v0, dv=5000, xSize=7.5, ySize=6, title='', show=True, savefig=False, imgname="plot.png"):
        ra = coord.ra.value
        dec = coord.dec.value

        fig = plt.figure(figsize=(xSize, ySize))        
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
        
        ax = plt.subplot(gs[0])
        ax.axis('equal')
        limangle = 1.15*theta
        ax.set_xlim((ra-limangle, ra+limangle))
        ax.set_ylim((dec-limangle, dec+limangle))
        
        # it is wrong if I draw a circle around (ra, dec) with radius theta
        # due to small circle in celestial sphere for DEC

        #circle = plt.Circle((ra, dec), theta, fc='none', ec='black')
        #ax.add_artist(circle)
        ax.plot(res['RA(deg)'], res['DEC(deg)'], 'b.')
        plt.gca().invert_xaxis() # RA from E to W
        ax.set_xlabel('RA(deg)')
        ax.set_ylabel('DEC(deg)')
        plt.title(title)
        
        ###
        
        ax2 = plt.subplot(gs[1])
        ax2.plot(res['RA(deg)'], res['Velocity'] - v0, 'b.')
        ax2.set_ylim((-dv, dv))
        ax2.set_xlim((ra-limangle, ra+limangle))
        ax2.yaxis.tick_right()
        plt.gca().invert_xaxis() 
        ax2.set_xlabel('RA(deg)')
        ax2.set_ylabel('Relative velocity')

        fig.tight_layout()

        if savefig:
            plt.savefig(imgname)

        if show:
            plt.show()

        plt.close()


    def plot_cone(self, coord, theta, res, xSize=7.5, ySize=6, title='', show=True, savefig=False, imgname="plot.png"):
        ra = coord.ra.value
        dec = coord.dec.value

        fig = plt.figure(figsize=(xSize, ySize))        
        gs = gridspec.GridSpec(1, 1)
        
        ax = plt.subplot(gs[0])
        ax.axis('equal')
        limangle = 1.15*theta
        ax.set_xlim((ra-limangle, ra+limangle))
        ax.set_ylim((dec-limangle, dec+limangle))
        
        # it is wrong if I draw a circle around (ra, dec) with radius theta
        # due to small circle in celestial sphere for DEC

        #circle = plt.Circle((ra, dec), theta, fc='none', ec='black')
        #ax.add_artist(circle)
        ax.plot(res['RA(deg)'], res['DEC(deg)'], 'b.')
        plt.gca().invert_xaxis() # RA from E to W
        ax.set_xlabel('RA(deg)')
        ax.set_ylabel('DEC(deg)')
        plt.title(title)
        circle1= plt.Circle((ra, dec), theta, color='blue', fill = False)
        
        ax.add_artist(circle1)
        
        fig.tight_layout()

        if savefig:
            plt.savefig(imgname)

        if show:
            plt.show()

        plt.close()


    def conedv_byname(self, objname, tangen_dist, dv, show=False, savefig=False, imgname="plot.png"):
        print "Query from name: ", objname
        z, v0, ra, dec = self.queryobject_byname(objname)
        coord = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
        print "RA, Dec (deg)  : ", ra, dec
        print "Redshift : ", z
        print "Velocity : ", v0

        
        if isinstance(z, np.float64) or isinstance(v0, np.float64): # if there is v or z in NED
            dA, theta = self.calc_dA_theta(z, tangen_dist)
            print "Angular-diameter distance : ", dA, " Mpc"
            print "Angular radius for "+str(tangen_dist)+" Mpc seen from distance "+str(dA)+" Mpc : ", theta, " deg"
        else:
            print objname + " does not currently have any redshift measurements in NED."
            print "Please provide your desired radial velocity or define a Metric Distance for the search center."
            
            v0 = input("Radial velocity? ")
            
            if v0 < 3000.0:
                print "For nearby objects (below 3000 km/s), you must supply a Metric Distance to the center of search."
                dist = input("Metric distance? ")
                print "User Defined Metric Distance: ", dist, " Mpc"

            else:
                z = self.z_from_vel(v0)
                print "Calculated z = ", z
                dist, theta = self.calc_dA_theta(z, tangen_dist)
                print "Angular-diameter distance : ", dist, " Mpc"

            
            theta = self.theta_from_tangen_dist_and_dist(tangen_dist, dist)
            print "Angular radius for "+str(tangen_dist)+" Mpc seen from distance "+str(dist)+" Mpc : ", theta, " deg"


        print "Searching objects in cone section... \ntheta = "+str(theta)+" deg; v = ["+str(v0-dv)+", "+str(v0+dv)+"]" 
        res = self.searchobject_in_cone_dv(coord, theta, v0, dv)
        #print "Number of objects : ", len(res)

        self.plot_env(coord, theta, res, v0, dv, title=str(tangen_dist)+" Mpc, "+str(dv)+" km/s, " 
                + "\n$H_0$ = " + str(self.CosmoModel.H0)
                + ", $\Omega_{m}$ = " + str(self.CosmoModel.Om0)[:5]
                + ", $\Omega_{\Lambda} = $"+ str(self.CosmoModel.Ode0)[:5], show=show, savefig=savefig, imgname=imgname)


    def conedv_byposvel(self, coord, v0, tangen_dist, dv, show=False, savefig=False, imgname="plot.png"):
        # position in in astropy.coordinates
        print "Query from position: ", coord
        print "Radial velocity : ", v0

        if v0 < 3000.0 :
            print "For nearby objects (below 3000 km/s), you must supply a Metric Distance to the center of search."
            dist = input("Metric distance? ")
            print "User Defined Metric Distance: ", dist, " Mpc"
            theta = self.theta_from_tangen_dist_and_dist(tangen_dist, dist)
            print "Angular radius for "+str(tangen_dist)+" Mpc seen from distance "+str(dist)+" Mpc : ", theta, " deg"
        else:
            z = self.z_from_vel(v0)
            print "Calculated z = ", z
            dA, theta = self.calc_dA_theta(z, tangen_dist)
            print "Angular-diameter distance : ", dA, " Mpc"
            print "Angular radius for "+str(tangen_dist)+" Mpc seen from distance "+str(dA)+" Mpc : ", theta, " deg"



        print "Searching objects in cone section...\ntheta = "+str(theta)+" deg; v = ["+str(v0-dv)+", "+str(v0+dv)+"]" 
        res = self.searchobject_in_cone_dv(coord, theta, v0, dv)
        #print "Number of objects : ", len(res)

        self.plot_env(coord, theta, res, v0, dv, title=str(tangen_dist)+" Mpc, "+str(dv)+" km/s, " 
            + "\n$H_0$ = " + str(self.CosmoModel.H0)
            + ", $\Omega_{m}$ = " + str(self.CosmoModel.Om0)[:5]
            + ", $\Omega_{\Lambda} = $"+ str(self.CosmoModel.Ode0)[:5], show=show, savefig=savefig, imgname=imgname)


    def cone_byname(self, objname, theta, show=False, savefig=False, imgname="plot.png"):
        print "Query from name: ", objname
        z, v0, ra, dec = self.queryobject_byname(objname)
        coord = coordinates.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
        print "RA, Dec (deg)  : ", ra, dec

        print "Searching objects in cone (only).. \ntheta = "+str(theta)+" deg"
        res = self.searchobject_in_cone(coord, theta, withz=False)

        self.plot_cone(coord, theta, res)


    def cone_bypos(self, coord, theta, show=False, savefig=False, imgname="plot.png"):
        print "Query from position ", coord

        print "Searching objects in cone (only).. \ntheta = "+str(theta)+" deg"
        res = self.searchobject_in_cone(coord, theta, withz=False)

        self.plot_cone(coord, theta, res)




if __name__ == '__main__':

    # objname = "ACO 13"
    # tangential_dist = 10.0 # Mpc
    # dv = 5000.0 # km/s

    # we can also use command line arguments
    # Example:
    # ./galenv.py -i ACO\ 13 -td 10 -dv 5000 -show -savefig Abell.png

    import sys
    ge = Galenv()

    show = False
    savefig = False
    imgname = "plot.png"

    while len(sys.argv) > 1:
        option = sys.argv[1]; del sys.argv[1]
        if option == '-i':
            objname = sys.argv[1]; del sys.argv[1]
        elif option == '-td':
            tangential_dist = float(sys.argv[1]); del sys.argv[1]
        elif option == '-dv':
            dv = float(sys.argv[1]); del sys.argv[1]
        elif option == '-show':
            show=True
        elif option == '-savefig':
            savefig=True; imgname = sys.argv[1]; del sys.argv[1]
        else:
            print sys.argv[0], ': invalid option', option
            sys.exit(1)

    ge.conedv_byname(objname, tangential_dist, dv, show=show, savefig=savefig, imgname=imgname)



