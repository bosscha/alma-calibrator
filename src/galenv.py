#!/home/ridlo/anaconda2/bin/python

# python version of http://ned.ipac.caltech.edu/forms/denv.html

"""
HISTORY:

2017.03.10:
    - first stab


RUN:

"""

# copyleft
__author__="Ridlo W. Wibowo @ ITB"
__version__="0.1.0@2017.03.10"


import numpy as np
import matplotlib.pyplot as plt

#from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import Planck15 # FlatLambdaCDM with parameters from Planck Collab 2015, Paper XIII
import astropy.units as u
from astroquery.ned import Ned
import matplotlib.gridspec as gridspec



def queryobject_by_name(objname):
    obj_table = Ned.query_object(objname)
    
    z = obj_table[0]['Redshift']
    v0 = obj_table[0]['Velocity']
    ra = obj_table[0]['RA(deg)']
    dec = obj_table[0]['DEC(deg)']

    return z, v0, ra, dec


def calculate_dA_theta(z, tangential_dist):
    angular_diameter_distance = Planck15.angular_diameter_distance(z) # in Mpc
    # Flat universe
    theta = np.degrees(tangential_dist/angular_diameter_distance * u.rad) # in degree
    
    return angular_diameter_distance.value, theta.value



def searchobject_in_cone_dv(objname, theta, v0, dv):
    # cone search, with theta
    result = Ned.query_region(objname, radius=theta*u.deg) 

    # select objects which has redshift information
    # return False if data is masked (no redshift) --> mainly from GALEXASC
    only_with_z = result[~result['Redshift'].mask] 

    # select objects which has velocity in range of +- dv
    vmax, vmin = v0+dv, v0-dv
    cut1 = only_with_z[only_with_z['Velocity'] > vmin]
    res = cut1[cut1['Velocity'] < vmax]
    
    return res



def plot_env(ra, dec, theta, res, v0, dv=5000, xSize=7.5, ySize=6, title='', show=True, savefig=False, imgname="plot.png"):
    fig = plt.figure(figsize=(xSize, ySize))
    
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
    
    ax = plt.subplot(gs[0])
    ax.axis('equal')
    limangle = 1.15*theta
    ax.set_xlim((ra-limangle, ra+limangle))
    ax.set_ylim((dec-limangle, dec+limangle))
    
    circle = plt.Circle((ra, dec), theta, fc='none', ec='black')
    ax.add_artist(circle)
    ax.plot(res['RA(deg)'], res['DEC(deg)'], 'b.')
    plt.gca().invert_xaxis() # RA from E to W
    plt.title(title)
    
    ###
    
    ax2 = plt.subplot(gs[1])
    ax2.plot(res['RA(deg)'], res['Velocity'] - v0, 'b.')
    ax2.set_ylim((-dv, dv))
    ax2.set_xlim((ra-limangle, ra+limangle))
    ax2.yaxis.tick_right()
    plt.gca().invert_xaxis() 
    
    fig.tight_layout()

    if savefig:
        plt.savefig(imgname)

    if show:
        plt.show()

    plt.close()



def run(objname, tangential_dist, dv, show=False, savefig=False, imgname="plot.png"):
    print "Query for object : ", objname
    z, v0, ra, dec = queryobject_by_name(objname)
    print "Redshift : ", z
    print "Velocity : ", v0
    print "RA, Dec (deg)  : ", ra, dec

    dA, theta = calculate_dA_theta(z, tangential_dist*u.Mpc)
    print "Angular diameter distance : ", dA
    print "Angular radius for "+str(tangential_dist)+" Mpc : ", theta, " deg"

    print "Searching objects in cone section..."
    res = searchobject_in_cone_dv(objname, theta, v0, dv)
    print "Number of objects : ", len(res)

    plot_env(ra, dec, theta, res, v0, dv, title=str(tangential_dist)+" Mpc, "+str(dv)+" km/s, " 
            + "\n$H_0$ = " + str(Planck15.H0)
            + ", $\Omega_{m}$ = " + str(Planck15.Om0)[:5]
            + ", $\Omega_{\Lambda} = $"+ str(Planck15.Ode0)[:5], show=show, savefig=savefig, imgname=imgname)



if __name__ == '__main__':

    # objname = "ACO 13"
    # tangential_dist = 10.0 # Mpc
    # dv = 5000.0 # km/s

    # run(objname, tangential_dist, dv, show=True, savefig=True, imgname="Abell.png")

    # we can also use command line arguments
    # Example:
    # ./galenv.py -i ACO\ 13 -td 10 -dv 5000 -show -savefig "Abell.png"

    import sys

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

    run(objname, tangential_dist, dv, show=show, savefig=savefig, imgname=imgname)