### Simulations of point sources in a field with a give array configuration
###

import numpy.random as rd
import shutil, os

rd.seed()
  

simDir = "./"
  
def saveSourceList(sourcefile, al, de, flux):
    "Write the source list with flux"
    
    with open(sourcefile,"a") as f: 
        f.write("%f  %f  %f \n"%(al,de,flux*1e3))
                
    
def createSource(nsource, fluxmin, fluxmax):
    "Create a CL with nsource located randomly"
    
    clfile = simDir+"/source.cl"
    
    if os.path.exists(clfile):
        shutil.rmtree(clfile, ignore_errors=True)

    sourcefile = simDir + "/sources.txt"
    with open(sourcefile,"w") as f:
        f.write("## Alpha   Delta  Flux (mJy) \n")

    al  = rd.uniform(-40.,40., nsource)
    de =  rd.uniform(-20.,20.,nsource)
    
    fluxsource = rd.uniform(fluxmin, fluxmax,nsource)
    diskSize = "0.001arcsec"   
    
    for i in range(nsource):
        
        alpha =  10. + al[i]/3600.
        delta = -50. + de[i]/3600.
        decString="J2000 %6.4fdeg %6.4fdeg"%(alpha,delta)
        minmaj = diskSize

        cl.addcomponent(dir = decString, flux = fluxsource[i], fluxunit = 'Jy', freq = "100GHz" ,shape = "Gaussian", majoraxis= minmaj, minoraxis = minmaj, positionangle = "0deg")
        saveSourceList(sourcefile, alpha, delta, fluxsource[i])
        
    ## ccentral source
    minmaj = "0.001arcsec"
    alpha =  10.
    delta = -50.
    decString="J2000 %6.4fdeg %6.4fdeg"%(alpha,delta)
    cl.addcomponent(dir = decString, flux = 1.0, fluxunit = 'Jy', freq = "100GHz" ,shape = "Gaussian", majoraxis= minmaj, minoraxis = minmaj,positionangle = "0deg")
    saveSourceList(sourcefile, alpha, delta, 1.0)   
        
        
    cl.rename(clfile)
    cl.done()

    
    
def simulation(antcfg, projectname, totime ,ovrwrt= True):
    
    simalma(
        overwrite      = ovrwrt, 
        dryrun         = False,
        project        = projectname,
        compwidth      = '800MHz' ,
        antennalist    =  antcfg ,
        totaltime      = totime,
        pwv            = 10.0,
        mapsize        =  "20arcsec",
        niter          =   2000,
        threshold      = '0.01mJy', 
        cell           = "0.1arcsec",
        imsize         = 500)


################################Main######################################################3
project = "psfield"
antcfg = "alma.cycle2.6.cfg"
inttime = "1200s"

fluxmin = 1e-4
fluxmax = 1e-3

createSource(50, fluxmin, fluxmax)    
simulation(antcfg, project, inttime, True)
