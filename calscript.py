# script for calibrator continuum imaging

import os
import sys
import glob
import numpy as np

sys.path.append('$DATADIR/analysis_scripts') #my directory of analysis_script # analysisUtils
import analysisUtils as aU
es = aU.stuffForScienceDataReduction()

# median freq of ALMA Band
ALMA_BAND = {1: 38.15e9, 
             2: 78.50e9, 
             3: 100.0e9,
             4: 144.0e9,
             5: 187.0e9,
             6: 243.0e9,
             7: 324.0e9,
             8: 442.5e9,
             9: 661.0e9,
             10: 868.5e9}

listband = { 1: "ALMA_RB_01",
             2: "ALMA_RB_02",
             3: "ALMA_RB_03",
             4: "ALMA_RB_04",
             5: "ALMA_RB_05",
             6: "ALMA_RB_06",
             7: "ALMA_RB_07",
             8: "ALMA_RB_08",
             9: "ALMA_RB_09",
             10: "ALMA_RB_10" }

tsys_nominal = {3: 75.0, 
                4: 86.0, 
                5: 120.0, 
                6: 90.0, 
                7: 150.0, 
                8: 387.0, 
                9: 1200.0, 
                10:1515.0}
                
sensitivities = {3: 0.20, 
                 4: 0.24, 
                 5: 0.37, 
                 6: 0.27, 
                 7: 0.50, 
                 8: 1.29, 
                 9: 5.32, 
                 10: 8.85}

# in SI
# k = 1.380658e-23

def is_power2(num): # num must be integer
    return (num & (num-1)) == 0


def estimate_freq(msName):
    """ only using Band, we can use central freq from spectral window information"""
    tb.open(msName + "/ASDM_RECEIVER")
    bnd = tb.getcol("frequencyBand")
    bandused = 0
    for b in listband:
        if np.any(bnd == listband[b]):
            bandused = b
            print "Band : ", bandused
    tb.close() 
    freq = ALMA_BAND[bandused]
    return freq # in Hz


def estimateSynthesizedBeam_perfield(msName, freq, field):
    """useL80method, modified from analysisUtils"""
    result = aU.getBaselineStats(msName, field=field, percentile=80) # field can be more than one... -_- choose field field with larger baseline?
    L80 = result[0] # meters
    arcsec = 3600*np.degrees(1)*0.574*aU.c_mks/(estimate_freq(msName)*L80)

    return arcsec


def estimateSynthesizedBeam(msName, freq, field=''):
    if len(field) > 1: # if there is more than one field (combine)
        fields = field.split(",")
        synthesizedbeam = []
        for fi in fields:
            synthesizedbeam.append(estimateSynthesizedBeam_perfield(msName, freq, fi))
        arcsec = min(synthesizedbeam) # choose the minimum one
    else:
        arcsec = estimateSynthesizedBeam_perfield(msName, freq, field)

    return arcsec


def estimate_cell_and_imsize(msName, field):
    print("Estimate frequency from Band ...")
    freq = estimate_freq(msName)
    print("Estimated Freq: %s GHz" % str(freq/1e9))

    print("Calculate spatial resolution ...")
    spatial_resolution = estimateSynthesizedBeam(msName, freq, field=field)
    print("Estimated spatial resolution: %f arcsec" % (spatial_resolution))
    
    cell = 0.2*spatial_resolution
    cell = '%.2f' % cell + 'arcsec' # set cell to spatRes/5
    print("Select cell size: {0}".format(cell))

    print("Calculate primary beam ...")
    primary_beam = aU.primaryBeamArcsec(frequency=freq, diameter=12.0, \
        taper=10.0, obscuration=0.75, verbose=False, showEquation=True, \
        use2007formula=True, fwhmfactor=None) # return primary beam in arcsec  
    print("Primary beam: %f" % (primary_beam) )

    imsize_onehalf = int(1.5*primary_beam/(0.2*spatial_resolution)) # set imsize to 2.0 primary beam
    imsize_two = int(2.0*primary_beam/(0.2*spatial_resolution))

    imsize = imsize_onehalf
    unfactorize = True
    while unfactorize:
        if is_power2(imsize) or (imsize % 10 == 0 and is_power2(imsize/10)): # casa want even number and factorizeable of 2,3,5,7 (?) -> 2^n x 10
             unfactorize = False
        else: imsize += 1

    
    imz = raw_input("The primary beam multiple -> 1.5x = {0}, 2x = {1}, suggested imsize (default) = [{2}]. Enter your imsize: ".format(imsize_onehalf, imsize_two, imsize))
    if imz == "":
        imz = imsize
    else:
        imz = int(imz)

    print("Select imsize: {0}".format(imz))

    return cell, imz


def shortresume_ms(msName):
    intentSources = es.getIntentsAndSourceNames(msName)
    for key in intentSources:
        intent = key
        fieldid = intentSources[key]['id']
        names = intentSources[key]['name']
        if len(fieldid) > 0:
            print "Intent: ", intent, "               Field id: ", fieldid, "                  Name: ", names


def find_number_and_center_of_antenna(msName):
    tb.open(msName + "/ANTENNA")
        
    dish_diameter = tb.getcol("DISH_DIAMETER") # check, all is 12 meter
    seven = dish_diameter[dish_diameter < 10.]
    if len(seven) > 0:
        print("There is 7m array!")
        return

    antenna = tb.getcol("NAME")
    nant = len(np.unique(antenna))

    pos = tb.getcol("POSITION")
    center = [pos[0].mean(), pos[1].mean(), pos[2].mean()]
    dist = np.sqrt((pos[0]-center[0])**2 + (pos[1]-center[1])**2 + (pos[2]-center[2])**2)
    centerantenna = antenna[dist.argmin()]
    tb.close()
    return nant, centerantenna


def get_spwedge(msName):
    tb.open(msName + "/SPECTRAL_WINDOW")
    nspw = tb.nrows()
    spwedge = []
    for i in range(nspw):
        chanfreq = tb.getcell("CHAN_FREQ", i)
        chanwidth = np.abs(tb.getcell("CHAN_WIDTH", i)[0]) # assume all chanwidth in one spw is same
        ## we found some weird spw: freq is increasing and chanwidth is positive
        low, high = 0, 0
        if chanfreq[-1] < chanfreq[0]:
            low = chanfreq[-1] - chanwidth/2.0
            high = chanfreq[0] + chanwidth/2.0 
        else:
            low = chanfreq[0] - chanwidth/2.0
            high = chanfreq[-1] + chanwidth/2.0 

        spwedge.append([low, high])

    tb.close()
    print "Spw edge:", spwedge

    return spwedge


def get_npolarizations(msName):
    tb.open(msName)
    nPolarizations = np.shape(tb.getcell("DATA",0))[0]
    tb.close()
    return nPolarizations


def estimate_sensitivity(msName):
    # find integration time
    time_os = au.timeOnSource(msName)
    print("-------------------------------------------------")
    int_time = time_os[0]['minutes_on_source'] # always in field 0 (from split)
    print "Integration time = ", int_time, " min"

    nant, centerantenna = find_number_and_center_of_antenna(msName) # find number of antenna
    npol = get_npolarizations(msName) # number of polarizations
    spwedge = get_spwedge(msName)

    bandunion = [] # union from that some intervals
    for begin,end in sorted(spwedge): # sort spwedge base on lower limit of spw
        if bandunion and bandunion[-1][1] >= begin: # check the intersect
            bandunion[-1][1] = max(bandunion[-1][1], end) # combine
        else:
            bandunion.append([begin, end]) #seperate interval

    print "Merged spw: ", bandunion
    print bandunion[0][0], bandunion[-1][1]
    print "Central freq: ", (bandunion[0][0] + bandunion[-1][1])/2.0
    totalbandwidth = 0.0
    for interval in bandunion:
        totalbandwidth += (interval[1] - interval[0])

    totalbandwidth = totalbandwidth*1.0e-9 # in GHz
    print "Bandwidth = ", totalbandwidth, " GHz"
    print "Number of antenna = ", nant
    print "Number of polarization = ", npol
    
    # find band
    tb.open(msName + "/ASDM_RECEIVER")
    bnd = tb.getcol("frequencyBand")
    bandused = 0
    for b in listband:
        if np.any(bnd == listband[b]):
            bandused = b
            print "Band : ", bandused
    tb.close() 


    try:
        tb.open(msName + "/ASDM_CALATMOSPHERE")
        tSys = tb.getcol("tSys")
        tb.close()
        # estimate Tsys
        tsys = tSys.flatten() # 1D array
        tsys0 = tsys[tsys > 0.0] # remove 0's part
        meantsys0 = np.mean(tsys0)
        std5 = 5.0*np.std(tsys0)
        tsys5 = tsys0[abs(tsys0 - meantsys0) < std5] # remove outliers, value > 5 sigma
        meanTsys = np.mean(tsys5)
        print "mean Tsys = ", meanTsys
    except:
        print "Using standard Tsys for this band"
        meanTsys = tsys_nominal[bandused]



    #  Assume apropri sensitivity for band 3,4,6,7,8,9,10
    #  Sensitivity units are in mJy and are normalized to:
    #    16 12-m antennas,
    #    8 GHz bandwidth
    #    Dual freq bandwidth
    #    for tsys_nominal given above
    #    Integration time of one minute
    apriori_sensitivity = sensitivities[bandused]

    # sigma = 2*k*Tsys/(aperture_eff * np.sqrt(nant*(nant-1) * npol * bandwidth * int_time))
    # scaling
    scalingTsys = meanTsys/tsys_nominal[bandused]
    scalingInttime_bandwidth_npol = 1.0/np.sqrt(int_time * totalbandwidth/8.0 * npol/2.0)
    scalingNant = np.sqrt(240.0/(nant*(nant-1))) # 240 = 16*(16-1)
    scalingfactor = scalingTsys * scalingInttime_bandwidth_npol * scalingNant

    sensitivity = apriori_sensitivity * scalingfactor

    print "Estimated sensitivity: ", sensitivity, " mJy"
    

    return sensitivity # in mJy


def print_last15lines_fromlastlog():
    # find more robust method! (e.g. using awk, grep or regex in python)
    lsf = sorted(os.listdir("./"))
    casalogg = []
    for l in lsf:
        if 'casa-' in l: casalogg.append(l)
    
    # open the last log
    with open(casalogg[-1], "r") as logg:
        lines = logg.readlines()
    
    print("Last 15 lines from CASA Log:")
    for line in lines[-15:]:
        print line[32:],   


def findrms_fromresidual(imgname):
    residualimg = imgname + ".residual"
    xstat = imstat(residualimg)
    print "RMS from residual-image " + imgname + " (using imstat): ", xstat['rms'][0]
    return xstat['rms'][0]


def object_split():
    print("\n......................................................................")
    print("Split the field from .clb ")
    print("List all *.ms.cal.clb files")
    listMS = []
    for ifile in glob.glob("*.ms.split.cal.clb"):
        print(ifile)
        listMS.append(ifile)

    for msName in listMS:
        print "Ready to split ", msName
        shortresume_ms(msName)
        
        ans = raw_input("Is there any interesting calibrator that you want to split? [y/n] ")
        if ans.lower()[0] == 'y':
            field = raw_input("Enter field id: ")
            objname = raw_input("Object name: ")
            outputvis = msName + '.field_'+ field + "." + objname +".ms"
            width = raw_input("Enter the width [default=32]: ")
            if width == "":
               width = 32
            else:
               width = int(width)
            print("Begin split...")
            split(vis=msName,outputvis=outputvis,keepmms=True,field=field,spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="data",keepflags=True,width=width,timebin="0s",combine="")
            print("End split. Result: {0}".format(outputvis))

            askagain = 1
            while askagain:
                ans = raw_input("Any other field in this MS? [y/n] ")
                if ans.lower()[0] == 'y':
                    field = raw_input("Enter field id: ")
                    objname = raw_input("Object name: ")
                    outputvis = msName + '.field_'+ field + "." + objname +".ms"
                    width = raw_input("Enter the width [default=32]: ")
                    if width == "":
                       width = 32
                    else:
                       width = int(width)
                    print("Begin split...")
                    split(vis=msName,outputvis=outputvis,keepmms=True,field=field,spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="data",keepflags=True,width=width,timebin="0s",combine="")
                    print("End split. Result: {0}".format(outputvis))
                else:
                    askagain = 0
        else:
            print("Next MS (if any)")
            continue


def continuum(msName, field='', imagename='', phasecenter='', niter=1000, threshold='0.0mJy', interactive=True):
    """imaging the continuum"""
    
    if imagename == '':
        imagename = msName + '.cont' # default is the msName.cont

    cell, imsize = estimate_cell_and_imsize(msName, field)
    sensitivity = estimate_sensitivity(msName)
    if threshold == '0.0mJy' or threshold == '':
        threshold = raw_input("Enter your threshold [0.0mJy]: ")

    if phasecenter == '':
        phasecenter = raw_input("Enter your phasecenter if you want ['']: ")


    print("Begin cleaning task...")
    clean(vis=msName,imagename=imagename, outlierfile="",field=field,spw="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",\
        intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,\
        painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",\
        niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",ftmachine="mosaic",\
        mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,interactive=interactive,\
        mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cell,\
        phasecenter=phasecenter,restfreq="",stokes="I",weighting="briggs",robust=0.5,uvtaper=False,outertaper=[''],\
        innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",\
        npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,\
        allowchunk=False)
    print("End cleaning task.")

    # print 15 last line in clean log
    print_last15lines_fromlastlog()

    # print rms from residual image, imstat
    findrms_fromresidual(imagename)



def selfcal_cycle(msName, cycle, field, niter, threshold, cell, imsize, phasecenter, interactive, solint, solnorm, refant, calmode):
    #clean
    if cycle == 1:
        vis = msName
        field = field
    else:
        vis = msName+'.self'+str(cycle-1)
        field = '0'

    imgname = msName+".cont"+str(cycle)

    print("Cleaning...")
    clean(vis=vis,imagename=imgname,outlierfile="",field=field,spw="",selectdata=True,timerange="",uvrange="",antenna="",
        scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",
        rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",
        interpolation="linear",niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",
        ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,
        interactive=interactive,mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",
        imsize=imsize,cell=cell,phasecenter=phasecenter,restfreq="",stokes="I",
        weighting="briggs",robust=0.5,uvtaper=False,outertaper=[''],innertaper=['1.0'],
        modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",
        npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,
        flatnoise=True,allowchunk=False)

    # print 15 last line in clean log
    print_last15lines_fromlastlog()

    # print rms from residual image, imstat
    findrms_fromresidual(imgname)

    #gaincal
    print("gaincal...")
    gaincal(vis=vis,caltable=msName+".G"+str(cycle),field=field,spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",
        scan="",observation="",msselect="",solint=solint,combine="",preavg=-1.0,refant=refant,
        minblperant=4,minsnr=3.0,solnorm=solnorm,gaintype="G",smodel=[],calmode=calmode,append=False,splinetime=3600.0,
        npointaver=3,phasewrap=180.0,docallib=False,callib="",gaintable=[],gainfield=[],interp=[],spwmap=[],parang=False)


    #check gain phase
    #print("plotcal... check")
    #plotcal(caltable=msName+".G"+str(cycle),xaxis="time",yaxis="phase",poln="",field=field,antenna="",spw="",timerange="",
    #    subplot=331,overplot=False,clearpanel="Auto",iteration="antenna",plotrange=[],showflags=False,plotsymbol="o",
    #    plotcolor="blue",markersize=5.0,fontsize=10.0,showgui=False,figfile="plotcal_"+str(cycle)+".png")

    #applycal
    print("applycal...")
    applycal(vis=vis,field=field,spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",
        msselect="",docallib=False,callib="",gaintable=msName+".G"+str(cycle),gainfield=[],
        interp="nearest",spwmap=[],calwt=[True],parang=False,applymode="",flagbackup=True)

    #split
    print("split...")
    split(vis=vis,outputvis=msName+".self"+str(cycle),keepmms=True,field=field,spw="",scan="",antenna="",correlation="",
        timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="corrected",keepflags=True,width=1,
        timebin="0s",combine="")

    print("Finish Cycle {0}".format(cycle))


############################################################################

########## SPLIT FOR ALL CALIBRATOR
# in case someday we have to remove the science target data
# width = 1
def cal_split():
    #os.chdir(".") # enter the directory
    print("List all *.ms.split.cal files in this directory.")
    listms = []
    for ifile in glob.glob("*.ms.split.cal"):
        print(ifile)
        listms.append(ifile)

    print("\n......................................................................")
    for msName in listms:
        print("Check calibrator field using listobs.. " + msName + ", please check the field in Log or look in short-resume below!")
        
        listobs(vis=msName,selectdata=True,spw="",field="",antenna="",uvrange="",timerange="",correlation="",scan="",intent="",feed="",array="",observation="",verbose=True,listfile="",listunfl=False,cachesize=50,overwrite=False)
        shortresume_ms(msName)
        
        print("Split only the calibrators from MS " + msName)
        field = raw_input("Please enter the calibrator field: ")
        print("field = {0}".format(field))
        print("outputvis = " + msName + ".clb")
        print("datacolumn = 'data'")
        print("width = 1")
        print("Begin split...")
        split(vis=msName,outputvis=msName+".clb",keepmms=True,field=field,spw="",scan="",antenna="",correlation="",timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="data",keepflags=True,width=1,timebin="0s",combine="")
        print("End split.")




########### SELFCAL
def cal_selfcal():
    """Self-cal procedure"""
    print("\n......................................................................")
    print("Selfcal")
    msName = raw_input("Enter the MS name/location: ")

    shortresume_ms(msName)

    field = raw_input("Please enter the field you want to selfcal [default = 0]: ")
    if field == "": field = '0'

    phasecenter = raw_input("Please enter phase center: ")

    cell, imsize = estimate_cell_and_imsize(msName, field)

    nant, centerantenna = find_number_and_center_of_antenna(msName)
    print "Number of antenna: ", nant, "\nAntenna center: ", centerantenna
    plotants(vis=msName,figfile="")
    refant = raw_input("Choose reference antenna [default: {0}]: ".format(centerantenna))
    if refant == "": refant = centerantenna

    sensitivity = estimate_sensitivity(msName)

    cycle = 1
    selfcal = True
    default_threshold = [str(40*sensitivity)+"mJy", str(20*sensitivity)+"mJy", str(4*sensitivity)+"mJy", str(3*sensitivity)+"mJy"]
    default_solint = ['2min', '30s', '30s', '30s']
    default_niter = [500, 500, 1000, 1000]

    while selfcal:
        # run a cycle (clean, gaincal, applycal, split)
        print("-------------------------------------------------")
        print("Cycle: {0}".format(cycle))

        interactive = raw_input("Do you want to do interactive cleaning? (y or [n]): ")
        if interactive == "" or interactive.lower()[0] == "n":
            interactive = False
            threshold = raw_input("For this cycle, set threshold in CLEAN algorithm [default = {0}]: ".format(default_threshold[cycle-1]))
            if threshold == "": 
                threshold = default_threshold[cycle-1]
        else:
            interactive = True
            threshold = str(4.0 * sensitivity) + "mJy" 

        niter = default_niter[cycle-1]

        solint = raw_input("Set the value of solint (GAINCAL, default={0}): ".format(default_solint[cycle-1]))
        if solint == "": 
            solint = default_solint[cycle-1]

        solnorm = raw_input("Set the value of solnorm (in GAINCAL) [y/n, default=false]: " )
        if solnorm == "" or solnorm[0] == "n": 
            solnorm = False
        else: 
            solnorm = True
        
        calmode = raw_input("Choose calmode [p or ap]: ")

        selfcal_cycle(msName, cycle, field, niter, threshold, cell, imsize, phasecenter, interactive, solint, solnorm, refant, calmode)

        ans = raw_input("Continue next cycle of SELFCAL? ([y] or n): ")
        if ans == "" or ans.lower()[0] == 'y':
            cycle += 1
        else: 
            selfcal = False
            print("End of selfcal cycle")


    # last cleaning
    ans = raw_input("Do you want to do last cleaning? ([y] or n): ")
    if ans == "" or ans.lower()[0] == "y":
        print("Do final Cleaning...")
        pbcor = raw_input("Set pbcor (t or [f]): ")
        if pbcor == "" or pbcor.lower()[0] == "f":
            pbcor = False
        else:
            pbcor = True

        imgname = msName+".cont"+str(cycle+1)
        clean(vis=msName+".self"+str(cycle),imagename=imgname, outlierfile="",field="0",spw="",selectdata=True,
            timerange="",uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False,
            gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,
            psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=1000,
            gain=0.1,threshold=str(2.0 * sensitivity) + "mJy",psfmode="clark",imagermode="csclean",ftmachine="mosaic",
            mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,interactive=True,
            mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cell,phasecenter=phasecenter,
            restfreq="",stokes="I",weighting="briggs",robust=0.5,uvtaper=False,outertaper=[''],innertaper=['1.0'],
            modelimage="",restoringbeam=[''],pbcor=pbcor,minpb=0.2,usescratch=False,noise="1.0Jy",
            npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,
            flatnoise=True,allowchunk=False)

        # print 15 last line in clean log
        print_last15lines_fromlastlog()
        
        # print rms from residual image, imstat
        findrms_fromresidual(imgname)

    else:
        print("No final cleaning\n")


    print("END OF SELFCAL.\n")



########## CONCAT
def cal_concat():
    print("\n......................................................................")
    print("Concatenate several fields of an object")
    listvis = input('Enter the list of MS in array form e.g. ["lala.ms", "popo.ms"]: ')
    concatvis = raw_input("Enter the name of output MS: ")
    dirtol = raw_input("Dirtol [0.1arcsec]? ")
    if dirtol == "": 
        dirtol = '0.1arcsec'
    weightscale = input('Enter visweightscale in array form e.g. [1, 3.2, 0.5]: ')

    if len(listvis) == len(weightscale):
        print("Begin concat task ...")
        concat(vis=listvis,concatvis=concatvis,freqtol="",dirtol=dirtol,respectname=False,timesort=False,copypointing=False,visweightscale=weightscale,forcesingleephemfield="")
        print("End concat task.")
    else:
        print("[Error] Length of vis and weightscale is different.")



########## SUBSTRACT UVMODEL
def cal_subtractuvmodel():
    print("\n......................................................................")
    print("Subtract uvmodelfit from vis.")
    msName = raw_input("Enter the ms name/loc: ")

    # copy, for safety reason
    print("Duplicate MS -> .tosubtract.ms")
    outputms = msName+".tosubtract.ms"
    os.system("cp -r "+msName+" "+outputms)
    print "Copied to ", outputms

    print "Run uvmodelfit, [pointsource model, P]"
    uvmodelfit(vis=outputms,field="0",spw="",selectdata=True,timerange="",uvrange="",antenna="",scan="",
        msselect="",niter=5,comptype="P",sourcepar=[1.0, 0.0, 0.0],varypar=[],outfile="uvmodelfit."+msName+".complist")
    print "Done."

    
    # Do the Fourier Transform of the model into the MODEL column
    print "Run Fourier transfrom of model"
    ft(vis=outputms, complist="uvmodelfit."+msName+".complist")
    print "Done."

    # subtract :: CORRECTED-MODEL
    print "Subtracting..."
    ms.open(outputms, nomodify=False)
    ms.uvsub()
    ms.close()
    print "Done."

    # cleaning again
    print("split...")
    split(vis=outputms,outputvis=msName+".substracted.ms",keepmms=True,field=field,spw="",scan="",antenna="",correlation="",
        timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="corrected",keepflags=True,width=1,
        timebin="0s",combine="")
    print "Output: ", msName+".substracted.ms"

    ans = raw_input("Do you want to do cleaning of substracted MS? ([y] or n): ")
    if ans == '' or ans.lower()[0] == 'y':
        continuum(msName+".substracted.ms", field="", imagename=msName+".substracted.cont", niter=1000, threshold='0.0mJy', interactive=True)
    
    print "Finish substracting."


if __name__ == '__main__':
    #cal_split()
    #object_split()
    #cal_selfcal()
    #cal_subtractuvmodel()
    #cal_concat()
    estimate_cell_and_imsize("all_withweight_J0423-0120.ms", field='0')
    estimate_sensitivity("all_withweight_J0423-0120.ms")
    
    #continuum("uid___A002_X9a03fb_X2ca.ms.split.cal.clb.field_0.J0423-0120.ms", field='0', imagename="X2ca_field0_reimage.cont", phasecenter="J2000 04:23:15.800730 -01.20.33.06550" ,niter=5000)