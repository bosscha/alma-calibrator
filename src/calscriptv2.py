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

class selfCal:
    """A class to do selfcal semi-automatic"""
    # automatic naming
    # automatic cell size
    # give a suggestion on imsize, threshold, antenna center

    def __init__(self, selfcalpar=""):
        if selfcalpar != "": # given filename
            pars = self.__read_selfcalinput(selfcalpar)
            self.msName       = pars[0]
            self.field        = pars[1]
            self.interactive  = pars[2]
            self.phasecenter  = pars[3] 
            self.threshold    = pars[4]
            self.niter        = pars[5]
            self.pbcor        = pars[6]
            self.calmode      = pars[7]
            self.solint       = pars[8]
            self.solnorm      = pars[9]
            self.autorefant   = pars[10]
            self.startcycle   = pars[11]

            print "Parameter for selfcal:"
            print "------------------------------------------------"
            for i in pars[:12]: 
                print i
            print "------------------------------------------------"

    def __read_selfcalinput(self, filename):
        with open(filename,"r") as ifile:
            for iline in ifile:
                line = iline.strip().split("=")

                if len(line) > 1:
                    var, value = line[0].strip(), line[1].strip()
                    
                    if var == 'MSNAME': 
                        msName = value.split()

                    if var == 'FIELD': 
                        field = value

                    if var == 'INTERACTIVE': 
                        interactive = value.split()
                        for i,par in enumerate(interactive):
                            if par == 'T':
                                interactive[i] = True
                            else:
                                interactive[i] = False

                    if var == 'PHASECENTER': 
                        phasecenter = value

                    if var == 'THRESHOLD': 
                        threshold = value.split()
                        for i,par in enumerate(threshold):
                            threshold[i] = float(par)

                    if var == 'NITER': 
                        niter = value.split()
                        for i,par in enumerate(niter):
                            niter[i] = int(par)

                    if var == 'PBCOR': 
                        pbcor = value.split()
                        for i,par in enumerate(pbcor):
                            if par == 'T':
                                pbcor[i] = True
                            else:
                                pbcor[i] = False

                    if var == 'CALMODE': 
                        calmode = value.split()

                    if var == 'SOLINT': 
                        solint = value.split()

                    if var == 'SOLNORM': 
                        solnorm = value.split()
                        for i,par in enumerate(solnorm):
                            if par == 'T':
                                solnorm[i] = True
                            else:
                                solnorm[i] = False

                    if var ==  'AUTOREFANT': 
                        autorefant = value
                        if autorefant == 'T':
                            autorefant = True
                        else:
                            autorefant = False

                    if var == 'STARTCYCLE':
                        startcycle = int(value)



        return msName, field, interactive, phasecenter, threshold, niter, pbcor, calmode, solint, solnorm, autorefant, startcycle


    def is_power2(self, num): # num must be integer
        return (num & (num-1)) == 0


    def estimate_freq(self, msName):
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


    def estimateSynthesizedBeam_perfield(self, msName, freq, field):
        """useL80method, modified from analysisUtils"""
        result = aU.getBaselineStats(msName, field=field, percentile=80) # field can be more than one... -_- choose field field with larger baseline?
        L80 = result[0] # meters
        arcsec = 3600*np.degrees(1)*0.574*aU.c_mks/(self.estimate_freq(msName)*L80)

        return arcsec


    def estimateSynthesizedBeam(self, msName, freq, field=''):
        if len(field) > 1: # if there is more than one field (combine)
            fields = field.split(",")
            synthesizedbeam = []
            for fi in fields:
                synthesizedbeam.append(self.estimateSynthesizedBeam_perfield(msName, freq, fi))
            arcsec = min(synthesizedbeam) # choose the minimum one
        else:
            arcsec = self.estimateSynthesizedBeam_perfield(msName, freq, field)

        return arcsec


    def estimate_cell_and_imsize(self, msName, field):
        print "Calculate the proper cell and imsize"

        ##### cell
        freq = self.estimate_freq(msName)
        print "Estimated Freq: %s GHz" % str(freq/1e9)

        spatial_resolution = self.estimateSynthesizedBeam(msName, freq, field=field)
        print "Estimated spatial resolution: %f arcsec" % (spatial_resolution)
        
        cell = 0.2*spatial_resolution # set cell to spatRes/5
        cell = '%.2f' % cell + 'arcsec' 
        print("Select cell size: {0}".format(cell))

        ### imsize
        primary_beam = aU.primaryBeamArcsec(frequency=freq, diameter=12.0, \
            taper=10.0, obscuration=0.75, verbose=False, showEquation=True, \
            use2007formula=True, fwhmfactor=None) # return primary beam in arcsec  
        print "Estimated primary beam size: %f" % (primary_beam)

        #imsize_onehalf = int(1.5*primary_beam/(0.2*spatial_resolution))
        imsize_oneseven = int(1.7*primary_beam/(0.2*spatial_resolution))
        imsize_two = int(2.0*primary_beam/(0.2*spatial_resolution))

        imsize = imsize_oneseven
        unfactorize = True
        while unfactorize:
            # casa want even number and factorizeable of 2,3,5,7 (?) -> 2^n x 10; sometime it will be fail
            if self.is_power2(imsize) or (imsize % 10 == 0 and self.is_power2(imsize/10)) or (imsize % 100 == 0 and (((imsize/100) % 2 == 0) or ((imsize/100) % 5 == 0))): 
                 unfactorize = False
            else: imsize += 1

        imz = raw_input("The primary beam multiple -> 1.7x = {0}, 2x = {1}, suggested imsize (default) = [{2}]. Enter your imsize: ".format(imsize_oneseven, imsize_two, imsize))
        if imz == "": imz = imsize
        else: imz = int(imz)

        print("Select imsize: {0}".format(imz))

        return cell, imz


    def shortresume_ms(self, msName):
        intentSources = es.getIntentsAndSourceNames(msName)
        for key in intentSources:
            intent = key
            fieldid = intentSources[key]['id']
            names = intentSources[key]['name']
            if len(fieldid) > 0:
                print "Intent: ", intent, "               Field id: ", fieldid, "                  Name: ", names


    def find_number_and_center_of_antenna(self, msName):
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


    def get_spwedge(self, msName):
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


    def get_npolarizations(self, msName):
        tb.open(msName)
        nPolarizations = np.shape(tb.getcell("DATA",0))[0]
        tb.close()
        return nPolarizations


    def estimate_sensitivity(self, msName):
        # find integration time
        time_os = aU.timeOnSource(msName)
        print("\n-------------------------------------------------")
        int_time = time_os[0]['minutes_on_source'] # always in field 0 (from split)
        print "Integration time = ", int_time, " min"

        nant, centerantenna = self.find_number_and_center_of_antenna(msName) # find number of antenna
        npol = self.get_npolarizations(msName) # number of polarizations
        spwedge = self.get_spwedge(msName)

        bandunion = [] # union from that some intervals
        for begin,end in sorted(spwedge): # sort spwedge base on lower limit of spw
            if bandunion and bandunion[-1][1] >= begin: # check the intersect
                bandunion[-1][1] = max(bandunion[-1][1], end) # combine
            else:
                bandunion.append([begin, end]) #seperate interval
        
        print "Merged spw: ", bandunion
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
            tb.open(msName + "/ASDM_CALATMOSPHERE") # sometime no table found
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
            meanTsys = tsys_nominal[bandused]
            print "Using standard Tsys for this band: ", meanTsys

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


    def print_last15lines_fromlastlog(self):
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


    def findrms_fromresidual(self, imgname):
        residualimg = imgname + ".residual"
        xstat = imstat(residualimg)
        print "RMS from residual-image " + residualimg+ " (using imstat): ", xstat['rms'][0]
        return xstat['rms'][0]


    def selfcal_cycle(self, msName, cycle, field, niter, threshold, cell, imsize, phasecenter, interactive, pbcor, solint, solnorm, refant, calmode):
        #clean
        if cycle == 1:
            vis = msName
            field = field
            if os.path.exists(msName + ".flagversions"):
                os.system("rm -rf "+msName + ".flagversions")
        else:
            vis = msName+'.self'+str(cycle-1)
            field = '0'

        print("-------------------------------------------------")
        print("Cycle: {0}".format(cycle))

        imgname = msName+".cont"+str(cycle)

        # check if image's exist already for this cycle
        for ifile in glob.glob(imgname+".*"):
            print "Removing previous image file:", ifile
            os.system("rm -rf "+ifile) # remove


        print("Cleaning...")
        clean(vis=vis,imagename=imgname,outlierfile="",field=field,spw="",selectdata=True,timerange="",uvrange="",antenna="",
            scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",
            rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",
            interpolation="linear",niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",
            ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,
            interactive=interactive,mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",
            imsize=imsize,cell=cell,phasecenter=phasecenter,restfreq="",stokes="I",
            weighting="briggs",robust=0.5,uvtaper=False,outertaper=[''],innertaper=['1.0'],
            modelimage="",restoringbeam=[''],pbcor=pbcor,minpb=0.2,usescratch=False,noise="1.0Jy",
            npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,
            flatnoise=True,allowchunk=False)

        # print 15 last line in clean log
        self.print_last15lines_fromlastlog()

        # print rms from residual image, imstat
        self.findrms_fromresidual(imgname)

        #gaincal
        caltab = msName+".G"+str(cycle)
        # check if gaintable exist already
        if os.path.exists(caltab):
            print "Removing gaintable:", caltab
            os.system("rm -rf "+caltab)

        print("gaincal...")
        gaincal(vis=vis,caltable=caltab,field=field,spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",
            scan="",observation="",msselect="",solint=solint,combine="",preavg=-1.0,refant=refant,
            minblperant=4,minsnr=3.0,solnorm=solnorm,gaintype="G",smodel=[],calmode=calmode,append=False,splinetime=3600.0,
            npointaver=3,phasewrap=180.0,docallib=False,callib="",gaintable=[],gainfield=[],interp=[],spwmap=[],parang=False)


        #check gain phase
        #print("plotcal... check")
        #plotcal(caltable=msName+".G"+str(cycle),xaxis="time",yaxis="phase",poln="",field=field,antenna="",spw="",timerange="",
        #    subplot=331,overplot=False,clearpanel="Auto",iteration="antenna",plotrange=[],showflags=False,plotsymbol="o",
        #    plotcolor="blue",markersize=5.0,fontsize=10.0,showgui=False,figfile="plotcal_"+str(cycle)+".png")
        
        outputv = msName+".self"+str(cycle)
        # check if split result already exists
        if os.path.exists(outputv):
            print "Removing ms: ", outputv
            os.system("rm -rf "+outputv)
            print "Removing .flagversions file"
            os.system("rm -rf "+outputv+".flagversions")

        #applycal
        print("applycal...")
        applycal(vis=vis,field=field,spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",
            msselect="",docallib=False,callib="",gaintable=caltab,gainfield=[],
            interp="nearest",spwmap=[],calwt=[True],parang=False,applymode="",flagbackup=True)

        #split
        print("split...")
        split(vis=vis,outputvis=outputv,keepmms=True,field=field,spw="",scan="",antenna="",correlation="",
            timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="corrected",keepflags=True,width=1,
            timebin="0s",combine="")

        print("Finish Cycle {0} \n".format(cycle))



    def cal_selfcal(self):
        """Self-cal procedure"""
        print("\n......................................................................")
        print("Selfcal")
        print("......................................................................")
        for msname in self.msName:
            print msname
            self.shortresume_ms(msname)

            field = self.field
            phasecenter = self.phasecenter
            cell, imsize = self.estimate_cell_and_imsize(msname, field)
            nant, centerantenna = self.find_number_and_center_of_antenna(msname)
            print "Number of antenna: ", nant, "\nAntenna center: ", centerantenna
            if self.autorefant:
                refant = centerantenna
            else:
                plotants(vis=msname, figfile="")
                refant = raw_input("Choose reference antenna [default: {0}]: ".format(centerantenna))
                if refant == "": 
                    refant = centerantenna

            sensitivity = self.estimate_sensitivity(msname)
            
            default_threshold = []
            for multiple in self.threshold:
                default_threshold.append(str(multiple * sensitivity) + "mJy")

            for cycle in range(self.startcycle, len(self.calmode)+1): # start from 1
                interactive = self.interactive[cycle-1]
                threshold = default_threshold[cycle-1]
                niter = self.niter[cycle-1]
                pbcor = self.pbcor[cycle-1]
                solint = self.solint[cycle-1]
                solnorm = self.solnorm[cycle-1]
                calmode = self.calmode[cycle-1]

                # run a cycle (clean, gaincal, applycal, split)
                self.selfcal_cycle(msname, cycle, field, niter, threshold, cell, imsize, phasecenter, interactive, pbcor, solint, solnorm, refant, calmode)


            # LAST CLEANING
            print "Last Cleaning"
            print "----------------------------------------------"
            imgname = msname+".cont"+str(cycle+1)

            # check if image's exist already for this cycle
            for ifile in glob.glob(imgname+".*"):
                print "Removing previous image file:", ifile
                os.system("rm -rf "+ifile) # remove
            
            clean(vis=msname+".self"+str(cycle),imagename=imgname, outlierfile="",field="0",spw="",selectdata=True,
                timerange="",uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False,
                gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,
                psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=self.niter[-1],
                gain=0.1,threshold=str(self.threshold[-1] * sensitivity) + "mJy",psfmode="clark",imagermode="csclean",ftmachine="mosaic",
                mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,interactive=self.interactive[-1],
                mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cell,phasecenter=phasecenter,
                restfreq="",stokes="I",weighting="briggs",robust=0.5,uvtaper=False,outertaper=[''],innertaper=['1.0'],
                modelimage="",restoringbeam=[''],pbcor=self.pbcor[-1],minpb=0.2,usescratch=False,noise="1.0Jy",
                npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,
                flatnoise=True,allowchunk=False)

            # print 15 last line in clean log
            self.print_last15lines_fromlastlog()
            
            # print rms from residual image, imstat
            self.findrms_fromresidual(imgname)

            print("END OF SELFCAL.\n")




    def continuum(self, msName, field='', imagename='', phasecenter='', niter=1000, threshold='0.0mJy', interactive=True, pbcor=False):
        """imaging the continuum"""
        
        if imagename == '':
            imagename = msName + '.cont' # default is the msName.cont

        cell, imsize = self.estimate_cell_and_imsize(msName, field)
        sensitivity = self.estimate_sensitivity(msName)
        if threshold == '': threshold = raw_input("Enter your threshold [0.0mJy]: ")
        if phasecenter == '': phasecenter = raw_input("Enter your phasecenter if you want ['']: ")


        print("Begin cleaning task...")
        clean(vis=msName,imagename=imagename, outlierfile="",field=field,spw="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",\
            intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,\
            painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",\
            niter=niter,gain=0.1,threshold=threshold,psfmode="clark",imagermode="csclean",ftmachine="mosaic",\
            mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,interactive=interactive,\
            mask=[],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cell,\
            phasecenter=phasecenter,restfreq="",stokes="I",weighting="briggs",robust=0.5,uvtaper=False,outertaper=[''],\
            innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=pbcor,minpb=0.2,usescratch=False,noise="1.0Jy",\
            npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,\
            allowchunk=False)
        print("End cleaning task.")

        # print 15 last line in clean log
        self.print_last15lines_fromlastlog()

        # print rms from residual image, imstat
        self.findrms_fromresidual(imagename)



    def cal_subtractuvmodel(self):
        print("\n......................................................................")
        print("Subtract uvmodelfit from vis.")
        print("\n......................................................................")

        for msname in self.msName:
            # copy, for safety reason
            msname = msname + '.self3'
            print "Substracting ", msname
            print("Duplicate MS")
            outputms = msname+".tosubtract.ms"
            if os.path.exists(outputms):
                print "Remove " + outputms
                os.system("rm -r " + outputms)
            
            os.system("cp -r "+msname+" "+outputms)
            print "Copied to ", outputms

            ofile = "uvmodelfit."+msname+".complist"
            if os.path.exists(ofile):
                print "Remove " + ofile
                os.system("rm -r " + ofile)
            print "Run uvmodelfit, [pointsource model, P]"
            uvmodelfit(vis=outputms,field="0",spw="",selectdata=True,timerange="",uvrange="",antenna="",scan="",
                msselect="",niter=5,comptype="P",sourcepar=[1.0, 0.0, 0.0],varypar=[],outfile=ofile)
            print "Done."

            
            # Do the Fourier Transform of the model into the MODEL column
            print "Run Fourier transform of model"
            ft(vis=outputms, complist=ofile)
            print "Done."

            # subtract :: CORRECTED-MODEL
            print "Subtracting..."
            ms.open(outputms, nomodify=False)
            ms.uvsub()
            ms.close()
            print "Done."

            # cleaning again
            print("Begin split...")
            if os.path.exists(msname+".substracted.ms"):
                print "Remove " + msname+".substracted.ms"
                os.system("rm -rf " + msname+".substracted.ms")

            split(vis=outputms,outputvis=msname+".substracted.ms",keepmms=True,field='0',spw="",scan="",antenna="",correlation="",
                timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="corrected",keepflags=True,width=1,
                timebin="0s",combine="")
            print "End split, output: ", msname+".substracted.ms"

            # remove temporary MS
            if os.path.exists(outputms):
                print "Remove " + outputms
                os.system("rm -r " + outputms)

            #ans = raw_input("Do you want to do cleaning of substracted MS? ([y] or n): ")
            #if ans == '' or ans.lower()[0] == 'y':
            self.continuum(msname+".substracted.ms", field="", imagename=msname+".substracted.cont", phasecenter=self.phasecenter, niter=1000, threshold='0.0mJy', interactive=True)
            
            print "Finish substracting."


if __name__ == '__main__':
    selfc = selfCal("selfcal.par")
    selfc.cal_selfcal()
    #selfc.cal_subtractuvmodel()
