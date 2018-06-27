import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('$DATADIR/analysis_scripts') #my directory of analysis_script # analysisUtils
import analysisUtils as aU

ARCSECTORAD = 4.84813681109536e-06


def find_image_parameter(casa_image):
    """find resolution, primary beam, freq, and center from CASA Image"""
    res = imhead(imagename=casa_image,mode="summary",hdkey="",hdvalue="",verbose=True)
    
    beam = res['restoringbeam']['major']['value'] # major axis
    print("Synthesized Beam: %f (arcsec)" % (beam))

    center = [res['refval'][0], res['refval'][1]] # RA Dec in radians
    print("Center (in rad): ", center)
    
    freq = res['refval'][3] # reference freq
    freq = freq/1e9 # in GHz
    primary_beam = aU.primaryBeamArcsec(frequency=freq, diameter=12.0, taper=10.0, obscuration=0.75, verbose=False, showEquation=True, use2007formula=True, fwhmfactor=None) # return primary beam in arcsec  
    print("Estimated primary beam size: %f (arcsec)" % (primary_beam))

    return center, beam, primary_beam


def calculate_rms_continuum(casa_image, n_sample, center_image, rmini, rmaxi, Rmin, Rmax):
    """Calculate RMS"""
    random.seed(1234)

    rms_array = []
    allrms = []

    # sampling
    for i in range(n_sample):
        # random radius of region
        r_region = random.uniform(Rmin, Rmax)

        rmin = rmini + r_region
        rmax = rmaxi - r_region

        u = random.random() # 0 <= u < 1
        v = random.random()
        # to get uniform distribution in Polar coord
        theta = 2*math.pi*u
        r = math.sqrt((rmax**2 - rmin**2)*v + rmin**2) # (ring)
        # in arcsec

        # change it back to cartesian
        r_rad = r*ARCSECTORAD
        dx = r_rad*math.cos(theta)
        dy = r_rad*math.sin(theta)
        x1 = center_image[0] + dx
        y1 = center_image[1] + dy
        x2 = center_image[0] - dx # symmetry
        y2 = center_image[1] - dy
        # print("Position: ", x1, y1, x2, y2, r, theta, r_region)


        regionrms1 = "circle[[" + str(x1) + "rad, " + str(y1) + "rad], " + str(r_region) + "arcsec" + "]"
        regionrms2 = "circle[[" + str(x2) + "rad, " + str(y2) + "rad], " + str(r_region) + "arcsec" + "]"

        stats1 = imstat(imagename=casa_image,axes=-1,region=regionrms1,box="",chans="",stokes="I",listit=True,verbose=True, 
            mask="",stretch=False,logfile='',append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,
            maxiter=-1,clmethod="auto")

        stats2 = imstat(imagename=casa_image,axes=-1,region=regionrms2,box="",chans="",stokes="I",listit=True,verbose=True, 
            mask="",stretch=False,logfile='',append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,
            maxiter=-1,clmethod="auto")

        # print("RMS: ", stats1['rms'][0], stats2['rms'][0])

        rms_pair = (stats1['rms'][0] + stats2['rms'][0])/2.
        allrms.append(rms_pair)
        rms_array.append([x1, y1, x2, y2, r, theta, r_region, stats1['rms'][0], stats2['rms'][0], rms_pair]) # numpy array
        


    # clipping
    allrms = np.array(allrms)
    stdrms = np.std(allrms)
    clip = 3*stdrms
    medianrms =np.median(allrms)

    __selectedrms = allrms[allrms < medianrms+clip]
    selectedrms = __selectedrms[__selectedrms > medianrms-clip]

    rms = np.median(selectedrms)

    return rms, rms_array


def plot_dist(rms_accepted, rms_array):
    arr = np.array(rms_array)
    r = arr[:,4]
    theta = arr[:,5]
    r_region =arr[:,6]
    rms = arr[:,9]

    plt.subplot(311)
    plt.plot(r, rms, 'r.')
    plt.axhline(y=rms_accepted, color='k', linestyle='-', linewidth=2)
    plt.xlabel("Radius from center")
    plt.ylabel("RMS")

    plt.subplot(312)
    plt.plot(theta%math.pi, rms, 'b.')
    plt.axhline(y=rms_accepted, color='k', linestyle='-', linewidth=2)
    plt.xlabel("Theta from center (folded)")
    plt.ylabel("RMS")

    plt.subplot(313)
    plt.plot(r_region, rms, 'm.')
    plt.axhline(y=rms_accepted, color='k', linestyle='-', linewidth=2)
    plt.xlabel("Radius of sample")
    plt.ylabel("RMS")

    plt.show()



def calc(casa_image, n_sample=100, plot=True):
    center, beam, primary_beam = find_image_parameter(casa_image)

    # size of region
    Rmin = max(primary_beam/12.0, 1.5*beam) # Dmin = max(PB/6, 3*beam)
    Rmax = primary_beam/6.0 # Dmax = PB/3

    # position of region
    rmin = 1.5*beam
    rmax = 0.75*primary_beam

    print("Calculate RMS from CASA image: "+casa_image)
    rms, rms_array = calculate_rms_continuum(casa_image, n_sample, center, rmin, rmax, Rmin, Rmax)

    # to check
    if plot:
        plot_dist(rms, rms_array)

    return rms, rms_array


if __name__ == '__main__':
    # casa_imagefile = "/mnt/sciops/data/rwibowo/projects/J1139-1350/Band3/uid___A002_Xbb5cfb_X2236.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.cont4.image"
    casa_imagefile = raw_input("CASA Image: ")
    rms, rms_array = calc(casa_imagefile)
    print(rms)

    # loc = '/mnt/sciops/data/rwibowo/projects/J1139-1350/Band3/'
    # imgfiles = ['uid___A002_Xbb5cfb_X2236.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.cont4.image', 'uid___A002_Xbb85b6_X2bec.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.cont4.image', 
    # 'uid___A002_Xbb85b6_X31a8.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.cont4.image', 'uid___A002_Xbbc4d2_X2616.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.cont4.image',
    # 'uid___A002_Xbbc4d2_X2928.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.cont4.image', 'uid___A002_Xbb5cfb_X2236.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.self3.substracted.cont.image', 
    # 'uid___A002_Xbb85b6_X2bec.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.self3.substracted.cont.image', 'uid___A002_Xbb85b6_X31a8.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.self3.substracted.cont.image',
    # 'uid___A002_Xbbc4d2_X2616.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.self3.substracted.cont.image', 'uid___A002_Xbbc4d2_X2928.ms.split.cal-CALIBRATE_PHASE-J1139-1350.ms.self3.substracted.cont.image']

    # for img in imgfiles:
    #     imgname = loc + img
    #     rms, rms_array = calc(imgname, plot=False)
    #     print("Result: ", img, rms, "\n\n")
