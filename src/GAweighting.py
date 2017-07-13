# script for weighting problem, with GA

import os
import sys
import glob
import time
import random
import numpy as np
import matplotlib.pyplot as plt

from deap import base, creator, tools, algorithms # import DEAP



def function_to_evaluate(weightlist):
    '''Function that we want to optimize: calculate rms'''
    mslist         = ["uid___A002_X9646fb_Xccd.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_X9646fb_Xccd.ms.split.cal.clb.field_1.J0423-0120.ms.self3.substracted.ms", "uid___A002_X97db9c_X159b.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_X95b353_X6f7.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_X9630c0_X81b.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_Xa916fc_X6d6b.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_Xa916fc_X6d6b.ms.split.cal.clb.field_1.J0423-0120.ms.self3.substracted.ms", "uid___A002_Xaaa05f_X15bb.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_Xaaa05f_X1af7.ms.split.cal.clb.field_0.J0423-0120.ms.self3.substracted.ms", "uid___A002_Xaaa05f_X1af7.ms.split.cal.clb.field_1.J0423-0120.ms.self3.substracted.ms"]
    concatresultMS = "concatresult.ms"
    imgname        = "concatresult.ms.cont"
    niter          = 500
    threshold      = "0.0mJy"
    psfmode        = 'clark'
    interactive    = False 
    mask           = "circle[[65.8158364deg, -1.3425182deg], 5arcsec]"
    imsize         = 1250
    cell           = "0.04arcsec"
    phasecenter    = "J2000 04:23:15.800730 -01.20.33.065501"
    weighting      = 'briggs'
    robust         = 0.5
    pbcor          = False
    regionrms      = "annulus[[65.8158364deg, -1.3425182deg], [5arcsec, 15arcsec]]" # region to calculate RMS
    weightlist     = [1] + weightlist # append, 1st MS always use "1" as weightscale
    

    # CONCAT
    if os.path.exists(concatresultMS):
        # print "Removing previous concat MS..."
        os.system('rm -rf '+concatresultMS)

    concat(vis=mslist,concatvis=concatresultMS,freqtol="",dirtol="0.1arcsec",
        respectname=False,timesort=False,copypointing=False,
        visweightscale=weightlist,forcesingleephemfield="")


    # CLEAN
    # if image already exist
    for ifile in glob.glob(imgname+".*"): # for safety reason
        # print "Removing previous image file:", ifile
        os.system("rm -rf "+ifile)

    # print "CLEANing..."
    clean(vis=concatresultMS,imagename=imgname,outlierfile="",field="0",spw="",selectdata=True,
        timerange="",uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",
        wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,
        wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=niter,gain=0.1,threshold=threshold,
        psfmode=psfmode,imagermode="csclean",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],
        negcomponent=-1,smallscalebias=0.6,interactive=interactive,mask=mask,nchan=-1,start=0,width=1,outframe="",
        veltype="radio",imsize=imsize,cell=cell,phasecenter=phasecenter,restfreq="",stokes="I",weighting=weighting,
        robust=robust,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=pbcor,
        minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,
        reffreq="",chaniter=False,flatnoise=True,allowchunk=False)


    # Calculate RMS
    # print "Calculate RMS..."
    stats = imstat(imagename=imgname+".image",axes=-1,region=regionrms,box="",chans="",stokes="I",listit=True,verbose=True,
        mask="",stretch=False,logfile='',append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,
        maxiter=-1,clmethod="auto")

    rms = stats['rms'][0]
    
    return rms,


numberofMS = 10 # 1st MS not included, so each individual has 9 parameters.
minweight = 0.0
maxweight = 4.0
indpb = 0.1


creator.create("FitnessMin", base.Fitness, weights=(-1.0,)) # minimize
creator.create("Individual", list, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

# function random.random as input 
# uniform random 0 to 4 as 'individu generator'
toolbox.register("attr_bool", random.uniform, minweight, maxweight) 
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=numberofMS-1) 
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", function_to_evaluate) # insert the function name here.
toolbox.register("mate", tools.cxTwoPoint) # how to do crossover
toolbox.register("mutate", tools.mutFlipBit, indpb=indpb) # probability of individu get a mutation
toolbox.register("select", tools.selTournament, tournsize=3) # tournament type "selection" before "mating process"


def main(number_of_ind=39, number_of_generation=10, cxpb=0.5, mutpb=0.2, verbose=False):
    # generate 0th population
    pop = toolbox.population(n=number_of_ind)
    # insert a "woho" -> best guess
    guess1 = creator.Individual([0.55, 0.82, 0.34, 0.28, 0.80, 0.37, 0.08, 0.05, 0.05])
    guess2 = creator.Individual([1, 1, 1, 1, 1, 1, 1, 1, 1])
    guess3 = creator.Individual([1, 1, 1, 1, 1, 1, 0, 0, 0])
    
    pop.append(guess1)
    pop.append(guess2)
    pop.append(guess3)
    # real number of individual in population is number_of_ind + number of guess

    # calculate fitness for the 0th generation
    invalids = [ind for ind in pop if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalids) # parallelization on solving the "objective function"
    for ind, fit in zip(invalids, fitnesses):
        ind.fitness.values = fit
    
    # initiate logging
    # Hall of Fame
    hof = tools.HallOfFame(1) # save the best of the best, only 1 individu
    # Statistics
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("min", np.min)
    stats.register("max", np.max)
    # Log book
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + stats.fields
    
    # logging
    record = stats.compile(pop)
    logbook.record(gen=0, nevals=len(invalids), **record)
    if verbose:
        print logbook.stream

    total_evals = len(invalids)
    best = []
    best.append(tools.selBest(pop, k=1)[0])
        
    for gen in range(1, number_of_generation+1): # iteration
        # Select the next generation individuals
        offspring = toolbox.select(pop, k=len(pop))
        
        # Vary the pool of individuals
        # inside: crossover AND mutation
        offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

        # some individual are "new" offspring, others are still the same parent 
        # e.g. not selected in crossover and also not a mutan
        # if fitness is "invalid", not calculated yet (real offspring) -> calculate fitness..
        invalids = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalids)
        for ind, fit in zip(invalids, fitnesses):
            ind.fitness.values = fit
        
        # Replace the current population by the offspring
        pop[:] = offspring
        
        # Update the hall of fame with the generated individuals
        hof.update(pop)

        # logging
        record = stats.compile(pop)
        logbook.record(gen=gen, nevals=len(invalids), **record)
        total_evals += len(invalids)
        if verbose:
            print logbook.stream
    
        best.append(tools.selBest(pop, k=1)[0])
    
    return pop, logbook, hof, best, total_evals


start_time = time.time()
pop, log, hof, best, total_evals = main(17, 50, cxpb=0.5, mutpb=0.1)
print("Best individual is: %s\nwith fitness: %s" % (hof[0], hof[0].fitness))
print "Total evaluation of objective function: ", total_evals
print "Running time: ", (time.time() - start_time)/3600.0, " hours"

logtxt = log.select("gen", "nevals", "avg", "min", "max")
np.savetxt("log3.txt", logtxt, delimiter=' ')
np.savetxt("onlythebest3.txt", best, delimiter=' ')
np.savetxt("bob3.txt", hof[0] + [hof[0].fitness.getValues()[0]], delimiter=' ')

gen, avg, min_, max_ = log.select("gen", "avg", "min", "max")
plt.plot(gen, avg, label="average")
plt.plot(gen, min_, label="minimum")
plt.plot(gen, max_, label="maximum")
plt.xlabel("Generation")
plt.ylabel("Fitness (rms)")
plt.legend(loc="upper right")
plt.show()


# b = np.array(best)
# pars = b.T
# iteration = np.arange(len(pars[0]))


# fig, axs = plt.subplots(3,3, figsize=(13,10), sharex=True, sharey=True)
# axs = axs.flat

# for i, par in enumerate(pars):
#     axs[i].plot(iteration, par, 'k.-', alpha=0.3)
#     axs[i].set_ylim([-0.1,4.1])
#     axs[i].set_title(str(i))
#     fig.tight_layout()

