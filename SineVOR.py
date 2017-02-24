# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 23:25:40 2014

@author: sonia_perez, julie_jung
"""
#packages to import
import pandas as pd
from pylab import *
from scipy import *
from scipy import optimize
import numpy

#variable is the location of the excel file to analyze

excelfilename = "/Users/juliejung/Desktop/Julie2015AugVOR.xlsx"

#get sheet names from excel file
sheetNames = pd.read_excel(excelfilename,sheetname='Sheet1')

#make empty numpy array, with sample name, rampdiff, r_r^2, lampdiff, l_r^2
results = []

#define a function which takes a dataframe (containing one tadpole's data)
#and returns a couple of key stats, in addition to saving out a picture of
#the tadpole's graph (w/ some data printed on there too)

def analyzeTadpole(tadpole,sampleName):
    myData = tadpole    
    num_points = len(myData["gamma'"])
    
    #set up some variables we'll need
    Tx = myData["gamma'"]
    Ty = Tx
    tX = myData['right']
    tY = myData['left']

    #fit function
    #fitfunc = lambda p, x: p[0]*sin(2*pi/p[1]*x+p[2]) + p[3] # Target function
    fitfunc = lambda p, x: p[0]*sin(2*pi/p[1]*(x-p[2])) + p[3] # Target function
    errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function
    
    
    #guess parameters
    guess_mean = np.mean(tX)
    guess_std = 3*np.std(tX)/(2**0.5)
    guess_phase = 0
    
    # Fit the right eye
    #p0 = [guess_std, 360, 0., guess_mean]
    #amp, wavelength, phase, something else
    p0 = [guess_std,360.,0, guess_mean] # Initial guess for the parameters
    p1,cov,infodict,mesg,ier= optimize.leastsq(errfunc, p0[:], args=(Tx, tX),full_output=True)

    # Right eye fitcurve R*2 calculations
    ss_err1=(infodict['fvec']**2).sum()
    ss_tot1=((tX-tX.mean())**2).sum()
    rsquared1=1-(ss_err1/ss_tot1)

    time = linspace(Tx.min(), Tx.max(), 100)
    # Plot of the data and the fit
    plot(Tx, tX, "ro", label = "right")
    plot(time,fitfunc(p1,time), "r-", label = "right fit")
    print numpy.polyfit(time, fitfunc(p1,time), 4) 
    #prints out coefficients of the multi-level polynomial! So the last number is the coefficient
    # 4 is the degrees that i've arbitrarily chosen- could choose 5 or soemthing else as well but we're only interested in the coefficient anyways
   
   
    # Fit the left eye
    #p0 = [guess_std2, 360, 0., guess_mean2]
    #amp, wavelength, phase, something else
    
    #guess parameters (left eye)
    guess_mean2 = np.mean(tY)
    guess_std2 = 3*np.std(tY)/(2**0.5)
    guess_phase2 = 0

    leftAmpGuess = abs(max(tY) - min(tY))
    
    p0 = [guess_std2,360., 0, guess_mean2] #Initial guess for the parameters
    p2,cov,infodict,mesg,ier= optimize.leastsq(errfunc, p0[:], args=(Tx, tY),full_output=True)

    # Left eye fitcurve R*2 calculations
    ss_err2=(infodict['fvec']**2).sum()
    ss_tot2=((tY-tY.mean())**2).sum()
    dft = num_points - 1
    dfe = num_points - 4
    rsquared2=1-((ss_err2/dfe)/(ss_tot2/dft))
    #print 'R^2 left =', rsquared2

    time = linspace(Ty.min(), Ty.max(), 100)
    #Plot of the data and the fit
    plot(Ty, tY, "b^", label = "left")
    plot(time,fitfunc(p2,time), "b-", label = "left fit")
    print numpy.polyfit(time, fitfunc(p2,time), 4) 
    #prints out coefficients of the multi-level polynomial! So the last number is the coefficient
    # 4 is the degrees that i've arbitrarily chosen- could choose 5 or soemthing else as well but we're only interested in the coefficient anyways
    
    #Legend the plot
    #title("Oscillations in the compressed trap")
    xlabel("Body angle (degrees)")
    ylabel("Eye angle (degrees)")
    legend(bbox_to_anchor = (1.05,1),loc = 2, borderaxespad = 0.)
    
    #limits of x and y axis
    ylim([-35,35])
    xlim([-200,200])
    
    ax = axes()

    #Get amplitude of fitfunc within time for 'right eye'
    rhighPoint = max(fitfunc(p1,time))
    rlowPoint = min(fitfunc(p1,time))

    rampDiff = abs(rhighPoint - rlowPoint)
 
    #Add label for the right eye fitcurve amplitude
    rmyText = "Right fit amplitude = %.2f" % rampDiff

    text(-160,25,rmyText, ha='left')

    #Get amplitude of fitfunc within time for 'left eye'
    lhighPoint = max(fitfunc(p2,time))
    llowPoint = min(fitfunc(p2,time))

    lampDiff = abs(lhighPoint - llowPoint)

    #Add label for the left eye fitcurve amplitude
    lmyText = "Left fit amplitude = %.2f" % lampDiff
    text(-160,17,lmyText, ha='left')
    #Add label for the R-squared values
    myText1 = "R-squared = %.3f" % rsquared1

    text(-160,21,myText1, ha='left')

    myText2 = "R-squared = %.3f" % rsquared2

    text(-160,13,myText2, ha='left')
    
    #save the figure
    savefig('%s.pdf' % sampleName, bbox_inches='tight')
    clf()
    return  sampleName, rampDiff, rsquared1, lampDiff, rsquared2,
    
    
# fot every tadple on the first page (list), add onto a table the name, RampDiff, 
    # R-R2, LampDiff, L-R2
for i in range(len(sheetNames)):
    print i, " ", sheetNames['Sheet1'][i]
    x = pd.read_excel(excelfilename, 
                     sheetname = sheetNames['Sheet1'][i])
    results.append(analyzeTadpole(x,str(sheetNames['Sheet1'][i])))
    
tastyResults = pd.DataFrame(results, columns = ["sample","RampDiff", "R-R2", 
                                                "LampDiff","L-R2"]) 
#                                                 index = sheetNames['sheetname'])

#location of the excel file to be created (contains results table)
resultsfile="/Users/juliejung/Desktop/VORindivsAug.xlsx"

#execute command to create excel file
tastyResults.to_excel(excel_writer = resultsfile)

 
    
###FOR ANALYZING A SINGLE TADPOLE###

#if our tadpole is C105T80

#mytadpole = pd.read_excel("/Users/sofia_perez/Desktop/C105T80.xlsx",
#                         sheetname = 'Sheet1')
   
#results.append(analyzeTadpole(mytadpole,'C105T80'))
#tastyResults = pd.DataFrame(results, columns = ["sample","RampDiff", "R-R2", "LampDiff","L-R2"]) 
                                               #    index = sheetNames['sheetname'])
#
#and then when you're done adding stuff#
#myData = pd.DataFrame(results, columns = ["sample", "RampDiff", "R-R2", 
  #                                                 "LampDiff","L-R2"])