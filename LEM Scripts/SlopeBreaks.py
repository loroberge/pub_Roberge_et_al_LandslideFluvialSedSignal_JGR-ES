# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 21:35:21 2021

@author: benja
"""
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import pwlf, scipy

plt.close('all')

matplotlib.rcParams["font.size"] = 10
matplotlib.rcParams["pdf.fonttype"] = 42

# # Optionally set font to Computer Modern to avoid common missing font errors
# matplotlib.rc('font', family='serif', serif='cm10')
# matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

#%% Slope breakpoints
def get_slopeBreaks(x,y, checkS=True, checkR2=False, plot=False, verbose=False, plot_Name='River', printBreak=True):
    '''
    Function to calcualte slope breaks in slope area plots

    Parameters
    ----------
    x : drainage area 
    y : slope
    checkS : check if slopes of two identified segments differ significantly. The default is True.
    checkR2 :  check if R2 of originql datacloud versus piecewise regression differ significantly. The default is False.
    plot : TYPE, optional
        plot data. The default is False.
    verbose : TYPE, optional
        print output. The default is False.
    plot_Name : STRING, optional
        name of the plot. The default is 'River'.

    Returns
    -------
    my_pwlf : instance of pwlf
        instance of pwlf.
    breaks : numpy array
        breaks of the piecewise regression.

    '''
    my_pwlf = pwlf.PiecewiseLinFit(x,y)
    breaks = my_pwlf.fit(2)
    s = abs(my_pwlf.slopes)
    
    x_p = breaks#np.linspace(breaks[0], breaks[2],20)
    y_p = my_pwlf.predict(x_p)

        
    #Evaluate performance: part 1 
    if verbose: 
        aa = (my_pwlf.predict(x[x<breaks[1]]))
        bb =(y[x<breaks[1]])
        res = scipy.stats.linregress(aa,bb)
        print(f"R-squared part 1: {res.rvalue**2:.6f}")   
        print("RMSE part 1: " + str(np.sqrt(np.mean(np.power(aa-bb,2)))))
        print("Slope part 1: " + str(s[0]))

        #Evaluate performance: part 2
        aa= my_pwlf.predict(x[x>=breaks[1]])
        bb= y[x>=breaks[1]]   
        res = scipy.stats.linregress(aa,bb)
        print(f"R-squared part 2: {res.rvalue**2:.6f}")   
        print("RMSE part 2: " + str(np.sqrt(np.mean(np.power(aa-bb,2)))))
        print("Slope part 2: " + str(s[1]))
    
    if checkS:
        if (s[0]>0.9*s[1] and s[0]<1.1*s[1]):
            breaks =None
            if verbose: 
                print('No intercept found')
        if s[1]<s[0]:
            breaks =None
            if verbose: 
                print('High drainage results in lower slope, indicating lack of breakpoint and fluvial dynamics')
                
                
    # Check increase in R2 by piecewise
    if checkR2:
        res = scipy.stats.linregress(x,y) 
        r2_1 = res.rvalue**2
        r2_2 = my_pwlf.r_squared()
        if verbose:         
            print(f"R-squared original: {r2_1:.6f}")  
            print(f"R-squared piecewise: {r2_2:.6f}") 
        if (r2_2<r2_1 or r2_2<1.02*r2_1):
            breaks =None
            if verbose: 
                print('No intercept found') 
    
    if plot:
        plt.figure()
        plt.plot(x,y, "o")
        plt.plot(x_p, y_p)
        # plt.xlabel('Drainage area (' + r'$^{10}log m^2$)',fontweight='bold')
        plt.xlabel('$^{10}log (Drainage area, m^2)$')
        plt.ylabel('$^{10}log (slope, m/m)$')
        ax = plt.gca()
        if np.any(breaks) and printBreak:
            plt.text(0.8, 0.85,f'Breakpoint: {breaks[1]:.2f}', 
            horizontalalignment='center',
            verticalalignment='center', 
            transform = ax.transAxes)
        
        if plot_Name:
            plt.title(plot_Name)
        # plt.show()
        
    return my_pwlf,breaks,x_p,y_p
