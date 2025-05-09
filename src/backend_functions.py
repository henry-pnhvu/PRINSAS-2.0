# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 12:56:22 2024

@author: NHUHA
"""

import csv
import matplotlib
import numpy as np
import scipy.optimize as sci_opt
import plot_formating as pf


# Function used to subtract flat background and limit the background-subtracted 
# data to Q-max
def subtract_background(QQ_original, IQ_original, dIQ_original, bkgrd, QQ_min, QQ_max):
    # Subtract background and find negative IQ values
    IQ_subtract = IQ_original - bkgrd
    idx_IQ_negative = np.where(IQ_subtract <= 0)
    dIQ_percent = dIQ_original/IQ_original
    
    # Remove negative IQ values
    try:
        QQ_subtract = QQ_original[:np.min(idx_IQ_negative)]
        IQ_subtract = IQ_subtract[:np.min(idx_IQ_negative)]
        dIQ_subtract = dIQ_percent[:np.min(idx_IQ_negative)]
    except ValueError:
        QQ_subtract = QQ_original
        dIQ_subtract = dIQ_percent
    
    # Remove data points at Q values outside the specified Q range and return
    select_pos = np.logical_and(QQ_min < QQ_subtract, QQ_subtract < QQ_max) 
    QQ_trim = QQ_subtract[select_pos]
    IQ_trim = IQ_subtract[select_pos]
    dIQ_trim = dIQ_subtract[select_pos]
    return QQ_trim, IQ_trim, dIQ_trim*IQ_trim


# Function used to read SAS data file and return the corresponding Q, IQ and 
# dIQ values
def read_SANS_data(dir_data):
    # Read file content
    with open(dir_data) as file:
        data_raw = file.readlines()
    
    # Extract all the potential dilimiters of the file
    delim_list = []
    for line in data_raw:
        try:
            delim_list.append(csv.Sniffer().sniff(line).delimiter)
        except csv.Error: pass
    
    # Ensure delim_list is not empty
    if not delim_list:
        raise ValueError("No valid delimiters found in the data.")
   
    try:
        # Identify the most appropriate delimiter from the list of obtained delimiters
        delim, count = np.unique(delim_list, return_counts=True)
        select_delim = delim[count.argmax()]
        
        # Obtained data entries using the identified delimiters
        data_splitted = []
        for line in data_raw:
            entry_selected = []
            for entry in line.split(select_delim):
                try: 
                    if float(entry) >= 0:
                        entry_selected.append(float(entry))
                except: pass 
            data_splitted.append(entry_selected)
        line_len = [len(line) for line in data_splitted]
        num_col = max(set(line_len), key=line_len.count)
        data_selected = np.array([line for line in data_splitted 
                                  if len(line) == num_col])
    except ValueError:
        raise ValueError('Cannot read file content')
        
    # Check whether return result contains any data
    if data_selected.size == 0:
        raise ValueError('Cannot read file content')
        
    # Assigned and return read data to Q, IQ, and dIQ
    if data_selected.shape[1] == 2:
        return data_selected[:,0], data_selected[:,1], 0
    elif data_selected.shape[1] > 2:        
        return data_selected[:,0], data_selected[:,1], data_selected[:,2]


# This function clear the plotted data in the figure while retain all the axes
# and legend settings. This is done so that only data points are replotted when
# required, saving execution time.
def clear_plot(fig, canvas):
    # Get axes from the figure
    ax = fig.get_axes()[0]
    
    # To completely remove errorbar elements, identify LineCollection objects
    # in the list of artists, then delete 4 artist, starting from i-1 to i+2,
    # where i is the idx of LineCollection objects.
    # To remove legend of error bar plot, remove it from ax.containers
    artist_list = ax.get_children()
    for i in range(len(artist_list)):
        if isinstance(artist_list[i], matplotlib.collections.LineCollection):
            [artist_list[j].remove() for j in range(i-1, i+3)]
            ax.containers.remove(ax.containers[0])
            
    # Clear the remaining unlabeled data set
    for line in ax.get_lines():
        line.remove()

    # Clear written text
    for text in ax.texts:
        text.remove()

    # Rescale axes
    ax.relim()
    ax.autoscale()
    
    # Replot the top axis for SAS plot
    if len(fig.get_axes()) > 1:
        ax2 = fig.get_axes()[1]
        ax2.set_xlim(0.1/np.array(ax.get_xlim()))
        ax2.set_xticklabels(['{:g}'.format(2.5*val)
                             for val in ax2.get_xticks()])
    try:
        ax.get_legend().remove()
    except: pass
    canvas.draw() # Update plotting canvas on the program interface
        
    
# This function used to plot the original SAS data to the 'SAS Data' plotting
# window
def plot_SANS_data(QQ, IQ, dIQ, fig, canvas):
    ax = fig.get_axes()[0]
    
    # Removing previously plotted SANS data, this is done since input data are 
    # replotted when new error values are input
    artist_list = ax.get_children()
    for i in range(len(artist_list)):
        if isinstance(artist_list[i], matplotlib.collections.LineCollection):
            [artist_list[j].remove() for j in range(i-1, i+3)]
            ax.containers.remove(ax.containers[0])
            
    ax.errorbar(QQ, IQ, yerr = dIQ, marker = '.', linestyle = '', 
                color = 'k', label = 'SAS Data', capsize=2, zorder = 1)
    ax.relim()
    ax.autoscale()
    ax2 = fig.get_axes()[1]
    ax2.set_xlim(0.1/np.array(ax.get_xlim()))
    ax2.set_xticklabels(['{:g}'.format(2.5*val)
                         for val in ax2.get_xticks()])
    if len(ax.get_lines()) > 3:
        pf.show_legend(ax, reverse_order = True)
    else:
        pf.show_legend(ax)
    canvas.draw()


# This function used to plot the background-subtracted SAS data to the 'SAS Data' 
# and 'SAS Data vs. Fitted Result' plotting window
def plot_SANS_subtract(QQ_trim, IQ_trim, QQ_origin, bkgrd, fig, canvas):
    ax = fig.get_axes()[0]
    line_list = ax.get_lines()
    if len(line_list) > 3:      # remove previously plotted result
        for line in line_list:
            if line.get_label() == 'Subtracted Data' or line.get_label() == 'Background':
                line.remove()
    if bkgrd > 0:
        ax.plot(QQ_origin, np.ones(QQ_origin.shape)*bkgrd, 
                '-b', fillstyle = 'none', label = 'Background', zorder = 6,
                linewidth = 1.5)
    ax.plot(QQ_trim, IQ_trim, 'sr', fillstyle = 'none', label = 'Subtracted Data', zorder = 5)
    ax.relim()
    ax.autoscale()
    
    # redraw the top inverted horizontal axis
    ax2 = fig.get_axes()[1]
    ax2.set_xlim(0.1/np.array(ax.get_xlim()))
    ax2.set_xticklabels(['{:g}'.format(2.5*val)
                         for val in ax2.get_xticks()])
    pf.show_legend(ax, reverse_order= True)
    canvas.draw()


# This function plots the PDSP fitted SAS data to the 'SAS Data vs. Fitted Result'
# window together with the background-subtracted SAS data for comparison
def plot_SANS_fit(QQ_trim, IQ_plot, fig, canvas, which):
    ax = fig.get_axes()[0]
    line_list = ax.get_lines()
    _, legend_labels = ax.get_legend_handles_labels()
    
    # remove the previously plotted data points
    if which == 'input':
        if len(line_list) > 0:
            line_list[legend_labels.index('Subtracted Data')].remove()
        ax.plot(QQ_trim, IQ_plot, 'sr', fillstyle = 'none', 
                label = 'Subtracted Data', zorder = 1)
    elif which == 'result':
        if len(line_list) > 1:
            line_list[legend_labels.index('Fitted Result')].remove()
        ax.plot(QQ_trim, IQ_plot, '.b', fillstyle = 'full', 
                label = 'Fitted Result', zorder = 5)
    ax.relim()
    ax.autoscale()            
    
    # redraw the top inverted horizontal axis
    ax2 = fig.get_axes()[1]
    ax2.set_xlim(0.1/np.array(ax.get_xlim()))
    ax2.set_xticklabels(['{:g}'.format(2.5*val)
                         for val in ax2.get_xticks()])
   
    # make sure the legend of subtracted data always appears first
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) == 2 and labels[0] == 'Fitted Result':
        pf.show_legend(ax, reverse_order=True)
    else:
        pf.show_legend(ax)
    canvas.draw()


# This function plots the dV/dr vs r result to the 'dV/dr Plot' window
def plot_dVdr(rr, dV_dr, fig, canvas):
    ax = fig.get_axes()[0]

    # Removing previously plotted result
    artist_list = ax.get_children()
    for i in range(len(artist_list)):
        if isinstance(artist_list[i], matplotlib.collections.LineCollection):
            [artist_list[j].remove() for j in range(i-1, i+3)]
            ax.containers.remove(ax.containers[0])

    ax.errorbar(rr, dV_dr[:,0], dV_dr[:,1]*dV_dr[:,0], marker = 's', color = 'r',
                linestyle = '-', fillstyle = 'none', capsize=2)
    ax.relim()
    ax.autoscale()
    canvas.draw()


# This function plots the f(r) and SSA(R) result to the 'f(r) vs. r || SSA(R) vs. R'
# plot window
def plot_fr_SSA(rr, fr, SSA, num_pts_SSA_extrapolate, r_SSA_extrapolate, fig, canvas):
    ax = fig.get_axes()[0]
    
    # Removing previously plotted result
    artist_list = ax.get_children()
    for i in range(len(artist_list)):
        if isinstance(artist_list[i], matplotlib.collections.LineCollection):
            [artist_list[j].remove() for j in range(i-1, i+3)]
            ax.containers.remove(ax.containers[0])

    line_list = ax.get_lines()
    if len(line_list) > 0:
        [line.remove() for line in ax.get_lines()]
        [text.remove() for text in ax.texts]
    
    # Re-calculate the SSA value extrapolated to r_SSA_extrapolate from the 
    # linear fit of the SSA data at the first num_pts_SSA_extrapolate lowest
    # r values
    log_rr = np.log10(rr)
    log_SSA = np.log10(SSA[:,0])
    log_rr_4_extrapolate = log_rr[np.log10(rr)>=np.log10(r_SSA_extrapolate)][:num_pts_SSA_extrapolate]
    log_SSA_4_extrapolate = log_SSA[np.log10(rr)>=np.log10(r_SSA_extrapolate)][:num_pts_SSA_extrapolate]
    line_fit = np.polyfit(log_rr_4_extrapolate, log_SSA_4_extrapolate, 1)
    SSA_extrapolate = 10**np.polyval(line_fit, np.log10(r_SSA_extrapolate))

    # Plot f(r) and SSA(R) curve
    ax.errorbar(rr, fr[:,0], yerr = fr[:,1]*fr[:,0], marker = 's', linestyle = '-', 
              color = 'r', fillstyle = 'none', capsize=2, label = 'f(r)')
    ax.errorbar(rr, SSA[:,0], yerr = SSA[:,1]*SSA[:,0], marker = 'o', linestyle = '', 
              color = 'b', fillstyle = 'none', capsize=2, label = 'SSA(R)')
    
    # Turn the SSA data points used for extrapolation into solid markers
    ax.plot(10**log_rr_4_extrapolate, 10**log_SSA_4_extrapolate, 'ob')
    
    # Plot the linear model fitted through the selected SSA data
    ax.plot([r_SSA_extrapolate*0.8, 10**(np.max(log_rr_4_extrapolate)+0.2)], 
            10**np.polyval(line_fit, 
                           np.array([np.log10(r_SSA_extrapolate*0.8), 
                                     np.max(log_rr_4_extrapolate)+0.2])), 
            '-g',  linewidth = 1, zorder = 20)
    
    # Plot the pore radius at which the SSA is extrapolated to
    ax.plot([r_SSA_extrapolate, r_SSA_extrapolate], [np.min(fr), SSA_extrapolate*1e5], 
             '--k', linewidth = 1, label = 'r = {:.2f} nm'.format(r_SSA_extrapolate))
    
    # Write the extrapolated SSA value to the plot
    ax.text(r_SSA_extrapolate*1.3, SSA_extrapolate*1.6, 
            '{:s}'.format(sci_num_dot(SSA_extrapolate)) + '$\mathrm{cm^2/cm^3}$',
            fontsize = 11, weight = 'bold',
            ha = 'left', va = 'bottom')
    pf.show_legend(ax)
    ax.relim()
    ax.autoscale()
    canvas.draw()
    

# This function execute the PDSP model fitting routine, detailed explanation
# and mathematical background of the fitting routine is explained in the accompanied
# paper
def fit_PDSP_model(QQ, IQ, dIQ, pts_per_dec, lambda_, contrast, density_solid, 
                   r_SSA_extrapolate, num_pts_SSA_extrapolate, major_phase):
    # Check the input Q range, if it is shorter than 1 decade, return error
    if len(QQ) < 5 or np.log10(np.max(QQ)/np.min(QQ)) < 1:
        raise ValueError('Input Q range must span at least 1 decade and contain at least 5 data points')
        
    # Determining the log distance between each r value based on the number
    # of values per decade required for the result
    logR_del = 1/pts_per_dec
    
    # Determining the r range of the result based on the relationship r = 2.5/Q
    # r_min is extended further to r_min = 0.5/Q_max to improve the smoothness 
    # of the final result, which is trimmed back to R_min_original once the fit
    # is completed
    R_max_original = 10**(np.ceil(np.log10(2.5/np.min(QQ))/logR_del)*logR_del)
    R_min_original = 10**(np.floor(np.log10(2.5/np.max(QQ))/logR_del)*logR_del)
    logR_min = np.floor(np.log10(0.5/np.max(QQ))/logR_del)*logR_del
    logR_max = np.ceil(np.log10(2.5/np.min(QQ))/logR_del)*logR_del
    logR_1D = np.arange(logR_min, logR_max+logR_del/2, logR_del)

    # Determine the fraction value in Equation (2) for each pair of Q and r_i
    eq4_fraction_2D = calc_eq4_fraction(logR_1D, logR_del, QQ)
    
    # Determination of the starting value of IQ0i and dIQ0_i by assuming that 
    # the intensity contribution to a particular Q value consist solely of the 
    # intensity from r_i = 2.5/Q 
    R_1D = 10**logR_1D
    R_Q_pair_diff = np.abs(R_1D[:,np.newaxis] - 2.5/QQ[:, np.newaxis].T)
    R_Q_corr_pos = np.argmin(R_Q_pair_diff, axis = 1)
    IQ0_guessed = (IQ[R_Q_corr_pos]/
                   (eq4_fraction_2D[range(len(R_Q_corr_pos)), R_Q_corr_pos]))
    
    # Bring IQ0_guessed closer to result prior to fit
    IQ_guessed = np.sum(IQ0_guessed[:,np.newaxis]*eq4_fraction_2D,0)
    div_factor = np.median(IQ_guessed/IQ)
    IQ0_guessed /= div_factor

    # Setting arbitrary ranges for result, required for the least square fit 
    log_IQ0_guessed = np.log10(IQ0_guessed)
    log_IQ0_range = np.max(log_IQ0_guessed) - np.min(log_IQ0_guessed)
    log_IQ0_upper_bound = log_IQ0_guessed + log_IQ0_range
    log_IQ0_lower_bound = log_IQ0_guessed - log_IQ0_range
    
    # Start the fitting procedure
    QQ_IQ_slope = np.polyfit(np.log10(QQ), np.log10(IQ), 1)[0]
    IQ0_fitted = 10**(sci_opt.minimize(calc_Xi, log_IQ0_guessed, 
                                       bounds = [(min_val, max_val) 
                                                 for min_val, max_val
                                                 in zip(log_IQ0_lower_bound, 
                                                        log_IQ0_upper_bound)],
                                       args = (logR_1D, eq4_fraction_2D, 
                                               QQ, IQ, dIQ, QQ_IQ_slope, lambda_)).x)
    
    # Calculate I(Q) from the fitted data using Equation (2)
    IQ_fitted = np.sum(IQ0_fitted[:,np.newaxis]*eq4_fraction_2D,0)
    
    # Remove the result values at r_i < 2.5/Q_max
    non_extplted_pos = np.logical_and(R_1D - R_min_original >= -1e-10,
                                      R_1D - R_max_original <= 1e-10)
    R_1D = R_1D[non_extplted_pos]
    logR_1D = logR_1D[non_extplted_pos]
    IQ0_fitted = IQ0_fitted[non_extplted_pos]
    
    # Estimate fit error, then profject to fit result        
    fit_diff = np.abs(IQ_fitted - IQ)
    sigma = np.where(fit_diff > dIQ, fit_diff, dIQ)
    
    rr_corres_QQ = 2.5/QQ[::-1]
    dIQ_percent_data = sigma/IQ
    dIQ_percent_fit = sigma/IQ_fitted
    dIQ_percent = np.where(dIQ_percent_fit < dIQ_percent_data,
                           dIQ_percent_fit, dIQ_percent_data)
    R_min_1D = 10**(logR_1D - logR_del/2)
    R_max_1D = 10**(logR_1D + logR_del/2)
    R_2D = np.linspace(R_min_1D[1:-1], R_max_1D[1:-1], 51)
    dIQ0_percent = np.median(np.interp(R_2D, rr_corres_QQ, dIQ_percent[::-1]), axis = 0)
    dIQ0_percent = np.concatenate(([dIQ0_percent[0]], dIQ0_percent, [dIQ0_percent[-1]]))

    # Calculate sample properties enabled by the PDSP fit result and return result
    dR_1D = R_max_1D - R_min_1D
    f_dash_r = np.array([IQ0_fitted/np.sum(IQ0_fitted), dIQ0_percent]).T
    f_r = np.array([f_dash_r[:,0]/dR_1D, dIQ0_percent]).T
    IQ0_fitted = np.array([IQ0_fitted, dIQ0_percent]).T
    
    rr, SSA, dV_dr, phi, Vpore_avg, phi_on_Vtotal, SSA_extrapolate = \
        calc_PDSP_result(R_1D, f_r, f_dash_r, IQ0_fitted, contrast, density_solid,
                         r_SSA_extrapolate, num_pts_SSA_extrapolate, major_phase)
    return rr, IQ_fitted, IQ0_fitted, f_r, f_dash_r, SSA, dV_dr,\
        phi, Vpore_avg, phi_on_Vtotal, SSA_extrapolate


# This function calculate the optimise function Ξ in eqn 9
def calc_Xi(log_IQ0, logR_1D, integral_2D, QQ, IQ, dIQ, QQ_IQ_slope, lambda_):
    IQ0_slope = np.polyfit(logR_1D, log_IQ0, 1)[0]
    log_IQ0_normalised = log_IQ0 - logR_1D*IQ0_slope
    fancy_R = -np.sum(np.diff(log_IQ0_normalised)**2)
    IQ_calc = np.sum(10**log_IQ0[:,np.newaxis]*integral_2D,0)
    
    chi2 = (np.sum(((IQ*QQ**-QQ_IQ_slope)-(IQ_calc*QQ**-QQ_IQ_slope))**2/
                   (dIQ*QQ**-QQ_IQ_slope)**2))/len(QQ)
    
    chi2 = np.sum((np.log10(IQ*QQ**-QQ_IQ_slope)
                   -np.log10(IQ_calc*QQ**-QQ_IQ_slope))**2/
                  (dIQ/IQ)**2)/len(QQ)

    return chi2 - lambda_*fancy_R


# This function calculate the term following IQ0i in equation (2) for all pairs
# of r_i and Q
def calc_eq4_fraction(logR_1D, logR_del, QQ):
    num_subintervals = 600 # number of intervals for integral calculation set to 600
   
    # Creating pairs of Rmin_i and Rmax_i corresponding to each value of r_i
    logR_min_integral_1D = logR_1D - logR_del/2
    logR_max_integral_1D = logR_1D + logR_del/2
    R_min_integral_1D = 10**logR_min_integral_1D
    R_max_integral_1D = 10**logR_max_integral_1D
    
    # Divide the each pair of Rmin_i and Rmax_i into 600 equal space for integral
    # calculation
    R_integral_2D = np.linspace(R_min_integral_1D, R_max_integral_1D, 
                                  num_subintervals+1)
    
    # For each pair of Rmin_i and Rmax_i, calculate 
    # (i) dr in the integral, 
    # (ii) the mid point r value at each sub-interval dr, and
    # (iii) Qr, Vr and F(Qr) for every sub-interval between Rmin_i and Rmax_i
    dR_2D = np.diff(R_integral_2D, axis = 0)
    R_mid_2D = R_integral_2D[:-1,:] + 1/2*dR_2D
    Qr_3D = R_mid_2D[:,:,np.newaxis]*QQ
    Vr_2D = calc_Vsph(R_mid_2D)
    Fsph_3D = calc_Fsph(Qr_3D)
    
    # Calculate 
    # (i) the term inside the integral, 
    # (ii) the integral, and
    # (iii) the entire fraction following IQ0i in equation (2) for each
    # pair of r_i and Q
    subinterval_area_3D = (Vr_2D[:,:,np.newaxis]**2 * 
                            Fsph_3D * dR_2D[:,:,np.newaxis])
    RHS_integral_2D = np.sum(subinterval_area_3D, axis = 0)
    RHS_fraction_2D = (RHS_integral_2D/
                    (R_max_integral_1D[:,np.newaxis] - 
                      R_min_integral_1D[:,np.newaxis]))
    
    return RHS_fraction_2D


# This function calculate the spherical form factor, used in the integral
# calculation of equation (2) 
def calc_Fsph(Qr):
    return (3*(np.sin(Qr) - Qr*np.cos(Qr))/Qr**3)**2


# This function calculate the spherical pore volume, used in the integral
# calculation of equation (2) 
def calc_Vsph(radius):
    return 4/3*np.pi*radius**3
    

# This function calculate the structural properties from the PDSP fit result
def calc_PDSP_result(rr, f_r, f_dash_r, IQ_0, contrast, density_solid, 
                     r_SSA_extrapolate, num_pts_SSA_extrapolate, major_phase):
    
    dIQ0_percent = IQ_0[:,1]
    
    # Average pore volume - Equation 
    Vpore_avg = np.array([np.sum(4/3*np.pi*(rr*1e-8)**3*f_dash_r[:,0]), 
                          np.mean(dIQ0_percent)])
    
    # Porosity and pore concentration/density number - Equation 
    phiTimes1minusPhi_on_Vavg = np.mean(IQ_0[:,0]/contrast**2/f_dash_r[:,0]*1e48)
    roots = np.roots([1, -1, phiTimes1minusPhi_on_Vavg*Vpore_avg[0]])
    if major_phase == 'solid':
        phi = np.array([np.min(roots), np.mean(dIQ0_percent)*3])
    elif major_phase == 'void':
        phi = np.array([np.max(roots), np.mean(dIQ0_percent)*3])
    phi_on_Vavg = np.array([phi[0]/Vpore_avg[0], np.mean(dIQ0_percent)*4])

    # Specific surface area and differitial pore volume distribution - Equation
    SSA = np.array([np.cumsum((4*np.pi*(rr*1e-8)**2 * 
                               f_dash_r[:,0] * phi_on_Vavg[0])[::-1])[::-1],
                    dIQ0_percent]).T
    dV_dr = np.array([phi[0]/density_solid*f_r[:,0]*(4/3*np.pi*(rr*1e-8)**3)/Vpore_avg[0],
                      dIQ0_percent]).T
    
    # Convert pore radius from Angstrom to nm 
    rr /= 10
    
    # Extrapolate SSA to r_SSA_extrapolate
    log_rr = np.log10(rr)
    log_SSA = np.log10(SSA[:,0])
    log_rr_4_extrapolate = log_rr[rr>=r_SSA_extrapolate][:num_pts_SSA_extrapolate]
    log_SSA_4_extrapolate = log_SSA[rr>=r_SSA_extrapolate][:num_pts_SSA_extrapolate]
    line_fit, err_fit_square = np.polyfit(log_rr_4_extrapolate, log_SSA_4_extrapolate, 1, full=True)[:2]
    log_SSA_extrapolate = np.polyval(line_fit, np.log10(r_SSA_extrapolate))
    SSA_extrapolate = np.array([10**log_SSA_extrapolate, np.mean(dIQ0_percent)*4])
    # SSA_extrapolate = (10**log_SSA_extrapolate, 10**(err_fit_square/2))
    
    return rr, SSA, dV_dr, phi, Vpore_avg, phi_on_Vavg, SSA_extrapolate    
    

# Reformat the default scientific number returned by Python
def sci_num_dot(num, dec_pts = 2):
    base = int(np.log10(num))
    val = num/10**base
    return (r'${:.' + str(dec_pts) + r'f} \cdot 10^{{{:.0f}}}$').format(val, base)

  