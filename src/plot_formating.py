# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 13:07:48 2024

@author: NHUHA
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

matplotlib.rcParams['font.serif'] = 'Times New Roman'# default sans-serif font Times New Roman"
matplotlib.rcParams['font.family'] = 'serif' # ALWAYS use sans-serif fonts

# Set up font for axis labels and legends
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Times New Roman'
matplotlib.rcParams['mathtext.it'] = 'Times New Roman:italic'
matplotlib.rcParams['mathtext.bf'] = 'Times New Roman:bold'

matplotlib.rc('image', cmap='turbo') # Colour scheme for colour map plot

# Commonly used axis labels
scat_vec = 'SCATTERING VECTOR Q (\u212B' + r'$\mathbf{^{-1}}$)'
scat_int = r'SCATTERING INTENSITY $\mathbf{(cm^{-1})}$'

# Default line and marker properties
plt.rc('lines', markeredgewidth = 0.9)
plt.rc('lines', linewidth = 0.9)
plt.rc('lines', markersize = 4)

# Set axis labels to be semi-bold
plt.rc('axes', labelweight = 'semibold')

# Label font size
value_label_font_size = 8.5
axis_label_font_size = 11
tick_label_size_x = 9
tick_label_size_y = 9
axis_title_size = 12
title_size = 12
legend_size = 11
ytick_pad = 1.5
tick_length = 3.5


def show_legend(ax, legend_loc = 0, reverse_order = False):
    
    handles, labels = ax.get_legend_handles_labels()
    
    if reverse_order:
        handles = handles[::-1]
        labels = labels[::-1]
        
    ax.legend(handles = handles, labels = labels,
              fontsize = legend_size, facecolor = 'w', loc = legend_loc, 
              framealpha = 1, title_fontsize = legend_size,
              borderpad = 0.4, handletextpad = 1e-5,
              columnspacing = 0.01)
    ax.get_legend().get_frame().set_edgecolor('k')
    ax.get_legend().get_frame().set_linewidth(0.75)


def set_SAS_plot(fig, scale, legend_loc = 0):
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(scat_vec, labelpad = -1.5)
    ax.set_ylabel(scat_int, labelpad = -0.5)

    format_plot(fig, ax, inv_x = True,
                legend_pos= legend_loc, xtick_pad = 4, scale = scale)
    ax2 = fig.get_axes()[1]
    ax2.format_coord = lambda x,y: ax.format_coord(0.1/x,y).split('|')[0]
    

def set_dVdr_plot(fig, scale):
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('PORE RADIUS r (nm)', labelpad = 0)
    ax.set_ylabel(r'PORE VOLUME DISTRIBUTION $\mathbf{(cm^{3}/g}$' + '\u212B)', 
                  labelpad = 1, va = 'bottom')
    ax.yaxis.set_label_coords(-0.09,0.445)

    format_plot(fig, ax, xtick_pad=3, scale = scale)
    
    
def set_fr_SSA_plot(fig, scale):
    ax = fig.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('PORE RADIUS r & PROBE RADIUS R (nm)', labelpad = 0.5)
    ax.set_ylabel(r'f(r) & SSA(R) $\mathbf{(cm^{2}/cm^{3})}$', labelpad = -0.5)

    format_plot(fig, ax, xtick_pad=3, scale = scale)


def format_plot(fig, ax, inv_x = False, xtick_num = -1, 
                xtick_pad = 2, legend_pos = 1, adj_xtick = True, 
                adj_ytick = True, scale = 1):
    
    fig.set_dpi(round(150*scale))
    fig.set_edgecolor('k')
    fig.set_facecolor('w')
    fig.set_size_inches(w = 4, h = 3.5)
    
    fig.subplots_adjust(wspace = 0.2)
    fig.set_tight_layout(True)

    ax.spines[:].set_linewidth(0.75)
    
    ax.title.set_fontsize(title_size)
    
    ax.xaxis.label.set_fontsize(axis_title_size)
    ax.yaxis.label.set_fontsize(axis_title_size)
    
    ax.xaxis.offsetText.set_fontsize(axis_title_size)
    ax.yaxis.offsetText.set_fontsize(axis_title_size)
    
    ax.tick_params(axis = 'x', labelsize = tick_label_size_x, 
                   pad = xtick_pad, width = 0.75)
    ax.tick_params(axis = 'y', which = 'both', labelsize = tick_label_size_y,
                   pad = ytick_pad, width = 0.75)
    ax.xaxis.offsetText.set_fontsize(tick_label_size_x)
    ax.yaxis.offsetText.set_fontsize(tick_label_size_y)
    
    ax.xaxis.major.formatter._useMathText = True
    ax.yaxis.major.formatter._useMathText = True

    ax.minorticks_on()
    ax.tick_params(which = 'both', direction = 'in', 
                    length = tick_length, width = 0.75)
    ax.yaxis.set_tick_params(right = 'on', which = 'both', width = 0.75)

   
    ax.grid(which = 'major', visible = True, 
            color = [0, 0, 0], linewidth = 0.75, linestyle = 'dotted')
    
    if ax.get_legend_handles_labels()[0] != [] :
        ax.legend(fontsize = legend_size, facecolor = 'w', loc = legend_pos, 
                  framealpha = 1, title_fontsize = legend_size,
                  borderpad = 0.4, handletextpad = 1e-5, 
                  columnspacing = 0.01)
        ax.get_legend().get_frame().set_edgecolor('k')
        ax.get_legend().get_frame().set_linewidth(0.75)
        
    
    if inv_x == True:
        ax2 = ax.twiny()
        ax.set_zorder(1)
        ax2.set_zorder(2)
        ax.patch.set_alpha(0)
        
        ax2.set_xscale(ax.get_xscale())
        ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
        ax2.xaxis.label.set_fontsize(axis_title_size)
        ax2.tick_params(which = 'both', direction = 'in', length = tick_length, 
                        labelsize = tick_label_size_x, pad = xtick_pad,
                        width = 0.75, zorder = 1)
        
        ax_xlim = ax.get_xlim()
        
        ax2.set_xlim(0.1/np.array(ax_xlim))
        ax2.set_xticklabels(['{:g}'.format(2.5*val)
                             for val in ax2.get_xticks()])
        ax2.set_xlabel('PORE RADIUS = 2.5/Q (nm)', labelpad = 5)
        
    else:
        ax.xaxis.set_tick_params(top = 'on', which = 'both')
        
        
    def adj_tick(ax, axis):
        if axis == 'x':
            xyscale = ax.get_xscale()
            xylim = ax.get_xlim()
        elif axis == 'y':
            xyscale = ax.get_yscale()
            xylim = ax.get_ylim()
            
        if xyscale == 'log':
            xylim = np.log10(xylim)

        xy_max = np.max(xylim)
        xy_min = np.min(xylim)
        xy_span = xy_max - xy_min
        
        if xyscale == 'linear' or xy_span > 5:
            xy_span_convert = xy_span/10**(np.floor(np.log10(xy_span)))
            xy_min_convert = xy_min/10**(np.floor(np.log10(xy_span)))
            
            if xy_span_convert < 1.6:
                xytick_space_convert = 0.2
            elif xy_span_convert < 2.5:
                xytick_space_convert = 0.4
            elif xy_span_convert < 4:
                xytick_space_convert = 0.5
            elif xy_span_convert < 8:
                xytick_space_convert = 1
            else:
                xytick_space_convert = 2
                
            if xy_min_convert % xytick_space_convert != 0:
                xy_min_used_convert = ((xy_min_convert + xytick_space_convert) - 
                                       (xy_min_convert % xytick_space_convert))
            else:
                xy_min_used_convert = xy_min_convert
                
            xytick_space = xytick_space_convert*10**(np.floor(np.log10(xy_span)))
            xy_min_used = xy_min_used_convert*10**(np.floor(np.log10(xy_span)))
            xyticks_used = np.arange(xy_min_used, xy_max+xytick_space/2, xytick_space)
            
            if xyscale == 'log':
                if xy_min_used > 0 or xy_max < 0:
                    tick_list = 10**(xyticks_used)
                else:
                    tick_list = np.append(10**(xyticks_used), 1)
                    
                if axis == 'x':
                    ax.set_xticks(tick_list)
                elif axis == 'y':
                    ax.set_yticks(tick_list)
            else:
                if axis == 'x':
                    ax.set_xticks(xyticks_used)
                elif axis == 'y':
                    ax.set_yticks(xyticks_used)
            
    if adj_xtick:
        adj_tick(ax, axis = 'x')
    if adj_ytick:
        adj_tick(ax, axis = 'y')
