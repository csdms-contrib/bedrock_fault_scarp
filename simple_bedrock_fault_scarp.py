#! /usr/env/python

#-----------------------------------------------------------------------
# simple_bedrock_fault_scarp.py
#
# This is a two-dimensional numerical model that computes the topographic
# evolution of the facet slope in the footwall of an active normal fault.
# The model is described and illustrated in the following journal article:
#
# Tucker, G. E., S. W. McCoy, A. C. Whittaker, G. P. Roberts, 
# S. T. Lancaster, and R. Phillips (2011), Geomorphic significance of 
# postglacial bedrock scarps on normal-fault footwalls, J. Geophys. Res., 
# 116, F01022, doi:10.1029/2010JF001861.
#
# Note: to use a different input file, edit the last line in this file.
#
# Written by Greg Tucker, Spring 2010
#
# Licensing information:
#
# simple_bedrock_fault_scarp: Computes formation of a normal-fault
#   facet and fault scarp.
#
# Copyright (C) 2011 Gregory E. Tucker
#
# Developer can be contacted at:
#   Cooperative Institute for Research in Environmental Sciences (CIRES)
#   University of Colorado Boulder
#   Campus Box 399
#   Boulder, CO 80309 USA
#   Phone: (+1) 303-492-6985
#   Email: gtucker@colorado.edu
#
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the 
# Free Software Foundation; either version 2 of the License, or (at your 
# option) any later version.
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
# for more details.
# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software 
# Foundation, Inc., 51 Franklin Street, Fifth Floor, 
# Boston, MA 02110-1301 USA. 
#-----------------------------------------------------------------------

import numpy
from numpy import *

import model_parameter_dictionary
from model_parameter_dictionary import ModelParameterDictionary

import math
from math import tan, atan, radians, degrees, sin

import pylab
from pylab import plot, axis, draw, show, gca, xlabel, ylabel, legend


#-----------------------------------------------------------------------
class SimpleBedrockFaultScarp:
    """Simulates hillslope profile on a normal-fault-bounded
       mountain front."""
    #-------------------------------------------------------------------
    def __init__( self ):
    
        self.status = 'created'    # (OpenMI)
        
    # end of __init__
    
    #-------------------------------------------------------------------
    def initialize( self, input_file_name = None ):
        
        # read user-defined inputs
        if input_file_name != None:
            self.read_parameters_from_file( input_file_name )
        else:
            self.read_parameters_from_command_line()
            
        # set internal parameters
        self.current_time = 0.0     # current time in simulation, yr
        self.dx = float( self.L ) / float( self.nn )  # node spacing, m
        self.heaverate = self.throwrate / self.faultdip
        self.event_throw = self.throwrate * self.seismic_interval
        self.event_heave = self.heaverate * self.seismic_interval

        # set up arrays, etc.
        self.x = arange( self.dx / 2.0, self.L, self.dx ) # cell x locs, m
        self.z = zeros( self.nn )  # cell heights, m
        self.slope = zeros( self.nn )  # centered slope gradient
        
        # set up tectonics:
        # We need to know how deep the fault plane is at each node, and
        # encode this in "faultplane.z". To start, we set the depth to
        # the fault at the left edge ("faultdep0") at roughly 1/5th of 
        # the initial slope length (a bit more or less, depending on the
        # fault dip, which will typically be in the neighborhood of
        # 45 to 60 degrees). To handle the fact that the footwall nodes
        # move laterally as well as vertically, we track how far apart
        # are the pair of nodes that represent the fault at the surface
        # ("gap"). When this gap gets bigger than dx times a tolerance
        # factor (say, 2.0), a new node is added.
        self.faultdep0 = self.dx / 2.0
        self.faultplanez = -self.faultdep0 + self.faultdip * self.x
        self.gap = self.dx
        self.tolerance_factor = 2.0
        self.next_earthquake = 0.0
        
        # set up erosion-rate variation
        if self.opt_evar.startswith( 'Oxy' ):    # oxygen isotope curve
            self.setup_delO18()
        
        # set up plotting
        if self.opt_plot:
            self.next_plot = self.plot_interval
            self.max_len = self.L + self.heaverate * self.run_duration
            self.max_height = self.throwrate * self.run_duration
            self.plot_hillslope_profile()
            
        # end of initialize
            
    #-------------------------------------------------------------------
    def read_parameters_from_file( self, file_name ):
        
        # read the parameters into a ModelParameterDictionary
        mpd = ModelParameterDictionary()
        mpd.read_from_file( file_name )
        self.L = mpd.read_float( 'HILLSLOPE_LENGTH' )
        self.nn = mpd.read_int( 'NUMBER_OF_NODES' )
        self.erorate = mpd.read_float( 'EROSION_RATE' )
        self.Sc = mpd.read_float( 'THRESHOLD_SLOPE' )
        self.dt = mpd.read_float( 'DT' )
        self.throwrate = mpd.read_float( 'THROW_RATE' )
        self.faultdip = tan( radians( mpd.read_float( 'FAULT_DIP' ) ) )
        self.seismic_interval = mpd.read_float( 'SEISMIC_INTERVAL' )
        self.run_duration = mpd.read_float( 'RUN_DURATION' )
        self.opt_evar = mpd.read_string( 'OPT_ERO_VAR' )
        if self.opt_evar.startswith( 'Sin' ):
            self.amplitude = mpd.read_float( 'AMPLITUDE' )
            self.period = mpd.read_float( 'PERIOD' )
            self.phase = radians( mpd.read_float( 'PHASE' ) )
            self.erorate_mean = self.erorate
        elif self.opt_evar.startswith( 'Oxy' ):
            self.del18O_filename = mpd.read_string( 'DEL18O_FILENAME' )
            self.delO18power = mpd.read_float( 'DEL18O_POWER' )
            self.erorate_mean = mpd.read_float( 'EROSION_RATE' )
            self.min_erorate = mpd.read_float( 'MIN_ERORATE' )
        self.opt_plot = mpd.read_bool( 'OPT_PLOT' )
        if self.opt_plot:
            self.plot_interval = mpd.read_float( 'PLOT_INTERVAL' )
        self.opt_eps_plot = mpd.read_bool( 'OPT_EPS_PLOT' )
            
    #-------------------------------------------------------------------
    def read_parameters_from_command_line( self ):
        
        mpd = ModelParameterDictionary()
        self.L = mpd.read_float_cmdline( 'HILLSLOPE_LENGTH' )
        self.nn = mpd.read_int_cmdline( 'NUMBER_OF_NODES' )
        self.erorate = mpd.read_float( 'EROSION_RATE' )
        self.Sc = mpd.read_float( 'THRESHOLD_SLOPE' )
        self.dt = mpd.read_float( 'DT' )
        self.throwrate = mpd.read_float_cmdline( 'THROW_RATE' )
        self.faultdip = tan( radians( mpd.read_float_cmdline( 'FAULT_DIP' ) ) )
        self.seismic_interval = mpd.read_float( 'SEISMIC_INTERVAL' )
        self.run_duration = mpd.read_float_cmdline( 'RUN_DURATION' )
        self.opt_evar = mpd.read_string_cmdline( 'OPT_ERO_VAR' )
        if self.opt_evar.startswith( 'Sin' ):
            self.amplitude = mpd.read_float_cmdline( 'AMPLITUDE' )
            self.period = mpd.read_float_cmdline( 'PERIOD' )
            self.erorate_mean = self.erorate
        elif self.opt_evar.startswith( 'Oxy' ):
            self.del18O_filename = mpd.read_string_cmdline( 'DEL18O_FILENAME' )
            self.delO18power = mpd.read_float_cmdline( 'DEL18O_POWER' )
            self.erorate_mean = mpd.read_float_cmdline( 'EROSION_RATE' )
            self.min_erorate = mpd.read_float_cmdline( 'MIN_ERORATE' )
        self.opt_plot = mpd.read_bool_cmdline( 'OPT_PLOT' )
        if self.opt_plot:
            self.plot_interval = mpd.read_float_cmdline( 'PLOT_INTERVAL' )
        self.opt_eps_plot = mpd.read_bool( 'OPT_EPS_PLOT' )
    
    #-------------------------------------------------------------------
    def setup_delO18( self ):
        
        data = loadtxt( self.del18O_filename, dtype = 'float64' )
        self.num_delO18_vals = data.size / 2
        data.shape = (self.num_delO18_vals,2)
        data = flipud(data)
        self.ero_time = ( max(data[:,0]) - data[:,0] )
        
        # In the following, we transform the oxygen isotope curve into
        # an erosion rate curve through the following transformations:
        #  1. Shift it so that the 
        delO18 = data[:,1]
        delO18 = ( max( delO18 ) - delO18 )**self.delO18power
        delO18 = delO18 / mean(delO18)
        delO18 = delO18 * ( self.erorate_mean - self.min_erorate )
        self.delO18 = delO18 + self.min_erorate   # shift
        self.delO18[where(self.delO18<0.0)] = 0.0   # avoid zero erosion rate
        print 'Mean erorate:',mean(self.delO18)
        print 'Min erorate:',min(self.delO18)
        print 'Max erorate:',max(self.delO18)
        print 'Erorate since 7200bp:',mean(self.delO18[where(self.ero_time>=(max(self.ero_time)-7200.0))])

    #-------------------------------------------------------------------
    def run_model( self, input_file_name = None ):
    
        self.initialize( input_file_name )
        while self.current_time < self.run_duration:
            self.run_step()
        self.finalize()

    #-------------------------------------------------------------------
    def run_step( self ):
         
        # Report the current time
        print 'Time', self.current_time
        
        # If it's time, do some tectonics
        if self.current_time >= self.next_earthquake:
            self.launch_earthquake()
            self.next_earthquake += self.seismic_interval
            
        # Set the current erosion rate
        if self.opt_evar.startswith( 'Sin' ):
            self.erorate = self.erorate_mean \
                           + self.amplitude \
                           * sin( (self.phase + 2.0*pi)*self.current_time \
                                 / self.period )
        elif self.opt_evar.startswith( 'Oxy' ):
            self.erorate = interp( self.current_time, \
                                   self.ero_time, \
                                   self.delO18 )
            #print 'Time',self.current_time,'ero rate',self.erorate
        
        # Erode the hillslope
        self.erode_hillslope()
        
        # Update current time
        self.current_time += self.dt
        remaining_time = self.run_duration - self.current_time
        if remaining_time < self.dt:
            self.dt = remaining_time
            
        # Plot
        if self.opt_plot and self.current_time >= self.next_plot:
            self.plot_hillslope_profile()
            self.next_plot += self.plot_interval
        
    #-------------------------------------------------------------------
    def launch_earthquake( self ):
        
        # First, find the bit that should be offset. To do this, we
        # update the fault plane altitude. Next, we identify nodes
        # where the fault plane altitude is higher than the ground
        # surface, meaning that we're now on the footwall. We then
        # remember the footwall node indices as ranging from the
        # first footwall node to the right-most node in the domain
        # (we could use "fault_above_land" but there's a slight chance
        # that the fault will be exposed downslope and buried upslope)
        self.faultplanez = -self.faultdep0 + self.faultdip * self.x
        fault_above_land = where( self.faultplanez >= self.z )
        if size(fault_above_land)==0:
            print 'There are no fault points above land'
            sys.exit(1)
        footwall = arange( fault_above_land[0][0], self.nn )
        
        # now, offset those points: this is the "earthquake" that
        # comes every dt (we could set a separate recurrence 
        # interval for earthquakes)
        self.x[footwall] += self.event_heave
        self.z[footwall] += self.event_throw

        # find the distance between points at the fault
        first_footwall = footwall[0];
        last_hangingwall = footwall[0]-1  # assumes footwall[0]>1
        gap = self.x[first_footwall] - self.x[last_hangingwall]
        
        # if it's too big, add a point
        while gap > self.dx * self.tolerance_factor:
            newx = self.x[first_footwall] - self.dx
            newz = self.faultdip * newx - self.faultdep0
            #print 'adding a point at',newx,newz
            #print 'xc before:',self.x
            self.x = insert( self.x, first_footwall, newx )
            self.z = insert( self.z, first_footwall, newz )
            self.slope = insert( self.slope, first_footwall, 0.0 )
            self.nn += 1
            gap = self.x[first_footwall] - self.x[last_hangingwall]

    
    #-------------------------------------------------------------------
    def erode_hillslope( self ):
        
        # calculate the centered slopes
        self.slope[1:(self.nn-1)] = ( self.z[2:] - self.z[:(self.nn-2)] ) \
                                    / ( self.x[2:] - self.x[:(self.nn-2)] )
        self.slope[self.nn-1] = ( self.z[self.nn-1] - self.z[self.nn-2] ) \
                                / ( self.x[self.nn-1] - self.x[self.nn-2] )
              
        # Here we take a weighted average of the adjacent pair and 
        # the pair beyond that.
        weight = 0.67
        if self.nn>=5:
            self.slope[2:(self.nn-2)] = weight * self.slope[2:(self.nn-2)] + \
                                    (1.0-weight) * ( self.z[4:] - self.z[:(self.nn-4)] ) \
                                    / ( self.x[4:] - self.x[:(self.nn-4)] )                                
        
        # calculate the normalized offsets in x and z
        dz = - sqrt( 1.0 / ( self.slope**2.0 + 1.0 ) )
        dx = sqrt( 1.0 - dz**2.0 )
        
        # turn these into actual offsets using erosion rate and time step
        dz = dz * self.erorate * self.dt
        dx = dx * self.erorate * self.dt
        
        # set offset to zero where slope isn't steep enough
        dz[ where( self.slope < self.Sc ) ] = 0.0
        dx[ where( self.slope < self.Sc ) ] = 0.0
        
        # apply the offsets
        self.x = self.x + dx
        self.z = self.z + dz
        
        # delete points that get too close
        tol_sq = (0.5*self.dx)**2.0
        delx = diff( self.x )
        delz = diff( self.z )
        dist_sq = delx**2.0 + delz**2.0
        #print 'mds',min(dist_sq)
        while min( dist_sq ) < tol_sq:
                    
            # find the first offending pair
            too_close = where( dist_sq < tol_sq )
            #print 'tc',too_close
            too_close = too_close[0][0]   # take the first one (need the [0][0] because it's a tuple)
            #print 'pair too close at',too_close
            #print self.z.size, self.x.size
            
            # set "too_close" to the average position
            self.x[too_close] = 0.5 * ( self.x[too_close] + self.x[too_close+1] )
            self.z[too_close] = 0.5 * ( self.z[too_close] + self.z[too_close+1] )
            #print 'after: ',self.z.size, self.x.size
            
            # delete the second node
            self.x = delete( self.x, too_close+1 )
            self.z = delete( self.z, too_close+1 )
            self.slope = delete( self.slope, too_close+1 )
            self.nn = self.nn - 1
            
            # re-check distances
            delx = diff( self.x )
            delz = diff( self.z )
            dist_sq = delx**2.0 + delz**2.0
            
        # Occasionally, if the erosion rate varies strongly in time,
        # we might get an overhang. If so, remove the second node (which
        # is the overhanging node if it is higher than its neighbor, or
        # the underhanging node if lower).
        min_delx = min( delx )
        while min_delx < 0.0:
        
            # find the first overhang
            overhang = where( delx < 0.0 )
            first_overhang = overhang[0][0]
            #print 'Found an overhang at',first_overhang
            
            # delete the overhanging node
            self.x = delete( self.x, first_overhang+1 )
            self.z = delete( self.z, first_overhang+1 )
            self.slope = delete( self.slope, first_overhang+1 )
            self.nn = self.nn - 1
            
            # re-check distances
            delx = diff( self.x )
            delz = diff( self.z )
            min_delx = min( delx )
            #raw_input( 'waiting ...' )

        # Occasionally, it is also possible to get slope reversals.
        # Remove these as well.
        min_delz = min( delz )
        while min_delz < 0.0:
        
            # find the first overhang
            downhill = where( delz < 0.0 )
            first_downhill = downhill[0][0]
            #print 'Found a slope reversal at',first_downhill
            
            # delete the downhill node
            self.x = delete( self.x, first_downhill+1 )
            self.z = delete( self.z, first_downhill+1 )
            self.slope = delete( self.slope, first_downhill+1 )
            self.nn = self.nn - 1
            
            # re-check distances
            delx = diff( self.x )
            delz = diff( self.z )
            min_delz = min( delz )
            #raw_input( 'waiting ...' )
           

    #-------------------------------------------------------------------
    def plot_hillslope_profile( self ):
        
        gca().cla()
        plot( self.x, self.z, 'g.-' )
        for tl in gca().get_xticklabels():
            tl.set_fontsize(14)
        for tl in gca().get_yticklabels():
            tl.set_fontsize(14)
        axis( [ 0.0, self.max_len, 0.0, 1.2 * self.max_height ] )
        xlabel( 'Distance (m)', fontsize = 14 )
        ylabel( 'Height (m)', fontsize = 14 )
        #axis('equal')
        draw()
        

    #-------------------------------------------------------------------
    def finalize( self ):
    
        # calculate and report the mean slope angle
        dx = self.x[ self.nn-1 ] - self.x[0]
        dz = self.z[ self.nn-1 ] - self.z[0]
        slope = dz/dx
        if self.opt_evar.startswith( 'Oxy' ):
            self.erorate = sum( (diff(self.ero_time)/max(self.ero_time)) * self.delO18[0:(size(self.delO18)-1)] )
            print 'Mean erosion rate:',self.erorate
        print 'Total relief:',dz
        print 'Length:',dx
        print 'Average slope:',slope
        print 'Slope angle:',(180.0/pi)*arctan( slope )
        gamma = atan(self.faultdip)-(self.erorate*sin(atan(self.faultdip)))/self.throwrate
        print 'Predicted slope angle:',degrees( gamma )
        
        if self.opt_plot: 
            
            self.plot_fault_and_pred_slope( gamma )            
            show()
            
        if self.opt_eps_plot:
        
            self.setup_eps_output()
            plot( self.x, self.z, 'g.-', label = 'model' )
            xlabel( 'Distance (m)' )
            ylabel( 'Height (m)' )
            self.plot_fault_and_pred_slope( gamma )            
            pylab.savefig( 'simple_scarp_plot.eps' )
            
        # Write the final x and z coordinates to file
        savetxt( 'bedrock_scarp_x.dat', self.x )
        savetxt( 'bedrock_scarp_z.dat', self.z )
            
    #-------------------------------------------------------------------
    def plot_fault_and_pred_slope( self, gamma ):
    
            # plot where the fault plane would be
            x1 = (0.5*self.dx)/self.faultdip
            z1 = 0.0
            x2 = x1 + self.heaverate * self.run_duration
            z2 = z1 + self.throwrate * self.run_duration
            plot( [x1,x2], [z1,z2], 'k--', label = 'projected fault plane' )
            
            # plot the predicted slope
            L = ( self.throwrate * self.run_duration ) / sin( atan( self.faultdip ) )
            x2 = L * cos( gamma )
            z2 = L * sin( gamma )
            plot( [x1,x2], [z1,z2], 'k', label = 'predicted slope' )
            
            # Add a legend
            legend()
            
    #-------------------------------------------------------------------
    def setup_eps_output( self ):

        fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
        inches_per_pt = 1.0/72.27               # Convert pt to inch
        golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
        fig_width = fig_width_pt*inches_per_pt  # width in inches
        fig_height = fig_width*golden_mean      # height in inches
        fig_size =  [fig_width,fig_height]
        params = {'backend': 'ps',
                  'axes.labelsize': 10,
                  'text.fontsize': 10,
                  'legend.fontsize': 10,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8,
                  'text.usetex': True,
                  'figure.figsize': fig_size}
        pylab.rcParams.update(params)




bfs = SimpleBedrockFaultScarp()
bfs.run_model( 'fault_scarp_inputs_gisp2_3.txt' )
