#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################################################################
#              Program to visualize the partciles output by 
#                                SIMP.f90.
#######################################################################

### Module importations ###
import numpy as np
import sys
import os
import plotting_functions as pf
import matplotlib.pyplot as plt

#Ask format of output
g_fmt = input("""Graphic display format (choose a number):
                    1 : Particle display
                    2 : Interpolation display

Enter desired value:""")

if not ( g_fmt == '1' or g_fmt == '2' ) :
    sys.exit("Wrong format value.")

### PATH to file: ###
#Ask name of run folder output
input_folder = input("""Enter the name of the run:""")
cwd = os.getcwd()
data_path = cwd +  '/fortran/data/'+ input_folder +'/run' 
output_folder = 'plots/' + input_folder

#Ask if the run has a boundary file:
boundary_file = input("""Does the run has boundary particles (T/F)?:""")


### Read momentum ###
momentum_data = open(data_path + '/particle_momentum.csv', mode='r')
particles_momentum = np.loadtxt( momentum_data, delimiter=',' )
p_tag, x, y, vx, vy = np.hsplit(particles_momentum,5)

#Get the number of time step:
ntimestep = sum(particles_momentum[:,0]==1) 

#Close and delete no longer useful data
momentum_data.close()
del(particles_momentum)

### Read state ###
state_data = open(data_path + '/particle_state.csv', mode='r')
particles_state = np.loadtxt( state_data, delimiter=',' )
p_tag, mass, p, hsml, A, h, epsilon_11, epsilon_12, epsilon_22, \
       sigma_11, sigma_12, sigma_22 = np.hsplit(particles_state, 12) 
p_tag.astype(int, copy=False)

#Close and delete no longer useful data
state_data.close()
del(particles_state)

### Read other ###
other_data = open(data_path + '/particle_other.csv', mode='r')
particles_other = np.loadtxt( other_data, delimiter=',')
npart, nbpart, time_step = np.hsplit(particles_other, 3)
npart  = int(npart[1])
nbpart = int(nbpart[1])
time_step = time_step[:,0]

#Close and delete no longer useful data
other_data.close()
del(particles_other)

### Read boundary ###
if (boundary_file == 'T'):
    boundary_data = open(data_path + '/boundary_particles.csv', mode='r')
    boundary_particles = np.loadtxt( boundary_data, delimiter=',' )
    bp_tag , xb, yb, massb = np.hsplit(boundary_particles, 4)
    bp_tag.astype(int, copy=False)

    #Close and delete no longer useful data
    boundary_data.close()
    del(boundary_particles)
    
    #Get space limit
    space_limits = [np.nanmin([np.nanmin(xb) - 0.1 * (np.nanmax(xb) - np.nanmin(xb)), 
                       np.nanmin(x) - 0.1 * (np.nanmax(x) - np.nanmin(x))]),
               np.nanmax([np.nanmax(xb) + 0.1 * (np.nanmax(xb) - np.nanmin(xb)),
                      np.nanmax(x) + 0.1 * (np.nanmax(x) - np.nanmin(x))]),
               np.nanmin([np.nanmin(yb) - 0.1 * (np.nanmax(yb) - np.nanmin(yb)), 
                       np.nanmin(y) - 0.1 * (np.nanmax(y) - np.nanmin(y))]),
               np.nanmax([np.nanmax(yb) + 0.1 * (np.nanmax(yb) - np.nanmin(yb)),
                      np.nanmax(y) + 0.1 * (np.nanmax(y) - np.nanmin(y))])]

elif(boundary_file == 'F'):

    #Create variables with empty values
    bp_tag , xb, yb, massb = np.nan, np.nan, np.nan, np.nan

    #Get space limit
    space_limits = [ np.nanmin(x) - 0.1 * (np.nanmax(x) - np.nanmin(x)),
                     np.nanmax(x) + 0.1 * (np.nanmax(x) - np.nanmin(x)), 
                     np.nanmin(y) - 0.1 * (np.nanmax(y) - np.nanmin(y)),
                     np.nanmax(y) + 0.1 * (np.nanmax(y) - np.nanmin(y))]

hsml_multiplicity = 3
#Get marker size
dummy, ax = plt.subplots(figsize=(20,20))
ax.set_xlim([space_limits[0], space_limits[1]])
ax.set_ylim([space_limits[2], space_limits[3]])
M = ax.transData.get_matrix()
xscale = M[0,0]
yscale = M[1,1]
m_size = hsml * xscale * hsml * yscale / (2*hsml_multiplicity*1000) ** 2
m_size[np.where(m_size> 3*m_size[0])] = 3*m_size[0]
plt.close()

#Get principal stress componant
sigma_I  = 0.5 * ( sigma_11 + sigma_22)
sigma_II = 0.5 * ( np.sqrt( (sigma_11-sigma_22) ** 2
           + 4 * sigma_12 ** 2) )

#Get principal strain rate componant
epsilon_I  = 0.5 * (epsilon_11 + epsilon_22 )
epsilon_II = 0.5 * (np.sqrt( (epsilon_22 - epsilon_11) ** 2 
             + 4 * epsilon_12 ** 2))

#Get velocity amplitude
v = np.sqrt(vx**2 +vy**2)

#Loop until stop
cdt = 1
while(cdt):
    # Ask the time step where the graphic are desired
    time_selected = input("""Enter the time steps of the desired graphics in seconds:
                    ALL : Create all possible figures
                    t1,t2,...,tn : Create figures for the closest times
                    EXIT : Quit the loop
                    
                    Enter desired value:""")

    if (time_selected == 'ALL'): 
        

        ### Ploting iteration dependant figure###
        for i in range(0, ntimestep) :
    
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/velocities/'
            if not os.path.exists(path) : os.makedirs(path)
            
            #Create velocity figures
            pf.speed_amplitude_evolution(
                                x[i * npart : (i+1) * npart], 
                                y[i * npart : (i+1) * npart],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                vx[i * npart : (i+1) * npart],
                                vy[i * npart : (i+1) * npart],
                                [np.min(v),np.max(v)],
                                time_step[i],
                                path,
                                space_limits)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/pressure/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create velocity figures
            pf.pressure_evolution(
                                x[i * npart : (i+1) * npart], 
                                y[i * npart : (i+1) * npart],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                p[i * npart : (i+1) * npart],
                                [np.min(p),np.max(p)],
                                time_step[i],
                                path,
                                space_limits)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/stresses/'
            if not os.path.exists(path) : os.makedirs(path)
    
            # Ice stress field figure creation function
            pf.stress_space(h[i * npart : (i+1) * npart],
                            A[i * npart : (i+1) * npart],
                            sigma_I[i * npart : (i+1) * npart],
                            sigma_II[i * npart : (i+1) * npart],
                            time_step[i],
                            path)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/stressfield/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create strain rate figures
            pf.stress_evolution(
                                x[i * npart : (i+1) * npart,0], 
                                y[i * npart : (i+1) * npart,0],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                (sigma_I[i * npart : (i+1) * npart,0],
                                sigma_II[i * npart : (i+1) * npart,0]),
                                np.array([np.min(sigma_I),np.max(sigma_I),
                                np.min(sigma_II),np.max(sigma_II)]),
                                time_step[i],
                                path,
                                space_limits,
                                g_fmt)
            
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/thickness/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create thickness figures
            pf.thickness_evolution(
                            x[i * npart : (i+1) * npart,0], 
                            y[i * npart : (i+1) * npart,0],
                            m_size[i * npart : (i+1) * npart],
                            xb, 
                            yb,
                            h[i * npart : (i+1) * npart,0],
                            np.array([np.min(h),np.max(h)]),
                            time_step[i],
                            path,
                            space_limits,
                            g_fmt)
    
            
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/concentration/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create concentration figures
            pf.concentration_evolution(
                            x[i * npart : (i+1) * npart,0], 
                            y[i * npart : (i+1) * npart,0],
                            m_size[i * npart : (i+1) * npart],
                            xb, 
                            yb,
                            A[i * npart : (i+1) * npart,0],
                            np.array([np.min(A),np.max(A)]),
                            time_step[i],
                            path,
                            space_limits,
                            g_fmt)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/strainrates/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create strain rate figures
            pf.strain_evolution(
                                x[i * npart : (i+1) * npart,0], 
                                y[i * npart : (i+1) * npart,0],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                (epsilon_I[i * npart : (i+1) * npart,0],
                                epsilon_II[i * npart : (i+1) * npart,0]),
                                np.array([np.min(epsilon_I),np.max(epsilon_I),
                                np.min(epsilon_II),np.max(epsilon_II)]),
                                time_step[i],
                                path,
                                space_limits,
                                g_fmt)
        

        break

    elif(time_selected == 'EXIT'):
        break 

    else:
        
        #Ask unit 
        unit = input("""What unit are the time selected (d,h,m,s):""")

        if (unit == 'd'):
            unit = 24 

        elif (unit == 'h'):
            unit = 1
        
        elif (unit == 'm'):
            unit = 1 / 60

        elif (unit == 's'):
            unit = 1/ 3600

        #Convert the time selected in time step
        time_selected = time_selected.split(',')

        time_selected = [float(i) * unit for i in time_selected]
        iteration_number = [np.argmin(abs(time_step-i))  for i in time_selected]

        ### Ploting iteration dependant figure###
        for i in iteration_number :

            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/velocities/'
            if not os.path.exists(path) : os.makedirs(path)

            #Create velocity figures
            pf.speed_amplitude_evolution(
                                x[i * npart : (i+1) * npart], 
                                y[i * npart : (i+1) * npart],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                vx[i * npart : (i+1) * npart],
                                vy[i * npart : (i+1) * npart],
                                [np.min(v),np.max(v)],
                                time_step[i],
                                path,
                                space_limits)
    
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/pressure/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create velocity figures
            pf.pressure_evolution(
                                x[i * npart : (i+1) * npart], 
                                y[i * npart : (i+1) * npart],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                p[i * npart : (i+1) * npart],
                                [np.min(p),np.max(p)],
                                time_step[i],
                                path,
                                space_limits)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/stresses/'
            if not os.path.exists(path) : os.makedirs(path)
    
            # Ice stress field figure creation function
            pf.stress_space(h[i * npart : (i+1) * npart],
                            A[i * npart : (i+1) * npart],
                            sigma_I[i * npart : (i+1) * npart],
                            sigma_II[i * npart : (i+1) * npart],
                            time_step[i],
                            path)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/stressfield/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create strain rate figures
            pf.stress_evolution(
                                x[i * npart : (i+1) * npart,0], 
                                y[i * npart : (i+1) * npart,0],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                (sigma_11[i * npart : (i+1) * npart,0],
                                sigma_12[i * npart : (i+1) * npart,0]),
                                np.array([np.min(sigma_11),np.max(sigma_11),
                                np.min(sigma_12),np.max(sigma_12)]),
                                time_step[i],
                                path,
                                space_limits,
                                g_fmt)
            
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/thickness/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create thickness figures
            pf.thickness_evolution(
                            x[i * npart : (i+1) * npart,0], 
                            y[i * npart : (i+1) * npart,0],
                            m_size[i * npart : (i+1) * npart],
                            xb, 
                            yb,
                            h[i * npart : (i+1) * npart,0],
                            np.array([np.min(h),np.max(h)]),
                            time_step[i],
                            path,
                            space_limits,
                            g_fmt)
    
            
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/concentration/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create concentration figures
            pf.concentration_evolution(
                            x[i * npart : (i+1) * npart,0], 
                            y[i * npart : (i+1) * npart,0],
                            m_size[i * npart : (i+1) * npart],
                            xb, 
                            yb,
                            A[i * npart : (i+1) * npart,0],
                            np.array([np.min(A),np.max(A)]),
                            time_step[i],
                            path,
                            space_limits,
                            g_fmt)
    
            #Create directory if it doesn't exist
            path = cwd + '/python/' + output_folder + '/strainrates/'
            if not os.path.exists(path) : os.makedirs(path)
    
            #Create strain rate figures
            pf.strain_evolution(
                                x[i * npart : (i+1) * npart,0], 
                                y[i * npart : (i+1) * npart,0],
                                m_size[i * npart : (i+1) * npart],
                                xb, 
                                yb,
                                (epsilon_I[i * npart : (i+1) * npart,0],
                                epsilon_II[i * npart : (i+1) * npart,0]),
                                np.array([np.min(epsilon_I),np.max(epsilon_I),
                                np.min(epsilon_II),np.max(epsilon_II)]),
                                time_step[i],
                                path,
                                space_limits,
                                g_fmt)
