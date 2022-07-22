#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################################################################
#              Module of plotting functions to use in the  
#                     data_visualization program
#######################################################################

### Module importations ###
import matplotlib.pyplot as plt
import matplotlib
#from labellines import labelLine, labelLines
import numpy as np
from scipy.interpolate import griddata
from cmcrameri import cm

### Function definition ###

# Ice field speed figure creation function
def speed_amplitude_evolution(x, y, size, xb, yb, vx, vy, v_max, time, output, limits):
    
    #Get speed magnitude
    v = np.sqrt(vx ** 2 + vy ** 2)

    #Chosing particle speed to show the direction
    N = len(x) / 500
    mask = np.where( np.random.randint(0,N, size = len(x)) == 0)

    #Field figure
    plt.figure(figsize=(20,20))
    plt.scatter(xb,yb, color='k', s=(8))
    plt.scatter(x, y, c = v, cmap = plt.get_cmap('viridis'), s=size, 
                vmin = v_max[0] , vmax= v_max[1])
    plt.quiver(x[mask],y[mask],vx[mask],vy[mask])
    clb = plt.colorbar( format='%.2e')
    clb.ax.set_title('[cm/s]', fontsize=(24))
    clb.ax.tick_params(labelsize=20)
    plt.clim(v_max[0], v_max[1])
    plt.xlabel('X coordinates [km]',fontsize=(24))
    plt.ylabel('Y coordinates [km]',fontsize=(24))
    plt.xlim([limits[0], limits[1]])
    plt.ylim([limits[2], limits[3]])
    plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
    plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
    plt.axis('scaled')
    plt.grid()

    plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
    
    plt.savefig(output + 'iceSpeedField_{:.8f}.jpg'.format(time), 
                bbox_inches='tight')
    plt.close()











def pressure_evolution(x, y, size, xb, yb, p, p_max, time, output, limits):

    #Field figure
    plt.figure(figsize=(20,20))
    plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
    plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
    plt.scatter(xb,yb, color='k', s=(8))
    plt.scatter(x, y, c = p, cmap = plt.get_cmap('viridis'), s=size, 
                vmin = p_max[0] , vmax= p_max[1])
    clb = plt.colorbar( format='%.2e')
    clb.ax.set_title('[Pa]', fontsize=(24))
    clb.ax.tick_params(labelsize=20)
    plt.xlabel('X coordinates [km]',fontsize=(24))
    plt.ylabel('Y coordinates [km]',fontsize=(24))

    plt.grid()
    plt.axis('scaled')

    plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
    
    plt.savefig(output + 'icePressureField_{:.8f}.jpg'.format(time), 
                bbox_inches='tight')
    plt.close()











# Ice stress field figure creation function
def stress_space(h,A, sigma_I, sigma_II, time, output):
    
    #Normalize stress invariant by pressure
    p = 27500 * h * np.exp(-20*(1-A))

    #Create ellipse coordinate for Viscous-Plastic
    a = 0.5 
    b = 0.25
    x = np.linspace(-1,0,100)
    y = np.sqrt(( 1 - (x + 0.5) ** 2 / a ** 2 ) * b ** 2) 
    
    #Stress figure
    plt.figure(figsize=(16,10))
    plt.scatter(sigma_I/p,sigma_II/p, color='r', s=(8))
    plt.plot(x,y, '--k')
    plt.plot(x,-y, '--k')

    #Create lines coordinate for Mohr-Coulomb
    x = np.linspace(-1,0.1,100)
    y = np.tan(np.pi/4)* (x) # - 0.1)

    plt.plot(x,y, '--k')
    plt.plot(x,-y, '--k')

    plt.xlim([-1,0])
    plt.ylim([-0.5,0.5])
    plt.xlabel(r'Principal stress $\sigma_I/P$ [Pa]',fontsize=(24))
    plt.ylabel(r'Principal stress $\sigma_{II}/P$ [Pa]',fontsize=(24))
    plt.grid()
    plt.axis('scaled')

    plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
    
    plt.savefig(output + 'iceStressField_{:.8f}.jpg'.format(time), 
                bbox_inches='tight')
    plt.close()












# Ice field mass distribution figure creation function
def mass_evolution(x, y, xb, yb, mass, m_max, time, output, limits, g_fmt):
    
    if (g_fmt == '1'):
        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter(x, y, c = mass, cmap = plt.get_cmap('Blues'), s=(24),
                    vmin = m_max[0], vmax = m_max[1])
        clb = plt.colorbar( format='%.2e' )
        clb.ax.set_title('[kg]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceMassField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()


    if (g_fmt == '2' ):
        #Create grid
        grid_x = np.linspace(limits[0], limits[1], 200)
        grid_y = -  np.linspace(limits[2], limits[3], 200)
        grid_x, grid_y = np.meshgrid(grid_x,grid_y)
        grid_z = griddata((x,y), mass, (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('Blues'), vmin = m_max[0],
            vmax = m_max[0], extent = limits)
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('Blues'), 
                    levels = np.linspace(m_max[0] , m_max[1],10))
        clb = plt.colorbar( format='%.2e' )
        clb.ax.set_title('[kg]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceMassField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()









############################ Thickness figures #############################

# Ice field thickness distribution figure creation function
def thickness_evolution(x, y, size, xb, yb, h, h_max, time, output, limits, g_fmt):
    
    if (g_fmt == '1'):
        #Thickness Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter(x, y, c = h, cmap = plt.get_cmap('Blues'), s=size,
                    vmin = h_max[0] , vmax = h_max[1])
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title('[m]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceThicknessField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Thickness anomaly Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter(x,y, c = h-1, cmap = plt.get_cmap('RdBu'), s=size,
                vmin = -max(abs(h_max-1))/1 , vmax = max(abs(h_max-1))/1)
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'$\Delta h$ [m]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceThicknessAnomalyField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()


    if (g_fmt == '2' ):
        #Create grid
        grid_x = np.linspace(limits[0], limits[1], 200)
        grid_y = - np.linspace(limits[2], limits[3], 200)
        grid_x, grid_y = np.meshgrid(grid_x, grid_y)

        # Thickness field
        grid_z = griddata((x,y), h, (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('Blues'), vmin = h_max[0],
            vmax = h_max[1], extent = limits)
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('Blues'), 
                    levels = np.linspace(h_max[0] , h_max[1],10))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title('[m]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceThicknessField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Thickness anomaly field
        grid_z = griddata((x,y), h - 1, (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('RdBu'), extent = limits,
            vmin = -max(abs(h_max-1)), vmax = max(abs(h_max-1)))
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('RdBu'), 
                levels = np.linspace(-max(abs(h_max-1)) , max(abs(h_max-1))))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'$\Delta h$ [m]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceThicknessAnomalyField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()







########################## Concentration figures ###########################

# Ice field concentration distribution figure creation function
def concentration_evolution(x, y, size, xb, yb, A, A_max, time, output, limits, g_fmt):
    
    if (g_fmt == '1'):
        #Concentration Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter(x, y, c = A, cmap = plt.get_cmap('Blues'), s=size,
                    vmin = A_max[0] , vmax = A_max[1])
        clb = plt.colorbar( format='%.2e')
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceConcentrationField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Concentration anomaly field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter( x, y, c = A-1, cmap = plt.get_cmap('RdBu'), s=size,
                vmin = -max(abs(A_max-1)) , vmax = max(abs(A_max-1)))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'$\Delta A$', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceConcentrationAnomalyField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

    if (g_fmt == '2' ):

        #Create grid
        grid_x = np.linspace(limits[0], limits[1], 200)
        grid_y = -  np.linspace(limits[2], limits[3], 200)
        grid_x, grid_y = np.meshgrid(grid_x,grid_y)

        #Concentration field
        grid_z = griddata((x,y), A, (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('Blues'), vmin = A_max[0],
            vmax = A_max[1], extent = limits)
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('Blues'), 
                    levels = np.linspace(A_max[0] , A_max[1],10))
        clb = plt.colorbar( format='%.2e')
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceConcentrationField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Concentration anomaly field
        grid_z = griddata((x,y), A - 1, (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('RdBu'),extent = limits,
            vmin = -max(abs(A_max-1)), vmax = max(abs(A_max-1)))
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('RdBu'), 
                levels = np.linspace(-max(abs(A_max-1)) , max(abs(A_max-1))))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'$\Delta A$', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'iceConcentrationAnomalyField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()








# Ice field strain distribution figure creation function
def strain_evolution(x, y, size, xb, yb, epsilon, e_max, time, output, limits, g_fmt):
    

    if (g_fmt == '1'):

        #Field figure for epsilon_II
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter( x, y, c = epsilon[1], cmap = cm.tokyo, s=size,
                    norm=matplotlib.colors.LogNorm( vmin = 10**-8 , vmax = 10**-4))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'[$s^{-1}$]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'e_IIiceStrainField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Field figure epsilon_I
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter(x, y, c = epsilon[0], cmap = cm.vik, s =size,
             norm=matplotlib.colors.SymLogNorm(vmin = -10**-5 , vmax = 10**-5, linthresh=10**-7, base=10))
        clb = plt.colorbar( format='%.2e', ticks = [-10**-5,-10**-6,-10**-7, 0,10**-7,10**-6,10**-5])
        clb.ax.set_title(r'[$s^{-1}$]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'e_IiceStrainField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

    if (g_fmt == '2' ):
        #Create grid
        grid_x = np.linspace(limits[0], limits[1], 200)
        grid_y = -  np.linspace(limits[2], limits[3], 200)
        grid_x, grid_y = np.meshgrid(grid_x,grid_y)

        #Figure for epsilon_II
        grid_z = griddata((x,y), epsilon[1], (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('GnBu'), vmin = e_max[2],
            vmax = e_max[3], extent = limits)
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('GnBu'), 
                    levels = np.linspace(e_max[2] , e_max[3],10))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'[$s^{-1}$]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'e_IIiceStrainField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Figure epsilon_I
        grid_z = griddata((x,y), epsilon[0], (grid_x, grid_y), method='cubic')
        grid_z = np.nan_to_num(grid_z,nan=0.)

        #Field figure
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        if (np.all(grid_z) == 0) :
            plt.imshow(grid_z, cmap = plt.get_cmap('RdBu'), 
            vmin = - max(abs(e_max[:1])), vmax = max(abs(e_max[:1])),
            extent = limits)
        else :
            plt.contourf(grid_x, grid_y, grid_z, cmap = plt.get_cmap('RdBu'), 
            levels = np.linspace( - max(abs(e_max[:1])) , 
            max(abs(e_max[:1])), 10))
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'[$s^{-1}$]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 'e_IiceStrainField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()







# Ice field strain distribution figure creation function
def stress_evolution(x, y, size, xb, yb, sigma, s_max, time, output, limits, g_fmt):
    
    if (g_fmt == '1'):

        #Field figure for sigma_II
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter( x, y, c = sigma[1], cmap = cm.lapaz, s=size,
                    norm=matplotlib.colors.LogNorm( vmin = 10**1 , vmax = 10**5))
                    #vmin = s_max[2]/1000 , vmax = s_max[3]/1000)
        clb = plt.colorbar( format='%.2e')
        clb.ax.set_title(r'[$N \cdot m^{-1}$]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 's_IIiceStressField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()

        #Field figure sigma_I
        plt.figure(figsize=(20,20))
        plt.scatter(xb,yb, color='k', s=(8))
        plt.scatter(x, y, c = sigma[0], cmap = cm.cork, s =size,
            norm=matplotlib.colors.SymLogNorm(vmin = -10**5 , vmax = 10**5, linthresh=10**3, base=10))
            #vmin = - max(abs(s_max[:1]))/1000, vmax = max(abs(s_max[:1]))/1000)
        clb = plt.colorbar( format='%.2e',ticks = [-10**5,-10**4,-10**3, 0,10**3,10**4,10**5])
        clb.ax.set_title(r'[$N \cdot m^{-1}$]', fontsize=(24))
        clb.ax.tick_params(labelsize=20)
        plt.xlabel('X coordinates [km]',fontsize=(24))
        plt.ylabel('Y coordinates [km]',fontsize=(24))
        plt.xlim([limits[0], limits[1]])
        plt.ylim([limits[2], limits[3]])
        plt.xticks(np.linspace(limits[0], limits[1], 5), fontsize=(24))
        plt.yticks(np.linspace(limits[2], limits[3], 10), fontsize=(24))
        plt.grid()
        plt.axis('scaled')

        plt.title('Simulation time [h] : {:.8f}'.format(time),fontsize=(32))
        
        plt.savefig(output + 's_IiceStressField_{:.8f}.jpg'.format(time), 
                    bbox_inches='tight')
        plt.close()


