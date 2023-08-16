import smuthi.simulation as sim
import smuthi.initial_field as init
import smuthi.layers as lay
import smuthi.particles as part
import smuthi.fields as flds
import smuthi.periodicboundaries as pb
import smuthi.periodicboundaries.post_processing as pbpost
import smuthi.utility.cuda
import numpy as np
from tqdm import tqdm
import sys
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import tempfile
import shutil
import imageio
import smuthi.postprocessing.far_field as ff


def only_tau(t_matrix, tau):
    t_matrix_shape = t_matrix.shape
    if tau == 0:
        t_matrix[t_matrix_shape[1]//2:] = 0.0
    if tau == 1:
        t_matrix[:t_matrix_shape[1]//2] = 0.0
    return t_matrix


def M103(t_matrix):
    t_matrix_shape = t_matrix.shape
    t_matrix[:t_matrix_shape[1]//2] = 0.0
    t_matrix[:10] = 0.0
    t_matrix[12:14] = 0.0
    return t_matrix
        
      
# in this file, all lengths are given in nanometers
#square lattice
a = 350 #lattice period
wl1 = a/0.46
wl2 = a/0.445 
WL = np.arange(wl1, wl2, 1)
polar_angles = np.arange(0, 90, 1) / 180 * np.pi
transmittance = np.zeros((len(WL),len(polar_angles)))
reflectance = np.zeros((len(WL),len(polar_angles)))
absorbance = np.zeros((len(WL),len(polar_angles)))
# extinction_coeffs = np.zeros(len(WL),len(polar_angles))






a1 = np.array([0, a, 0], dtype=float)
a2 = np.array([a, 0, 0], dtype=float)
pb.default_Ewald_sum_separation = pb.set_ewald_sum_separation(a1, a2)


for idx1, wl in enumerate(WL):
    neffmax = 3
    neffimag = 0.01
    waypoints = [0, 0.8, 0.8 - 1j * neffimag, 2.1 - 1j * neffimag, 2.1, neffmax]
    neff_discr = 1e-3  
    flds.default_Sommerfeld_k_parallel_array = flds.reasonable_Sommerfeld_kpar_contour(vacuum_wavelength=wl,
                                                                                   neff_waypoints=waypoints,
                                                                                   neff_resolution=neff_discr)
    flds.angular_arrays(angular_resolution=np.pi/360)
    # Scattering particle
    sphere = smuthi.particles.Sphere(position=[a/2, a/2, 0],
                                        refractive_index=np.sqrt(50),
                                        radius=100,
                                        l_max=3)

    t_matrix = sphere.compute_t_matrix(wl, 1)
    # print(t_matrix)
    mod_t_matrix = M103(t_matrix)
    # print(mod_t_matrix) 


    for idx2, beta in enumerate(tqdm(polar_angles, desc='Polar angle iterations  ', file=sys.stdout,
                                bar_format='{l_bar}{bar}| elapsed: {elapsed} ' 'remaining: {remaining}')):
        # suppress consol output
        # old_stdout = sys.stdout
        # sys.stdout = open(os.devnull, 'w')
        # old_stderr = sys.stderr
        # sys.stderr = open(os.devnull, 'w')

        initial_field = init.PlaneWave(vacuum_wavelength=wl,
                                    polar_angle=beta,
                                    azimuthal_angle=0,
                                    polarization=0,
                                    amplitude=1,
                                    reference_point=[0,0,0]) #TODO what is reference point?

        #square lattice
        a = 300 #lattice period
        a1 = np.array([0, a, 0], dtype=float)
        a2 = np.array([a, 0, 0], dtype=float)
        pb.default_Ewald_sum_separation = pb.set_ewald_sum_separation(a1, a2)


        # 1 homogeneous layer 
        layer_system = lay.LayerSystem(thicknesses=[0, 200, 0],
                                        refractive_indices=[1, 1, 1])


        # list of all scattering particles (only one in this case)
        sphere.t_matrix = mod_t_matrix
        one_sphere = [sphere]

        simulation = sim.Simulation(initial_field=initial_field,
                                    layer_system=layer_system,
                                    particle_list=one_sphere,
                                    periodicity=(a1, a2),
                                    ewald_sum_separation_parameter='default',
                                    number_of_threads_periodic=-2)

        simulation.run()
        print(sphere.t_matrix)
        # Since this moment we calculate extinction
        # Calculate total extinction of dipole (another multipoles are cut by only_l=1)
        # general_ecs = ff.extinction_cross_section(initial_field=initial_field,
                                            # particle_list=one_sphere,
                                            # layer_system=layer_system, only_l=1, only_tau=0, only_pol=0)

        # extinction_coeffs[idx] = general_ecs
        # print(general_ecs)
        # plane wave expansion of total transmitted field
        pwe_total_T = pbpost.transmitted_plane_wave_expansion(initial_field,
                                                                one_sphere,
                                                                layer_system,
                                                                a1, a2)
        # plane wave expansion of total reflected field
        pwe_total_R = pbpost.reflected_plane_wave_expansion(initial_field,
                                                            one_sphere,
                                                            layer_system,
                                                            a1, a2)


        # farfield objects       
        ff_T = pbpost.periodic_pwe_to_ff_conversion(pwe_total_T,
                                                    simulation.initial_field,
                                                    simulation.layer_system)
        ff_R = pbpost.periodic_pwe_to_ff_conversion(pwe_total_R,
                                                    simulation.initial_field,
                                                    simulation.layer_system)
            
        # power flow per area
        initial_power = pbpost.initial_plane_wave_power_per_area(initial_field, layer_system)
        transmitted_power = pbpost.scattered_periodic_ff_power_per_area(ff_T)
        reflected_power = pbpost.scattered_periodic_ff_power_per_area(ff_R)
        # T, R, A
        transmittance[idx1, idx2] = transmitted_power / initial_power
        reflectance[idx1, idx2] = reflected_power / initial_power
        absorbance[idx1, idx2] = 1 - transmittance[idx1, idx2] - reflectance[idx1, idx2]



np.savetxt('tmatrix.txt', sphere.t_matrix)
np.savetxt('T.txt', transmittance)
np.savetxt('R.txt', reflectance)
np.savetxt('A.txt', absorbance)
np.savetxt('beta.txt', polar_angles / np.pi * 180)
np.savetxt('WL.txt', WL)

# np.savetxt('magnetic_dipole_coeffs.txt', extinction_coeffs)

# np.savetxt('test.txt', sphere.compute_t_matrix(wl, 2))



