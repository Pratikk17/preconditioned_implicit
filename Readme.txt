 
These matlab codes are intended to reproduce the results from the paper

"Preconditioned Implicit Time Integration Schemes For Maxwell's Equations On Locally Refined Grids" by
 Marlis Hochbruck, Jonas Koehler, and Pratik M. Kumbhar.

 These codes are tested with Matlab R2022a.
 
 To regenerate figures in this paper, open the matlab and run the following scripts:
 
 1) Run "fig_5_2_fov.m" to create the plot on fields of values. The eigenvalues and fields of values are already computed (to save time) and stored in the folder "matrices/gmsh_square_meshes". If these stored files are deleted, this script can then also compute eigenvalues and field of values.
 
 2) Run "fig_5_3_optimization.m" to generate the plot on dependence of 1/phi_0 on gamma for fixed values of alpha and tau. This script uses already stored data in "optimization_data_fig_5_3.m". However, this data can also be generated using "data_optimization_fig_5_3.m".
 
 3) Run "fig_5_4_PQMR_QMR_dt_5em2.m" to plot the relative error of QMR and pQMR against number of iterations for different meshes.
 
 4) Run "fig_5_5_errorbounds_dt_1em2.m" creates figure 5.5, which numerically verifies the error bound in Theorem 4.5.
 
 We thank Andreas Sturm for the generation of meshes in Figure 5.1, and implementation of the functions to compute the (full and split) stiffness, mass and inverse mass matrices on these meshes.
 
 We have also added codes which generates a locally refined mesh on a unit square such that fine meshes are constructed along rectangular strips. To use these meshes, one needs to make some modifications. These meshes may contain hanging nodes and hence use "Startup_hanging.m" and "Startup_nohanging.m" instead of "Startup_gmsh_square_meshes.m" wherever required. Moreover, store and use matrices from the folder "matrices/rectangular_meshes" instead of "matrices/gmsh_square_meshes". 

