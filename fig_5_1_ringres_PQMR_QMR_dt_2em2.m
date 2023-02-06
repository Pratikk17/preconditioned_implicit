
% This script generates plots for Fig. 5.4 for Ringresonator. It plots the
% number of QMR and PQMR iterations required at each time step 

clear all;
close all;
format long

%====== Add path to folders containing meshes and routines
addpath(genpath(pwd));

flag_mesh='ringres';
mesh_parameters.mesh_level = '1';                                       % outer mesh level: '1', '2', '3', '4'

%====  Global parameters
Globals2D;
N=4;                                                                    % polynomial degree
alpha_stab=0;                                                           % 0 for central flux 
dt=2e-2;                                                                % time steps tau     
tol=dt.^2;                                                               % tolerance for Krylov subspace methods, here it is for qmr and pqmr
maxit=1000;                                                              % maximum number of iterations for Krylov subspace methods
tol_break=1e-14;                                                        % breakdown of QMR
options.plot_mesh=1;                                                    % To plot mesh
options.info = 1;

%======construction of matrices
[A, M, invM, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee, B_i, ind_ei_ie_ii] ...
        = StartUp_ringres (mesh_parameters.mesh_level, N, options.plot_mesh, 0, options.info);             
       
 A_E  = -A_E;                     % adapt E stiffness matrices to our sign convention 
 A_Ei = -A_Ei;
 A_Ee = -A_Ee;
   
%=== dimension of problem
d_H=size(A_E,1);      d_E=size(A_E,2);
%figure(11); PlotMesh2D();
    
%======= construction of discrete curl operators
CH=invM_E*A_H;         CE=invM_H*A_E;
CH_i=invM_E*A_Hi;      CE_i=invM_H*A_Ei;
CH_e=invM_E*A_He;      CE_e=invM_H*A_Ee;
    
[norm_CE_e,norm_CE_i,sqrt_MH,sqrt_ME,sqrt_ME_inv]=compute_norm_CE_e_CE_i_ringres(mesh_parameters,M_E,M_H,invM_E,CE_e,CE_i);
    
iterations_qmr=[];
iterations_pqmr=[];

alpha=(1/24 + i*sqrt(3)/24)*(dt^2);
gamma=real(alpha);
        
A = speye(size(invM_E)) + alpha*CH*CE ;
B = speye(size(invM_E)) + gamma*CH_i*CE_i;
[B_L,B_U,B_P,B_Q]=lu(B); 
     
%========== building linear system
rng(1);
x0=rand(d_E,1);
b=rand(d_E,1);
    

E_ex=A\b;                                                                  %====exact solution 
[~,x_qmr,iter_qmr] = qmr_weighted(A, b, x0,M_E,tol,maxit,tol_break);
iter_qmr
[err_qmr]=Compute_A_error(x_qmr,E_ex,M_E);
figure(3);
semilogy(err_qmr./err_qmr(1),'-'); hold on;
          
[B1, B2] = lu(B);
[B1_L,B1_U,B1_P,B1_Q] = lu(B1);
[~,x_pqmr,iter_pqmr] = Pqmr_weighted(A,b,x0,M_E,B_L,B_U,B_P,B_Q,B1_L,B1_U,B1_P,B1_Q,tol,maxit,tol_break);
iter_pqmr
[err_pqmr]=Compute_A_error(x_pqmr,E_ex,M_E);

figure(3); 
semilogy(err_pqmr./err_pqmr(1),'--');
 
figure(3); hold on
xlabel('# iterations')
ylabel('relative error')
set(gca,'FontSize',16);
legend('qmr','pqmr','Location','NorthEast')
 %xlim([1 maxit])
% tex_filename=sprintf('tikz_plots/ringres_qmr_pqmr_dt_%s.tex',dt);
% matlab2tikz(tex_filename);
