
% This script generates plots for Fig. 5.4 which shows the dependence of QMR and
% PQMR error on fine mesh and coarse mesh respectively. 
%



clear all;
close all;
format long
%====== Add path to folders containing meshes and routines
addpath(genpath(pwd));

flag_mesh=2;                                                            % 1 for rectangular_strip_meshes and 2 for gmsh_square_meshes
mesh_parameters.mesh_level = '1';                                       % outer mesh level: '1', '2', '3', '4'

%====  Global parameters
Globals2D;
N=4;                                                                    % polynomial degree
alpha_stab=0;                                                           % 0 for central flux 
dt=5e-2;                                                                % time steps tau     
tol=1e-3;                                                               % tolerance for Krylov subspace methods, here it is for qmr and pqmr
maxit=100;                                                              % maximum number of iterations for Krylov subspace methods
tol_break=1e-14;                                                        % breakdown of QMR

for k=1:4
    mesh_parameters.inner_level = num2str(k)                                      % inner mesh level: '1', '2', '3', '4'
    triangluar_parameters=0;
  
    %======construction of matrices
    [~, ~, invM, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee] ...
                 = StartUp_gmsh_square_meshes(mesh_parameters.mesh_level,  mesh_parameters.inner_level, N, 0, 0);
    A_E  = -A_E;                     % adapt E stiffness matrices to our sign convention 
    A_Ei = -A_Ei;
    A_Ee = -A_Ee;
    
    %=== dimension of problem
    d_H=size(A_E,1);      d_E=size(A_E,2);
    figure(11); PlotMesh2D();
    
    %======= construction of discrete curl operators
    CH=invM_E*A_H;         CE=invM_H*A_E;
    CH_i=invM_E*A_Hi;      CE_i=invM_H*A_Ei;
    CH_e=invM_E*A_He;      CE_e=invM_H*A_Ee;
    
    [norm_CE_e,norm_CE_i,sqrt_MH,sqrt_ME,sqrt_ME_inv]=compute_norm_CE_e_CE_i(mesh_parameters,flag_mesh,M_E,M_H,invM_E,CE_e,CE_i);
    
    iterations_qmr=[];
    iterations_pqmr=[];
    
    alpha=(1/24 + i*sqrt(3)/24)*(dt^2);
    gamma=real(alpha);
        
    A = speye(size(invM_E)) + alpha*CH*CE ;
    B = speye(size(invM_E)) + gamma*CH_i*CE_i;
    [B_L,B_U,B_P,B_Q]=lu(B); 
    Binv=B_Q*(B_U\(B_L\(B_P*eye(d_E))));                                                % Binv=inv(B) using lu decomposition
    Bsqrt_inv=sqrtm(full(Binv));
    At=Bsqrt_inv*A*Bsqrt_inv;
    
    %========== building linear system
    rng(1);
    x0=rand(d_E,1);
    b=rand(d_E,1);
    
    E_ex=A\b;                                                                  %====exact solution 
    [~,x_qmr,iter_qmr] = qmr_weighted(A, b, x0,M_E,tol,maxit,tol_break);
    iter_qmr
    [err_qmr]=Compute_A_error(x_qmr,E_ex,M_E);
    figure(3);
    semilogy(err_qmr./err_qmr(1),'--'); hold on;
         
    [B1, B2] = lu(B);
    [B1_L,B1_U,B1_P,B1_Q] = lu(B1);
    [~,x_pqmr,iter_pqmr] = Pqmr_weighted(A,b,x0,M_E,B_L,B_U,B_P,B_Q,B1_L,B1_U,B1_P,B1_Q,tol,maxit,tol_break);
    iter_pqmr
    [err_pqmr]=Compute_A_error(x_pqmr,E_ex,M_E);

    figure(3); 
    semilogy(err_pqmr./err_pqmr(1),'-');

end

figure(3); hold on
xlabel('# iterations')
ylabel('relative error')
set(gca,'FontSize',16);
legend('qmr-lev-1','pqmr-lev-1','qmr-lev-2','pqmr-lev-2','qmr-lev-3','pqmr-lev-3','qmr-lev-4','pqmr-lev-4','Location','SouthEast')

