
% This script generates Fig. 5.5. We plot: 
% Error of PQMR for different iterations for fixed time step dt=1e-2.
% Note that maxit is fixed to maxit.


clear all;
close all;
format long
%====== Add path to folders containing meshes and routines
addpath(genpath(pwd));

flag_mesh=2;                                                               % 1 for rectangular_strip_meshes and 2 for gmsh_square_meshes
mesh_parameters.mesh_level = '1';                                          % outer mesh level: '1', '2', '3', '4'
mesh_parameters.inner_level = '1';                                         % inner mesh level: '1', '2', '3', '4'
triangluar_parameters=0;

%====  Global parameters
Globals2D;
N=4;                                                                       % polynomial degree
alpha_stab=0;                                                              % 0 for central flux 
tol=1e-5;                                                                  % tolerance for Krylov subspace methods, here it is for qmr and pqmr
maxit=15;                                                                  % maximum number of iterations for Krylov subspace methods
tol_break=1e-14;                                                            % breakdown of QMR

plotIt = true;
dt = 1e-2;
alpha=(1/24 + i*sqrt(3)/24)*(dt^2);
gamma=real(alpha);

%========= construction of matrices
[~, M, invM, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee] ...
             = StartUp_gmsh_square_meshes(mesh_parameters.mesh_level,  mesh_parameters.inner_level, N, 0,0);
A_E=-A_E; A_Ei=-A_Ei; A_Ee=-A_Ee;                                          % adapt E stiffness matrices to our sign convention

%====== dimension of problem
%figure(11); PlotMesh2D();
d_E=size(A_E,2);

%===construct discrete curl operators
CH=invM_E*A_H;         CE=invM_H*A_E;
CH_i=invM_E*A_Hi;      CE_i=invM_H*A_Ei;
CH_e=invM_E*A_He;      CE_e=invM_H*A_Ee;

[norm_CE_e,norm_CE_i,sqrt_MH,sqrt_ME,sqrt_ME_inv]=compute_norm_CE_e_CE_i...
                      (mesh_parameters,flag_mesh,M_E,M_H,invM_E,CE_e,CE_i);

sq_norm_CE_e=norm_CE_e^2;
sq_norm_CE_i=norm_CE_i^2;

rng(1);
b = 0*rand(length(invM_E),1);
x0 = rand(length(invM_E),1);


A = speye(size(invM_E)) + alpha*CH*CE ;
B = speye(size(invM_E)) + gamma*CH_i*CE_i;
E_ex=A\b;  
   
%=======fov are already saved and available. If not, then first compute them
filename_exist = sprintf('/matrices/gmsh_square_meshes/fov/FOV_At_4RK_gamma_real_alpha_outer_%d_inner_%d_polydeg_%d_dt_%s.mat',...
            str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N,num2str(dt));
filename = sprintf('matrices/gmsh_square_meshes/fov/FOV_At_4RK_gamma_real_alpha_outer_%d_inner_%d_polydeg_%d_dt_%s.mat',...
            str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N,num2str(dt));
if exist(filename_exist,'file')==2
    load([filename])
    disp('FOV loaded');
else
    [L,U,P,Q]=lu(B); B=Q*(U\(L\(P*eye(d_E))));                        % B=inv(Binv) using lu decomposition
    Bsqrt_inv=sqrtm(full(B));
    At=Bsqrt_inv*A*Bsqrt_inv;
    At_M=sqrt_ME*At*sqrt_ME_inv;                                  % equivalent of A_titde in M norm
    At_FOV = wber3(At_M,10,'--b');                                % calculating FOV
    save([filename],'At_FOV')
    disp('FOV Saved');
end

[rho_SC] = calcRoh(alpha,gamma,sq_norm_CE_e, sq_norm_CE_i,At_FOV,1);
phi_fov=1/rho_SC(1);
phi_Q=1/rho_SC(2);
phi_R=1/rho_SC(3);

rho_S=min(rho_SC(2),rho_SC(3));
phi_S=1/rho_S;

[B1, B2] = lu(B);
[B1_L,B1_U,B1_P,B1_Q] = lu(B1);
[B_L,B_U,B_P,B_Q] = lu(B);
[~, ~, iter_pre, errors_pre, errorBounds_const] = Pqmr_weighted_projectionMatrix(A,b,x0,M_E,B_L,B_U,B_P,B_Q,B1_L,B1_U,B1_P,B1_Q,tol,maxit,tol_break, E_ex, sqrt_ME, sqrt_ME_inv);

lin = linspace(1, iter_pre, iter_pre);

errorBounds_fov=errorBounds_const.*min(3./(phi_fov.^lin),2*(phi_fov.^lin)./((phi_fov.^lin)-1));
errorBounds_Q=errorBounds_const.*min(3./(phi_Q.^lin),2*(phi_Q.^lin)./((phi_Q.^lin)-1));
errorBounds_R=errorBounds_const.*min(3./(phi_R.^lin),2*(phi_R.^lin)./((phi_R.^lin)-1));
errorBounds_S=errorBounds_const.*min(3./(phi_S.^lin),2*(phi_S.^lin)./((phi_S.^lin)-1));

figure
semilogy(lin, errors_pre,'r')
xlabel('Number of iterations')
ylabel('Error')
hold on
semilogy(lin, errorBounds_fov,'b','LineWidth',1.5)
%semilogy(lin, errorBounds_Q,'k','LineWidth',1.5)
%semilogy(lin, errorBounds_R,'g','LineWidth',1.5)
semilogy(lin, errorBounds_S,'--m','LineWidth',1.5)

legend('PQMR','Bound-$\mathcal{F}(\widetilde{\mathbf{A}}$', 'Bound-S','Interpreter','latex')     
set(gca,'FontSize',15)
% tex_filename=sprintf('tikz_plots/error_plot_theorem_outer_%d_inner_%d_dt_%s.tex',str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),dt);
% matlab2tikz(tex_filename);

