%
% This script creates data for dependence of 1/phi_0 (see theorem 4.5)
% on gamma. This data is stored in "optimization_data_fig_5_3.dat",
% and can be plotted to produce fig 5.3 by running "fig_5_3_optimization.m'
%


clear all;
close all;
%====== Add path to folders containing meshes and routines
addpath(genpath(pwd));

%======= for square meshes with fine mesh in center without hanging nodes
flag_mesh=2;  
mesh_parameters.mesh_level = '1';                                       % outer mesh level: '1', '2', '3', '4'
mesh_parameters.inner_level = '1';                                      % inner mesh level: '1', '2', '3', '4'
triangluar_parameters=0;

%====  Global parameters
Globals2D;
N=4;                                                                        % polynomial degree
alpha_stab=0;                                                               % 0 for central flux 
dT=[5e-2];                                                        % time steps tau      

%====== construction of mass and stiffness matrices
[~, ~, ~, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee] ...
             = StartUp_gmsh_square_meshes(mesh_parameters.mesh_level,  ...
               mesh_parameters.inner_level, N, 0, 0);
A_E=-A_E; A_Ei=-A_Ei; A_Ee=-A_Ee;                                           % adapt E stiffness matrices to our sign convention 

% dimension of problem
d_H=size(A_E,1);        d_E=size(A_E,2);
figure(11);             PlotMesh2D();

CH=invM_E*A_H;         CE=invM_H*A_E;
CH_i=invM_E*A_Hi;      CE_i=invM_H*A_Ei;
CH_e=invM_E*A_He;      CE_e=invM_H*A_Ee;

[norm_CE_e,norm_CE_i,sqrt_MH,sqrt_ME,sqrt_ME_inv]=compute_norm_CE_e_CE_i...
                      (mesh_parameters,flag_mesh,M_E,M_H,invM_E,CE_e,CE_i);

sq_norm_CE_e=norm_CE_e^2;
sq_norm_CE_i=norm_CE_i^2;

N_gamma=1e2+1;
phi0_Q=zeros(length(dT),N_gamma);
phi0_R=zeros(length(dT),N_gamma);
phi0_S=zeros(length(dT),N_gamma);

gamma=zeros(length(dT),N_gamma);
gamma_opt=zeros(length(dT),1);
alpha=zeros(length(dT),1);

iter=linspace(1,100,100);

gamma_opt_S_m=zeros(length(dT),length(iter));
phi0_S_opt_m=zeros(length(dT),length(iter));

gamma_opt_R_m=zeros(length(dT),length(iter));
phi0_R_opt_m=zeros(length(dT),length(iter));

gamma_opt_Q_m=zeros(length(dT),length(iter));
phi0_Q_opt_m=zeros(length(dT),length(iter));

for p=1:length(dT)
    dt=dT(p)   
    alpha(p)=(1/24 + i*sqrt(3)/24)*(dt^2);                                 % corresponds to lambda_i^2  of RK4 matrix
%    A = speye(size(invM_E)) + alpha * CH*CE ;                             % A       
    gamma(p,:) = linspace(1e-2*abs(alpha(p)),1*abs(alpha(p)),N_gamma);
 %   B = speye(size(invM_E)) + gamma* CH_i*CE_i;
    for j = 1:length(gamma(p,:))   
        [rho_SC] = calcRoh(alpha(p),gamma(p,j),sq_norm_CE_e, sq_norm_CE_i,0,0);
        rho_S=min(rho_SC(2),rho_SC(3));

        phi0_Q(p,j)=1/rho_SC(2);
        phi0_R(p,j)=1/rho_SC(3);
        phi0_S(p,j)=1/rho_S;
    end
    [min_val,index_g]=min(1./phi0_S(p,:));
    gamma_opt(p)=gamma(p,index_g);
    phi0_S_opt(p)=phi0_S(p,index_g);

    for m=1:length(iter)
        [min_val,index_g]=min(min(3./(phi0_Q(p,:).^m),2*(phi0_Q(p,:).^m)./((phi0_Q(p,:).^m)-1)));
        gamma_opt_Q_m(p,m)=gamma(p,index_g);
        phi0_Q_opt_m(p,m)=phi0_Q(p,index_g);

        [min_val,index_g]=min(min(3./(phi0_R(p,:).^m),2*(phi0_R(p,:).^m)./((phi0_R(p,:).^m)-1)));
        gamma_opt_R_m(p,m)=gamma(p,index_g);
        phi0_R_opt_m(p,m)=phi0_R(p,index_g);

        [min_val,index_g]=min(min(3./(phi0_S(p,:).^m),2*(phi0_S(p,:).^m)./((phi0_S(p,:).^m)-1)));
        gamma_opt_S_m(p,m)=gamma(p,index_g);
        phi0_S_opt_m(p,m)=phi0_S(p,index_g);        
    end

end
filename='optimization_data_fig_5_3';
save(filename,'dT','alpha','sq_norm_CE_e','gamma','iter',...
    'phi0_Q','phi0_R','phi0_S','phi0_S_opt','phi0_Q_opt_m',...
    'phi0_R_opt_m','phi0_S_opt_m','gamma_opt_Q_m','gamma_opt_R_m',...
    'gamma_opt_S_m','gamma_opt')

figure(1); 
plot(gamma(1,:),1./phi0_S(1,:),'r');hold on
xlabel('\gamma'), ylabel('1/\phi_0')
legend('S')
