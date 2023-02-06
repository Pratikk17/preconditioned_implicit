%========= This is the main script file 

clear all;
close all;
format long
%====== Add path to folders containing meshes and routines
addpath(genpath(pwd));

flag_mesh=2;                                                                % 1 for rectangular_strip_meshes and 2 for gmsh_square_meshes
mesh_parameters.mesh_level = '1';                                           % outer mesh level: '1', '2', '3', '4'
mesh_parameters.inner_level = '1';                                          % inner mesh level: '1', '2', '3', '4'
triangluar_parameters=0;

%====  Global parameters
Globals2D;
N=8;                                                                        % polynomial degree
alpha_stab=0;                                                               % 0 for central flux 
FinalTime=2;                                                                % final time
dT=[1e-2;5e-3;2e-3;1e-3];                               % time steps tau     
%dT=[1e-3];
%tol=1e-5;                                                                   % tolerance for Krylov subspace methods, here it is for qmr and pqmr
maxit=1000;                                                                  % maximum number of iterations for Krylov subspace methods
tol_break=1e-14;                                                            % breakdown of QMR

%======= Options regarding the time integrator, initial values, plots,
options.InitValue = 'sin';                                                  % initial value: 'sin'- Homogoneos,'sinexp'-with source term
options.N_save =10;                                                         % save solution after N.save time steps 
options.plot_solution =1;                                                   % if 1: plot solution every N_save time steps

[A, M, invM, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee] ...
             = StartUp_gmsh_square_meshes(mesh_parameters.mesh_level,  mesh_parameters.inner_level, N, 0, 0);
A_E=-A_E; A_Ei=-A_Ei; A_Ee=-A_Ee;                                           % adapt E stiffness matrices to our sign convention 

% define exact solutions and  set initial conditions according to options.InitValue
[Hx_ex,Hy_ex,Ez_ex,J_t,J_X,~,~,~,u0] = exact_initial_solution(options);

% dimension of problem
dim_prb = size(A_E,1)+size(A_E,2);
% figure(11)
% PlotMesh2D();

time_RK4=zeros(length(dT),1);
time_RK4_PQMR=zeros(length(dT),1);
time_DIRK4=zeros(length(dT),1);
time_DIRK4_PQMR=zeros(length(dT),1);

errvec_RK4=zeros(length(dT),1);
errvec_RK4_PQMR=zeros(length(dT),1);
errvec_DIRK4=zeros(length(dT),1);
errvec_DIRK4_PQMR=zeros(length(dT),1);

for p=1:length(dT)
    dt=dT(p)
    %tol=1e-5;
    tol=dt^4;

    if options.plot_solution==0                    % solution is saved only at finaltime
        options.N_save=0; 
    end

    disp('Implicit RK4');
    tic;              
    [Tout_RK4,Uout_RK4] = RK4_implicit(A_H,A_E,invM_H,invM_E, J_t, J_X,FinalTime, u0, dt, options.N_save);
    time_RK4(p)=toc;
    errvec_RK4(p)=Compute_L2_error(Uout_RK4,Hx_ex,Hy_ex,Ez_ex,M);
    [errvec_RK4(p) time_RK4(p)]
   %Plot_Solution(Tout_RK4,Uout_RK4,options.plot_solution,dim_prb);        % plot solutio
  
    disp('Implicit RK4_PQMR');
    tic;              
    [Tout_RK4_PQMR,Uout_RK4_PQMR,iter_PQMR] = RK4_PQMR(A_H,A_E,A_Hi,A_Ei,M_E,invM_H,invM_E,J_t, J_X,...
        FinalTime,u0,dt,options.N_save,tol,tol_break,maxit);
    time_RK4_PQMR(p)=toc;
    errvec_RK4_PQMR(p)=Compute_L2_error(Uout_RK4_PQMR,Hx_ex,Hy_ex,Ez_ex,M);
    [errvec_RK4_PQMR(p) time_RK4_PQMR(p) iter_PQMR]
    %Plot_Solution(Tout_RK4_PQMR,Uout_RK4_PQMR,options.plot_solution,dim_prb);        % plot solution

    disp('Diagonally Implicit RK4');
    tic;              
    [Tout_DIRK4,Uout_DIRK4] = DIRK4(A_H,A_E,invM_H,invM_E, J_t, J_X,FinalTime, u0, dt, options.N_save);
    time_DIRK4(p)=toc;
    errvec_DIRK4(p)=Compute_L2_error(Uout_DIRK4,Hx_ex,Hy_ex,Ez_ex,M);
    [errvec_DIRK4(p) time_DIRK4(p)]
   %Plot_Solution(Tout_DIRK4,Uout_DIRK4,options.plot_solution,dim_prb);        % plot solutio

    disp('DIRK4_PQMR');
    tic;              
    [Tout_DIRK4_PQMR,Uout_DIRK4_PQMR,iter_DIRK_PQMR] = DIRK4_PQMR(A_H,A_E,A_Hi,A_Ei,M_E,invM_H,invM_E,J_t, J_X,...
        FinalTime,u0,dt,options.N_save,tol,tol_break,maxit);
    time_DIRK4_PQMR(p)=toc;
    errvec_DIRK4_PQMR(p)=Compute_L2_error(Uout_DIRK4_PQMR,Hx_ex,Hy_ex,Ez_ex,M);
    errvec_DIRK4_PQMR(p)    
    [errvec_DIRK4_PQMR(p) time_DIRK4_PQMR(p) iter_DIRK_PQMR(p)]
    %Plot_Solution(Tout_DIRK4_PQMR,Uout_DIRK4_PQMR,options.plot_solution,dim_prb);        % plot solution
end

figure; loglog(dT,errvec_RK4,'*-r','linewidth', 1.5); hold on;
loglog(dT,errvec_RK4_PQMR,'d--c','linewidth', 1.5);
loglog(dT,errvec_DIRK4,'s-k','linewidth', 1.5); hold on;
loglog(dT,errvec_DIRK4_PQMR,'d--g','linewidth', 1.5);
loglog(dT,10*dT.^4,'--b','linewidth', 1.5);
ylim([10^(-9) 10^(-4)])
xlabel('time step (\tau)');
ylabel('error');
set(gca,'FontSize',15)
legend('implicit RK4','RK4-PQMR','DIRK4','DIRK4-PQMR','location','NorthWest');
%legend('implicit RK4','DIRK4','location','NorthWest');

tex_filename=sprintf('tikz_plots/time_integration_error.tex');
 matlab2tikz(tex_filename);



 figure; loglog(dT,time_RK4,'*-r','linewidth', 1.5); hold on;
loglog(dT,time_RK4_PQMR,'d--c','linewidth', 1.5);
loglog(dT,time_DIRK4,'s-k','linewidth', 1.5); hold on;
loglog(dT,time_DIRK4_PQMR,'d--g','linewidth', 1.5);
%loglog(dT,10*dT.^4,'--b','linewidth', 1.5);
%ylim([10^(-9) 10^(-4)])
xlabel('time step (\tau)');
ylabel('time (s)');
set(gca,'FontSize',15)
legend('implicit RK4','RK4-PQMR','DIRK4','DIRK4-PQMR','location','NorthWest');
%legend('implicit RK4','DIRK4','location','NorthWest');

tex_filename=sprintf('tikz_plots/time_integration_time.tex');
 matlab2tikz(tex_filename);