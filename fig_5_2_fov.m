
% These script plots the eigenvalues, region bounding FOV and a
% superset (defined in Theorem 4.3) containing FOV. We consider the 
% fourth order RK implicit time integration scheme on mesh 
% \Thf^{(1)}.  In this case, \alpha is given by eq(5.2).
% The eigenvalues and FOV are already calculated, if not,
% then this code can calculate them. For the preconditioner,
% we choose \gamma=real(alpha)=\tausq/24.
% tau={1e-1,1e-2,1e-3}
%


clear all;
close all;
format long

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
dT=[1e-1;1e-2;1e-3];                                                        % time steps tau      

%=options
options.plot_FOV=1;

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

for p=1:length(dT)
    dt=dT(p)   
    alpha=(1/24 + i*sqrt(3)/24)*(dt^2);
    gamma=real(alpha);

    %======== works only for Andreas' meshes. check run_FOV.m for other meshes
    filename_exist = sprintf('/matrices/gmsh_square_meshes/eig/eig_At_4RK_gamma_real_alpha_outer_%d_inner_%d_polydeg_%d_dt_%s.mat',...
                     str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N,num2str(dt));
    filename = sprintf('matrices/gmsh_square_meshes/eig/eig_At_4RK_gamma_real_alpha_outer_%d_inner_%d_polydeg_%d_dt_%s.mat',...
               str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N,num2str(dt));
    if exist(filename_exist,'file')==2
        load([filename])
        disp('Eigenvalues loaded');
    else
        A = speye(size(invM_E)) + alpha * CH*CE ;                           % A       
        B = speye(size(invM_E)) + gamma* CH_i*CE_i;
        [L,U,P,Q]=lu(B); B=Q*(U\(L\(P*eye(d_E))));                          % B=inv(Binv) using lu decomposition
        Bsqrt_inv=sqrtm(full(B));
        At=Bsqrt_inv*A*Bsqrt_inv;
        At_eigen=eig(At);
        % Note that eigen values of At and At_M=sqrt_ME*At*sqrt_ME_inv are same. Hence we consider only eigenvalues of At.
        save([filename],'At_eigen')
        disp('Eigenvalues Saved');
    end
    %figure; plot(real(At_eigen),imag(At_eigen),'*r'); hold on
    %xlabel('real  axis'); ylabel('imaginary axis')
    
    if options.plot_FOV==1
        filename_exist = sprintf('/matrices/gmsh_square_meshes/fov/FOV_At_4RK_gamma_real_alpha_outer_%d_inner_%d_polydeg_%d_dt_%s.mat',...
            str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N,num2str(dt));
        filename = sprintf('matrices/gmsh_square_meshes/fov/FOV_At_4RK_gamma_real_alpha_outer_%d_inner_%d_polydeg_%d_dt_%s.mat',...
            str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N,num2str(dt));
        if exist(filename_exist,'file')==2
            load([filename])
            disp('FOV loaded');
        else
            A = speye(size(invM_E)) + alpha* CH*CE ;                          % A       
            B = speye(size(invM_E)) + gamma* CH_i*CE_i;
            [L,U,P,Q]=lu(B); B=Q*(U\(L\(P*eye(d_E))));                        % B=inv(Binv) using lu decomposition
            Bsqrt_inv=sqrtm(full(B));
            At=Bsqrt_inv*A*Bsqrt_inv;
            At_M=sqrt_ME*At*sqrt_ME_inv;                                     % equivalent of A_titde in M norm
            At_FOV = wber3(At_M,10,'--b');                                % calculating FOV
            save([filename],'At_FOV')
            disp('FOV Saved');
        end
        plot(real(At_FOV),imag(At_FOV),'--b'); hold on
         xlabel('real  axis'); ylabel('imaginary axis')   
    end

    % ========= Box containing FOV
    alpha_r=real(alpha);    alpha_im=imag(alpha);
    Gamma_e_g=1+gamma*sq_norm_CE_e;
    Gamma_i_g=1+gamma*sq_norm_CE_i;
    Gamma_e_r=1+alpha_r*sq_norm_CE_e;

    if alpha_im ~= 0
    % for quadrilateral Q
        q1r=1;                                             q1i=0;
        q2r=Gamma_e_r;                                     q2i=alpha_im*sq_norm_CE_e;
        q3r=(alpha_r/gamma)*(Gamma_e_g);                   q3i=(alpha_im/gamma)*(Gamma_e_g); 
        q4r=alpha_r/gamma;                                 q4i=alpha_im/gamma;
    
        Q = polyshape([q1r q2r q3r q4r],[q1i q2i q3i q4i]);
        plot(Q)

        % for quadrilateral R
    
        beta= Gamma_e_g-1/Gamma_i_g;
        r1r=1;                                                 r1i=0;                  
        r2r=Gamma_e_g;                                         r2i=0;
        r3r=Gamma_e_g+(alpha_r-gamma)*beta/gamma;              r4i=alpha_im*beta/gamma;
        r4r=1+(alpha_r-gamma)*beta/gamma;   r3i=alpha_im*beta/gamma; 
    
        R = polyshape([r1r r2r r3r r4r],[r1i r2i r3i r4i]);
        plot(R)
    end
end    
