function [M,M_H,M_E,A_H,A_E,invM_H,invM_E, A_Hi, A_He, A_Ei, A_Ee] ...
    = StartUp_nohanging(mesh_parameters,options,triangluar_parameters)
% Calculates the mass matrix, stiffness matrix for the Maxwell's equation
% on the considered mesh with no hanging nodes 


% path for saving matrices (about 600 MB required per mesh)
path_mass = 'matrices/rectangular_strip_meshes/mass';
path_central = 'matrices/rectangular_strip_meshes/central';

Globals2D;
mesh_name=generate_mesh_script(mesh_parameters);                           % Generate mesh
[Nv, VX, VY, K, EToV] = MeshReaderGambit2D(mesh_name);                     % read mesh
StartUp2D;                                                                 % construct grid and metric

%====== creating mat file to store the matrices
% file name of matrices
if alpha_stab==0
    mat_file_stiff = strcat('matrices/rectangular_strip_meshes/central/Central_Nc_', num2str(mesh_parameters.Nc),'_c_', num2str(mesh_parameters.Nrefine_coarse),...
    '_f_', num2str(mesh_parameters.Nrefine_fine),'_r_', num2str(mesh_parameters.ratio_x),'_polydeg_', num2str(N)...
    ,'_stiff_', num2str(alpha_stab), '.mat');
else
    mat_file_stiff = strcat('matrices/rectangular_strip_meshes/upwind/Upwind_Nc_', num2str(mesh_parameters.Nc),'_c_', num2str(mesh_parameters.Nrefine_coarse),...
    '_f_', num2str(mesh_parameters.Nrefine_fine),'_r_', num2str(mesh_parameters.ratio_x),'_polydeg_', num2str(N)...
    ,'_stiff_', num2str(alpha_stab), '.mat');
end
mat_file_mass = strcat('matrices/rectangular_strip_meshes/mass/Mass_Nc_', num2str(mesh_parameters.Nc),'_c_', num2str(mesh_parameters.Nrefine_coarse),...
    '_f_', num2str(mesh_parameters.Nrefine_fine),'_r_', num2str(mesh_parameters.ratio_x),'_polydeg_', num2str(N)...
    ,'.mat');

%=========== compute stiffness matrix and mass matrix and its inverse
if exist(mat_file_stiff, 'file') == 2                                      % File exists
    load(mat_file_stiff);
    fprintf('Stiffness matrix loaded \n');
    clear Ss Sf
else   
    [Ss, Sf] = Maxwell_Stif_Matrix2D_nohanging();
    A = Ss + Sf;
    save(mat_file_stiff, 'A');
    fprintf('Stiffness matrix saved \n');
    clear Ss Sf
end
if exist(mat_file_mass, 'file') == 2                                       % File exists
    load(mat_file_mass);
    fprintf('Mass and inverse mass matrices loaded \n');
else   
    M = Maxwell_MassMatrix2D;
    invM = Maxwell_invMassMatrix2D; 
    save(mat_file_mass,'M','invM');
    fprintf('Mass and inverse mass matrices saved \n');
end


%=========== indices of fields Hx, Hy, Ez; ordering [Hx; Hy; Ez]
a=1:K*3*Np;
a=ceil(a/Np);
a=Modified_mod(a,3);
indHx=find(a==1);
indHy=find(a==2);
indEz=find(a==3);

%========== Stiffness, mass and inverse mass matrix: A, M, invM
% Order A, M, invM corresponding to [Hx; Hy; Ez]
A = A([indHx indHy indEz], [indHx indHy indEz]);
M = M([indHx indHy indEz], [indHx indHy indEz]);
invM = invM([indHx indHy indEz], [indHx indHy indEz]);

%=========== indices corresponding to this ordering
indHx = 1:K*Np;
indHy = K*Np+1:2*K*Np;
indEz = 2*K*Np+1:3*K*Np;

% Separating different fields 
% M_H d/dt [Hx; Hy] = A_E Ez
% M_E d/dt Ez = A_H [Hx; Hy]
% ( note that : A = [0 A_E; A_H, 0] and A_E = -A_H' )
A_H = A(indEz, [indHx indHy]);
A_E=-A_H';
M_H = M([indHx indHy], [indHx indHy]);
M_E = M(indEz, indEz);
invM_E = invM(indEz, indEz);
invM_H = invM([indHx indHy], [indHx indHy]);

[A_Hi,A_He,A_Ei,A_Ee]=extract_LI_matrices(triangluar_parameters, A_H, A_E,options);


end
