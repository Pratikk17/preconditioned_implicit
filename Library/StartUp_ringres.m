function [A, M, invM, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee, B_i, ind_ei_ie_ii] ...
                    = StartUp_ringres (mesh_level, pol_deg, plot_mesh, plot_mat, info)
                
%% Computes the (full and split) stiffness, mass and inverse mass matrices;
% plot_mesh = 1 -> plot mesh with explicit and implicit elements
% plot_mat  = 1 -> spy plot of stiffness matrix and its splitting

%% initialize

% path for saving matrices
%mat_path = '/local/sturm/Maxwell2D_Matrices/central_ringres/';
% path for saving matrices (about 600 MB required per mesh)
path_mass = 'matrices/ringres_mesh/mass';
path_central = 'matrices/ringres_mesh/central';


% mesh path
mesh = strcat('mesh/ringresonator_level_', mesh_level, '.mesh');

% get implicit elements (el_i): fine elements (innner radius<=0.05) and their neighbours
% get impl elements with only impl neigh (el_i_i)
% get expl elements with only expl neigh (el_e_e)
[~, ~, el_i, ~, el_i_i, ~, el_e_e, ~] = fine_coarse_elements (mesh, 0.05, 1e-5, ...
                                                              'in_rad', plot_mesh);

%
Globals2D;
N = pol_deg;

% read in mesh
[Nv, VX, VY, K, EToV] = MeshGen_extension_mesh(mesh);

% construct grid and metric
StartUp2D;

%% Stiffness, mass and inverse mass matrix: A, M, invM

% file name of matrices
mat_file_mass = strcat('matrices/ringres_mesh/mass/mass_mat_level_', mesh_level, '_polydeg_', num2str(N), '.mat');

mat_file_central = strcat('matrices/ringres_mesh/central/central_mat_level_', mesh_level, '_polydeg_', num2str(N), '.mat');
%
tic;
if exist(mat_file_mass, 'file') == 2
    load(mat_file_mass);
    if info, fprintf('Mass matrix loaded \n'); end
    initialization_time = toc;
    if info, fprintf('time to initialize = %8.3e\n', initialization_time); end
else 
    % mass matrix and inverse of mass matrix
    M = Maxwell_MassMatrix2D;
    invM = Maxwell_invMassMatrix2D; 
    % save matrices
    save(mat_file_mass,  'M', 'invM');
    if info, fprintf('Mass matrix saved \n'); end
    initialization_time = toc;
    if info, fprintf('time to initialize = %8.3e\n', initialization_time); end
end

if exist(mat_file_central, 'file') == 2
    load(mat_file_central);
    if info, fprintf('Stiffness matrix loaded \n'); end
    initialization_time = toc;
    if info, fprintf('time to initialize = %8.3e\n', initialization_time); end
else 
% As contains curl part of central fluxes discretization
    % Af contains fluxes part of central fluxes discretization
    [Ss, Sf] = Maxwell_Stif_Matrix2D_nohanging();
    % A contains central fluxes discretization
    A = Ss + Sf;
    save(mat_file_central,  'A');
   if info, fprintf('Stiffness matrix saved \n'); end
end
    
% indices of fields Hx, Hy, Ez; ordering [Hx; Hy; Ez]
a=1:K*3*Np;
a=ceil(a/Np);
a=Modified_mod(a,3);
indHx=find(a==1);
indHy=find(a==2);
indEz=find(a==3);

% Order A, M, invM corresponding to [Hx; Hy; Ez]
A = A([indHx indHy indEz], [indHx indHy indEz]);
M = M([indHx indHy indEz], [indHx indHy indEz]);
invM = invM([indHx indHy indEz], [indHx indHy indEz]);

% indices corresponding to this ordering
indHx = 1:K*Np;
indHy = K*Np+1:2*K*Np;
indEz = 2*K*Np+1:3*K*Np;

% Separate different fields 
% M_H d/dt [Hx; Hy] = -A_E Ez
% M_E d/dt Ez = A_H [Hx; Hy]
% ( note that : A = [0 -A_E; A_H, 0] and A_E = -A_H' )
A_H = A(indEz, [indHx indHy]);
A_E = -A_H';
M_H = M([indHx indHy], [indHx indHy]);
M_E = M(indEz, indEz);
invM_E = invM(indEz, indEz);
invM_H = invM([indHx indHy], [indHx indHy]);


% structure of A
% ___________________________________
% |                    |    ||* * * |
% |                    |    ||* * * | 
% |                  K*Np   ||* * * |
% |                    |    ||* * * |
% |                    |    ||* * * |
% |                  =======||===== |   A_E
% |                         ||* * * |
% |                         ||* * * | 
% |                         ||* * * | 
% |<---K*Np--->||           ||* * * | 
% |            ||           ||* * * |
% |=========================||===== |
% | * * * * * *||* * * * * *||      |
% | * * * * * *||* * * * * *||      |
% | * * * * * *||* * * * * *||      |
% ------------------------------------
%      A_Hx         A_Hy     
%             A_H


%% Decomposition of A, M, invM

% count of implicit elements
count_el_i = length(el_i);
% count of explicit elements
count_el_e = K - count_el_i;


% splitting: [Hx; Hy; Ez] = [Hx^e; Hx^i; Hy^e; Hy^i; Ez^e; Ez^i]
% number of dof in each element is Np = (N+1)(N+2)/2
% ind_i = indices of components of Hx/Hy/Ez treated implcitly
ind_i = zeros(count_el_i*Np,1);
j=0;
for k=el_i'
    j=j+1;
    ind_i((j-1)*Np+1:j*Np, 1)= (k-1)*Np+1:k*Np;
end

% ind_e = indices of components of Ez/Hx/Hy treated explicitly
ind_e = setdiff((1:K*Np)',ind_i);



% A_Hi, A_Ei : implicitly treated parts
% A_Hi, A_Ee : explicitly treated parts
A_Hi = A_H;
A_Hi(:, [ind_e ind_e + K*Np]) = 0;
A_He = A_H;
A_He(:, [ind_i ind_i + K*Np]) = 0;
%
A_Ei = A_E;
A_Ei([ind_e ind_e + K*Np], :) = 0;
A_Ee = A_E;
A_Ee([ind_i ind_i + K*Np], :) = 0;

% decomposition of A
% _____________________________________
% |                     |     ||* * * |
% |                     |     ||* * * |
% |                     |     ||* * * | --------- A_E^e
% |                   K*Np    ||* * * |          |
% |                     |     ||------|          |
% |                     |     ||* * * |          |
% |                     |     ||* * * | ~~~~~~~~~~~~~~~~ A_E^i
% |                   ========||===== |          |      |       
% |                           ||* * * |          |      |
% |                           ||* * * | ---------|      |
% |                           ||* * * |                 |
% |                           ||* * * |                 |
% |                           ||------|                 |
% |<---K*Np---->||            ||* * * | ~~~~~~~~~~~~~~~~|
% |             ||            ||* * * |
% |===========================||===== |
% | * * * * |* *||* * * * |* *||      |
% | * * * * |* *||* * * * |* *||      |
% | * * * * |* *||* * * * |* *||      |
%---------------------------------------
%    A_Hx^e        A_Hy^e     
%           A_Hx^i        A_Hx^i
%       |     '      |       '
%       |     '      |       '
%       ---------------------------------A_H^e
%             '              '
%             '              '
%             ---------------------------A_H^i



if plot_mat
    %
    figure; spy(A);
    title('Stiffness matrix A = [0, A_E; A_H, 0]');
    figure; spy(A_Ei,'r'); hold on; spy(A_Ee,'b');
    title('Split E stiffness matrix A_E; red: implicit, blue: explicit');
    figure; spy(A_Hi,'r'); hold on; spy(A_He,'b');
    title('Split H stiffness matrix A_H; red: implicit, blue: explicit');
    %
end



% count of implicit elements with only implicit neighbours
count_el_i_i = length(el_i_i);
% count of explicit elements with only explicit neigbours
count_el_e_e = K - count_el_i;

% indices of implicit elements with only implicit neighbours
ind_i_i = zeros(count_el_i_i*Np,1);
j=0;
for k=el_i_i'
    j=j+1;
    ind_i_i((j-1)*Np+1:j*Np, 1)= (k-1)*Np+1:k*Np;
end

% indices of implicit elements with explicit neighbours
ind_i_e = setdiff(ind_i,ind_i_i);

%
ind_e_e = zeros(count_el_e_e*Np,1);
j=0;
for k=el_e_e'
    j=j+1;
    ind_e_e((j-1)*Np+1:j*Np, 1)= (k-1)*Np+1:k*Np;
end

%
ind_e_i = setdiff(ind_e,ind_e_e);

%
ind_ei_ie_ii = [ind_e_i; ind_i_e; ind_i_i];


%
A_Hei_ie_y = A_H(ind_e_i, ind_i_e);
A_Hii_ee_y = A_H(ind_i_e, ind_i_e);
A_Hii_ie_y = A_H(ind_i_i, ind_i_e);
A_Hii_ei_y = A_H(ind_i_e, ind_i_i);
A_Hii_ii_y = A_H(ind_i_i, ind_i_i);
%
A_Hei_ie_x = A_H(ind_e_i, ind_i_e + K*Np);
A_Hii_ee_x = A_H(ind_i_e, ind_i_e + K*Np);
A_Hii_ie_x = A_H(ind_i_i, ind_i_e + K*Np);
A_Hii_ei_x = A_H(ind_i_e, ind_i_i + K*Np);
A_Hii_ii_x = A_H(ind_i_i, ind_i_i + K*Np);
%
A_Eie_ei_y = A_E(ind_i_e, ind_e_i);
A_Eii_ee_y = A_E(ind_i_e, ind_i_e);
A_Eii_ei_y = A_E(ind_i_e, ind_i_i);
A_Eii_ie_y = A_E(ind_i_i, ind_i_e);
A_Eii_ii_y = A_E(ind_i_i, ind_i_i);
%
A_Eie_ei_x = A_E(ind_i_e + K*Np, ind_e_i);
A_Eii_ee_x = A_E(ind_i_e + K*Np, ind_i_e);
A_Eii_ei_x = A_E(ind_i_e + K*Np, ind_i_i);
A_Eii_ie_x = A_E(ind_i_i + K*Np, ind_i_e);
A_Eii_ii_x = A_E(ind_i_i + K*Np, ind_i_i);


%
invM_Hii_ee_y = invM_H(ind_i_e, ind_i_e);
invM_Hii_ii_y = invM_H(ind_i_i, ind_i_i);
%
invM_Hii_ee_x = invM_H(ind_i_e + K*Np, ind_i_e + K*Np);
invM_Hii_ii_x = invM_H(ind_i_i + K*Np, ind_i_i + K*Np);
%
invM_Eee_ii = invM_E(ind_e_i, ind_e_i);
invM_Eii_ee = invM_E(ind_i_e, ind_i_e);
invM_Eii_ii = invM_E(ind_i_i, ind_i_i);

%
B_Hi = [invM_Eee_ii*A_Hei_ie_y, sparse(length(ind_e_i),length(ind_i_i)), invM_Eee_ii*A_Hei_ie_x, sparse(length(ind_e_i),length(ind_i_i));...
            invM_Eii_ee*A_Hii_ee_y, invM_Eii_ee*A_Hii_ei_y, invM_Eii_ee*A_Hii_ee_x, invM_Eii_ee*A_Hii_ei_x;...
            invM_Eii_ii*A_Hii_ie_y, invM_Eii_ii*A_Hii_ii_y, invM_Eii_ii*A_Hii_ie_x, invM_Eii_ii*A_Hii_ii_x];
        
B_Ei = [invM_Hii_ee_y*A_Eie_ei_y, invM_Hii_ee_y*A_Eii_ee_y, invM_Hii_ee_y*A_Eii_ei_y; ...
            sparse(length(ind_i_i),length(ind_e_i)), invM_Hii_ii_y*A_Eii_ie_y, invM_Hii_ii_y*A_Eii_ii_y;
            invM_Hii_ee_x*A_Eie_ei_x, invM_Hii_ee_x*A_Eii_ee_x, invM_Hii_ee_x*A_Eii_ei_x; ...
            sparse(length(ind_i_i),length(ind_e_i)), invM_Hii_ii_x*A_Eii_ie_x, invM_Hii_ii_x*A_Eii_ii_x];
        
B_i = B_Hi*B_Ei;
        

end