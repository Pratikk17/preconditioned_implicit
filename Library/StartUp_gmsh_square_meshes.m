function [A, M, invM, A_H, A_E, M_H, M_E, invM_H, invM_E, A_Hi, A_He, A_Ei, A_Ee] ...
                    = StartUp_gmsh_square_meshes(mesh_level, inner_level, pol_deg, plot_mesh, plot_mat)

%% Computes the (full and split) stiffness, mass and inverse mass matrices;
% plot_mesh = 1 -> plot mesh with explicit and implicit elements
% plot_mat  = 1 -> spy plot of stiffness matrix and its splitting

% path for saving matrices (about 660 MB required)
mat_path = 'matrices/gmsh_square_meshes/mass_stiffness/';

Globals2D;
alpha = 0;
N = pol_deg;

mesh = strcat('mesh/gmsh_square_meshes/mesh_level_', mesh_level, '_inner_level_', inner_level, '.mesh');
[Nv, VX, VY, K, EToV] = MeshGen_extension_mesh(mesh);
   

% construct grid and metric
StartUp2D;

% file name of matrices
mat_name = strcat('mat_outer_', mesh_level, '_inner_', inner_level, ...
    '_polydeg_', num2str(N), '.mat');
mat_file = strcat(mat_path, mat_name);

%
tic;
if exist(mat_file, 'file') == 2
    load(mat_file);
    fprintf('Mass and stiffness matrices loaded \n');
    clear Ss Sf
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
else   
    [Ss, Sf] = Maxwell_Stif_Matrix2D_nohanging();
    A = Ss + Sf;
    M = Maxwell_MassMatrix2D;
    invM = Maxwell_invMassMatrix2D; 
    save(mat_file, 'A', 'M', 'invM', 'Ss', 'Sf');
    fprintf('Mass and stiffness matrices saved \n');
    clear Ss Sf
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
end

% indices of fields Hx, Hy, Ez; ordering [Hx; Hy; Ez]
a=1:K*3*Np;
a=ceil(a/Np);
a=Modified_mod(a,3);
indHx=find(a==1);
indHy=find(a==2);
indEz=find(a==3);

%% Decompose mesh into implicit and explicit part; get indices of 
%% vertices of implicit/ explicit elements

% implicit elements are in the square [-0.05, 0.05]^2
% all remaining elements are explicit
% the numbering of vertices is first vert_e then vert_i

% indices of vertices of implicit elements
vert_i = find( (sum(VX(EToV),2)/3 < 0.05) & ...
               (sum(VX(EToV),2)/3 > -0.05) & ...
               (sum(VY(EToV),2)/3 < 0.05) & ...
               (sum(VY(EToV),2)/3 > -0.05) );

% indices of vertices of explicit elements
% (K = total number of vertices of all elements)
vert_e = setdiff((1:K)',vert_i);

% number of implicit vertices
num_vert_i = length(vert_i);
% number of explicit vertices
num_vert_e = K - num_vert_i;

% Plot implicit (red) and explicit elements (blue)
if plot_mesh
    %
    figure;
    PlotMesh2D();
    hold on
    plot(x(:, vert_i), y(:, vert_i), '.r');
    plot(x(:, vert_e), y(:, vert_e), 'ob');
    tit = sprintf('implicit elements: red, #=%g\nexpicit elements: blue, #=%g', ...
        num_vert_i, num_vert_e);
    title(tit);
    hold off
    %
end

% splitting: [Hx; Hy; Ez] = [Hx^e; Hx^i; Hy^e; Hy^i; Ez^e; Ez^i]
% number of dof in each element is Np = (N+1)(N+2)/2
% ind_i = indices of components of Hx/Hy/Ez treated implcitly
ind_i = zeros(num_vert_i*Np,1);
j=0;
for k=vert_i'
    j=j+1;
    ind_i((j-1)*Np+1:j*Np, 1)= (k-1)*Np+1:k*Np;
end

% ind_e = indices of components of Ez/Hx/Hy treated explicitly
ind_e = setdiff((1:K*Np)',ind_i);

%% Stiffness, mass and inverse mass matrix: A, M, invM

% Order A, M, invM corresponding to [Hx; Hy; Ez]
A = A([indHx indHy indEz], [indHx indHy indEz]);
M = M([indHx indHy indEz], [indHx indHy indEz]);
invM = invM([indHx indHy indEz], [indHx indHy indEz]);

% indices corresponding to this ordering
indHx = 1:K*Np;
indHy = K*Np+1:2*K*Np;
indEz = 2*K*Np+1:3*K*Np;

% Separating different fields 
% M_H d/dt [Hx; Hy] = A_E Ez
% M_E d/dt Ez = A_H [Hx; Hy]
% ( note that : A = [0 A_E; A_H, 0] and A_E = -A_H' )
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
% |                     |     ||* * * | --------- A_E^e
% |                     |     ||* * * |          |
% |                   K*Np    ||------|          |
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
    figure; spy(A_E,'g*'); hold on; spy(A_Ei,'ro'); spy(A_Ee,'bo');
    title('Split E stiffness matrix A_E; red: implicit, blue: explicit');
    figure; spy(A_H,'g*'); hold on; spy(A_Hi,'ro'); spy(A_He,'bo');
    title('Split H stiffness matrix A_H; red: implicit, blue: explicit');
    %
end




end