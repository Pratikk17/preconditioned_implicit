function [sqrt_MH,sqrt_ME,sqrt_ME_inv]=save_sqrt_matrices_gmsh_square_meshes(M_E,M_H,invM_E,mesh_parameters)

% path for saving matrices (about 660 MB required)
mat_path = 'matrices/gmsh_square_meshes/sqrt/';

Globals2D;

% file name of matrices
mat_name = strcat('sqrt_outer_', mesh_parameters.mesh_level, '_inner_', mesh_parameters.inner_level, ...
    '_polydeg_', num2str(N), '.mat');
mat_file = strcat(mat_path, mat_name);

%
tic;
if exist(mat_file, 'file') == 2
    load(mat_file);
    fprintf('Square root of matrices loaded \n');
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
else
    blockSize = 1 + 2*N + N*(N-1)/2
    sqrt_ME = M_E;
    sqrt_ME_inv = M_E;
    for j = 1:blockSize:length(M_E)-blockSize+1
        sqrt_ME(j:j+blockSize-1,j:j+blockSize-1) = sqrtm(full(M_E(j:j+blockSize-1,j:j+blockSize-1)));
%        inv_sqrt_ME(j:j+blockSize-1,j:j+blockSize-1) = inv(sqrt_ME(j:j+blockSize-1,j:j+blockSize-1));
        sqrt_ME_inv(j:j+blockSize-1,j:j+blockSize-1) = sqrtm(full(invM_E(j:j+blockSize-1,j:j+blockSize-1)));
    end
    sqrt_MH = M_H;
     for j = 1:blockSize:length(M_H)-blockSize+1
         sqrt_MH(j:j+blockSize-1,j:j+blockSize-1) = sqrtm(full(M_H(j:j+blockSize-1,j:j+blockSize-1)));
     end
    save(mat_file, 'sqrt_MH', 'sqrt_ME','sqrt_ME_inv', '-v7.3');
    fprintf('Square root of matrices saved \n');
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
end
fprintf('check: I - inv_sqrt_ME ^2 * M_E = 0    : %e \n',max(max(abs(eye(length(M_E)) - sqrt_ME_inv*sqrt_ME_inv * M_E))));
fprintf('check: inv_sqrt_ME ^2 - inv(full(M_E)) = 0    : %e \n',max(max(abs(sqrt_ME_inv*sqrt_ME_inv - inv(full(M_E))))));
fprintf('check: inv_sqrt_ME ^2 - invM_E = 0    : %e \n',full(max(max(abs(sqrt_ME_inv*sqrt_ME_inv - invM_E)))));
fprintf('check: inv(full(M_E)) - invM_E = 0    : %e \n',max(max(abs(invM_E - inv(full(M_E))))));
% fprintf('difference sqrt M_H  = %e \n', norm(M_H-sqrt_MH^2,'fro')/norm(M_H,'fro'));
% fprintf('difference sqrt M_E  = %e \n', norm(M_E-sqrt_ME^2,'fro')/norm(M_E,'fro'));
% fprintf('difference sqrt invM_E  = %e \n', norm(invM_E-inv_sqrt_ME^2,'fro')/norm(invM_E,'fro'));
end


%     end