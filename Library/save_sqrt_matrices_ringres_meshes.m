function [sqrt_MH,sqrt_ME,sqrt_ME_inv]=save_sqrt_matrices_ringres_meshes(M_E,M_H,invM_E,mesh_parameters)

% path for saving matrices (about 660 MB required)
mat_path = 'matrices/ringres_mesh/sqrt/';

Globals2D;

% file name of matrices
mat_name = strcat('sqrt_mesh_level_', mesh_parameters.mesh_level,'_polydeg_', num2str(N), '.mat');
mat_file = strcat(mat_path, mat_name);

%
tic;
if exist(mat_file, 'file') == 2
    load(mat_file);
    fprintf('Square root of matrices loaded \n');
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
else
    blockSize = 1 + 2*N + N*(N-1)/2;
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

end
