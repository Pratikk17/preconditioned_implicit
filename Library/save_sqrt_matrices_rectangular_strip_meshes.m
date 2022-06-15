function [sqrt_MH,sqrt_ME,sqrt_ME_inv]=save_sqrt_matrices_rectangular_strip_meshes(M_E,M_H,invM_E,mesh_parameters)

Globals2D;

% file name of matrices

mat_file = strcat('matrices/Pratik/sqrt/sqrt_Nc_', num2str(mesh_parameters.Nc),'_c_', num2str(mesh_parameters.Nrefine_coarse),...
    '_f_', num2str(mesh_parameters.Nrefine_fine),'_r_', num2str(mesh_parameters.ratio_x),'_polydeg_', num2str(N)...
    ,'.mat');

%
tic;
if exist(mat_file, 'file') == 2
    load(mat_file);
    fprintf('sqrt of matrices loaded \n');
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
else
    blockSize = 1 + 2*N + N*(N-1)/2;
    sqrt_ME = M_E;
    sqrt_ME_inv = M_E;
    for j = 1:blockSize:length(M_E)-blockSize+1
        sqrt_ME(j:j+blockSize-1,j:j+blockSize-1) = sqrtm(full(M_E(j:j+blockSize-1,j:j+blockSize-1)));
        %inv_sqrt_ME(j:j+blockSize-1,j:j+blockSize-1) = inv(sqrt_ME(j:j+blockSize-1,j:j+blockSize-1));
        sqrt_ME_inv(j:j+blockSize-1,j:j+blockSize-1) = sqrtm(full(invM_E(j:j+blockSize-1,j:j+blockSize-1)));
    end
    sqrt_MH = M_H;
    for j = 1:blockSize:length(M_H)-blockSize+1
        sqrt_MH(j:j+blockSize-1,j:j+blockSize-1) = sqrtm(full(M_H(j:j+blockSize-1,j:j+blockSize-1)));
    end
    save(mat_file, 'sqrt_MH', 'sqrt_ME','sqrt_ME_inv');
    fprintf('sqrt of matrices saved \n');
    initialization_time = toc;
    fprintf('time to initialize = %8.3e\n', initialization_time);
end

end



    % inv_sqrt_ME = sqrtm(full(invM_E));
    % min(min(invM_E * M_E - inv_sqrt_ME*inv_sqrt_ME * M_E))
    % 
    % min(min(inv_sqrt_ME*inv_sqrt_ME - inv(full(M_E))))


    % return

