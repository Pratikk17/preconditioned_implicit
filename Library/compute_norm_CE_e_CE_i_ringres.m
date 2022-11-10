function [norm_CE_e,norm_CE_i,sqrt_MH,sqrt_ME,sqrt_ME_inv]...
         =compute_norm_CE_e_CE_i_ringres(mesh_parameters,M_E,M_H,invM_E,CE_e,CE_i)

Globals2D;

filename_CE_e = sprintf('matrices/ringres_mesh/norm_CEe_CEi/norm_CE_e_mesh_level_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),N);
filename_CE_i = sprintf('matrices/ringres_mesh/norm_CEe_CEi/norm_CE_i_mesh_level_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),N);
filename_sqrt = sprintf('matrices/ringres_mesh/sqrt/sqrt_mesh_level_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),N);

file_exist  = sprintf('/matrices/ringres_mesh/norm_CEe_CEi/norm_CE_i_mesh_level_%d_polydeg_%d.mat',...
       str2num(mesh_parameters.mesh_level),N);
file_exist_CEe = sprintf('/matrices/ringres_mesh/norm_CEe_CEi/norm_CE_e_mesh_level_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),N);
    
if exist(file_exist,'file')==2
    load([filename_CE_e])
    disp('norm_CE_e loaded');
    load([filename_CE_i])
    disp('norm_CE_i loaded');
    load([filename_sqrt])        
    disp('sqrt of M_H, M_E and invM_E loaded');
else
    [sqrt_MH,sqrt_ME,sqrt_ME_inv]=save_sqrt_matrices_ringres_meshes(M_E,M_H,invM_E,mesh_parameters);
     if exist(file_exist_CEe,'file')==2
         load([filename_CE_e])
            disp('norm_CE_e loaded');
        else
            %norm_CE_e = norm(sqrt_MH*CE_e*sqrt_ME_inv);
            CE_e_M_norm=sqrt_MH*CE_e*sqrt_ME_inv;
            svds_CE_e=svds(CE_e_M_norm'*CE_e_M_norm,1);
            norm_CE_e=sqrt(svds_CE_e);
            save([filename_CE_e],'norm_CE_e')
            disp('norm_CE_e saved');
        end

        norm_CE_i = MnormMatrix(sqrt_MH, CE_i, sqrt_ME_inv);
        save([filename_CE_i],'norm_CE_i')
        disp('norm_CE_i saved');
    end
end