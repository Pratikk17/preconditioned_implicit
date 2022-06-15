function [norm_CE_e,norm_CE_i,sqrt_MH,sqrt_ME,sqrt_ME_inv]...
         =compute_norm_CE_e_CE_i(mesh_parameters,flag_mesh,M_E,M_H,invM_E,CE_e,CE_i)
% flag_mesh=1 for rectangular_strip_meshes and flag_mesh=2 for gmsh_square_meshes
% meshes
Globals2D;
if flag_mesh==1
    filename_CE_e = sprintf('matrices/rectangular_strip_meshes/norm_CEe_CEi/norm_CE_e_Nc_%d_c_%d_r_%d_polydeg_%d.mat',...
        mesh_parameters.Nc,mesh_parameters.Nrefine_coarse,mesh_parameters.ratio_x,N);
    filename_CE_i = sprintf('matrices/rectangular_strip_meshes/norm_CEe_CEi/norm_CE_i_Nc_%d_c_%d_f_%d_r_%d_polydeg_%d.mat',...
        mesh_parameters.Nc,mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine,mesh_parameters.ratio_x,N);
    filename_sqrt = sprintf('matrices/rectangular_strip_meshes/sqrt/sqrt_Nc_%d_c_%d_f_%d_r_%d_polydeg_%d.mat',...
        mesh_parameters.Nc,mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine,mesh_parameters.ratio_x,N);

    file_exist  = sprintf('/matrices/rectangular_strip_meshes/norm_CEe_CEi/norm_CE_i_Nc_%d_c_%d_f_%d_r_%d_polydeg_%d.mat',...
        mesh_parameters.Nc,mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine,mesh_parameters.ratio_x,N);
    file_exist_CEe = sprintf('/matrices/rectangular_strip_meshes/norm_CEe_CEi/norm_CE_e_Nc_%d_c_%d_r_%d_polydeg_%d.mat',...
        mesh_parameters.Nc,mesh_parameters.Nrefine_coarse,mesh_parameters.ratio_x,N);
    
    if exist(file_exist,'file')==2
        load([filename_CE_e])
        disp('norm_CE_e loaded');
        load([filename_CE_i])
        disp('norm_CE_i loaded');
        load([filename_sqrt])        
        disp('sqrt of M_H, M_E and invM_E loaded');
    else
        [sqrt_MH,sqrt_ME,sqrt_ME_inv]=save_sqrt_matrices_rectangular_strip_meshes(M_E,M_H,invM_E,mesh_parameters);
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

if flag_mesh==2
    filename_CE_e = sprintf('matrices/gmsh_square_meshes/norm_CEe_CEi/norm_CE_e_outer_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),N);
    filename_CE_i = sprintf('matrices/gmsh_square_meshes/norm_CEe_CEi/norm_CE_i_outer_%d_inner_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N);
    filename_sqrt = sprintf('matrices/gmsh_square_meshes/sqrt/sqrt_outer_%d_inner_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N);

    file_exist  = sprintf('/matrices/gmsh_square_meshes/norm_CEe_CEi/norm_CE_i_outer_%d_inner_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),str2num(mesh_parameters.inner_level),N);
    file_exist_CEe = sprintf('/matrices/gmsh_square_meshes/norm_CEe_CEi/norm_CE_e_outer_%d_polydeg_%d.mat',...
        str2num(mesh_parameters.mesh_level),N);

    
    if exist(file_exist,'file')==2
        load([filename_CE_e])
        disp('norm_CE_e loaded');
        load([filename_CE_i])
        disp('norm_CE_i loaded');
        load([filename_sqrt])        
        disp('sqrt of M_H, M_E and invM_E loaded');
    else
        [sqrt_MH,sqrt_ME,sqrt_ME_inv]=save_sqrt_matrices_gmsh_square_meshes(M_E,M_H,invM_E,mesh_parameters);
        if exist(file_exist_CEe,'file')==2
            load([filename_CE_e])
            disp('norm_CE_e loaded');
        else
            CE_e_M_norm=sqrt_MH*CE_e*sqrt_ME_inv;
            svds_CE_e=svds(CE_e_M_norm'*CE_e_M_norm,1);
            norm_CE_e=sqrt(svds_CE_e);
            save([filename_CE_e],'norm_CE_e')
            disp('norm_CE_e saved');
        end
        CE_i_M_norm=sqrt_MH*CE_i*sqrt_ME_inv;
        svds_CE_i=svds(CE_i_M_norm'*CE_i_M_norm,1);
        norm_CE_i=sqrt(svds_CE_i);
        save([filename_CE_i],'norm_CE_i')
        disp('norm_CE_i saved');
    end
end