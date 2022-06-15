function [mesh_para_start,ref_times]=start_mesh(mesh_parameters,K);
% For meshes with hanging nodes, we first need to start with a conforming
% mesh and then refine it to add hanging nodes. We generate
% mes_parameters_start for this starting confirming mesh
% ref_times is the number of times one needs to refine ref_flag triangles

switch (mesh_parameters.ratio_x)
    case 1
        mesh_para_start=mesh_parameters;
        mesh_para_start.Nrefine_coarse=min(mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine);
        mesh_para_start.Nrefine_fine=min(mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine);
        ref_times=abs(mesh_parameters.Nrefine_coarse-mesh_parameters.Nrefine_fine);
    case 2
        if mesh_parameters.Nrefine_coarse==0
            error('Error: First refine coarse mesh atleast once');
%             mesh_para_start=mesh_parameters;
%             mesh_para_start.Nrefine_coarse=min(mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine);
%             mesh_para_start.Nrefine_fine=mesh_para_start.Nrefine_coarse;
%             ref_times=abs(mesh_parameters.Nrefine_fine-mesh_parameters.Nrefine_coarse);
        else
            mesh_para_start=mesh_parameters;
            mesh_para_start.Nrefine_coarse=min(mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine+1);
            mesh_para_start.Nrefine_fine=mesh_para_start.Nrefine_coarse-1;
            if mesh_para_start.Nrefine_coarse<mesh_parameters.Nrefine_coarse
                ref_times=mesh_parameters.Nrefine_coarse-mesh_para_start.Nrefine_coarse;
            else
                ref_times=mesh_parameters.Nrefine_fine-mesh_para_start.Nrefine_fine;
            end
        end
    case 4
        if mesh_parameters.Nrefine_coarse<=1
            error('Error: First refine coarse mesh atleast once');
        else
            mesh_para_start=mesh_parameters;
            mesh_para_start.Nrefine_coarse=min(mesh_parameters.Nrefine_coarse,mesh_parameters.Nrefine_fine+2);
            mesh_para_start.Nrefine_fine=mesh_para_start.Nrefine_coarse-2;
            if mesh_para_start.Nrefine_coarse<mesh_parameters.Nrefine_coarse
                ref_times=mesh_parameters.Nrefine_coarse-mesh_para_start.Nrefine_coarse;
            else
                ref_times=mesh_parameters.Nrefine_fine-mesh_para_start.Nrefine_fine;
            end
        end
    otherwise
        disp('Error')
    end
end
