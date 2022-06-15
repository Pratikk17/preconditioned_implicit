function flag_h=hanging_node_flag(mesh_parameters)
% flag_h=0 if there are no hanging nodes in the mesh and 1 if there is
% atleast one hanging node in mesh
switch (mesh_parameters.ratio_x)
    case 1
        if mesh_parameters.Nrefine_coarse==mesh_parameters.Nrefine_fine
            flag_h=0;
        else
            flag_h=1;
        end
    case 2
        if mesh_parameters.Nrefine_coarse==mesh_parameters.Nrefine_fine+1
            flag_h=0;
        else
            flag_h=1;
        end
    case 4
        if mesh_parameters.Nrefine_coarse==mesh_parameters.Nrefine_fine+2
            flag_h=0;
        else
            flag_h=1;
        end
    otherwise
        disp('Error')
    end
end