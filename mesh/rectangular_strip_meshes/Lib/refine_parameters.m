function [ref_array]=refine_parameters(mesh_parameters,mesh_parameters_start,ini_elem,K);
% ref_array gives the triangle numbers which needs to be further refined to
% generate the required mesh

switch (mesh_parameters.ratio_x)
    case 1
        if mesh_parameters.Nrefine_coarse<mesh_parameters.Nrefine_fine
            ref_array=[4*(mesh_parameters.Nc^2)*(4^mesh_parameters.Nrefine_coarse)+1:K];
        else
            refine_start=1;
            refine_end=4*(mesh_parameters.Nc^2)*(4^mesh_parameters.Nrefine_fine);
            ref_array=[refine_start:refine_end];
            ref_array=[ref_array ini_elem+1:K];
        end
    case 2
        if mesh_parameters_start.Nrefine_fine<mesh_parameters.Nrefine_fine
            ref_array=[4*(mesh_parameters.Nc^2)*(4^mesh_parameters.Nrefine_coarse)+1:K];
        else
            refine_start=1;
            refine_end=4*(mesh_parameters.Nc^2)*(4^mesh_parameters_start.Nrefine_coarse);
            ref_array=[refine_start:refine_end];
            ref_array=[ref_array ini_elem+1:K];
        end
    case 4
        if mesh_parameters_start.Nrefine_fine<mesh_parameters.Nrefine_fine
            ref_array=[4*(mesh_parameters.Nc^2)*(4^mesh_parameters.Nrefine_coarse)+1:K];
        else
            refine_start=1;
            refine_end=4*(mesh_parameters.Nc^2)*(4^mesh_parameters_start.Nrefine_coarse);
            ref_array=[refine_start:refine_end];
            ref_array=[ref_array ini_elem+1:K];
        end
    otherwise
        disp('Error')
    end
end