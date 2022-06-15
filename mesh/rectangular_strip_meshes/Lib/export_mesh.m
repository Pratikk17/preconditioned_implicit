function mesh_name=export_mesh(path,mesh_parameters,Coord,Elem)
% this function exports matlab mesh into .neu format so that it can be read
% for DG discretization and LI
    
%========== writing into .neu file ( in order to import into matlab suitable for codes using Hestaven book)
filename=sprintf('%s%s_%d_%s_%d_%s_%d_%s_%d.%s',path,'Nc',mesh_parameters.Nc,'c',...
    mesh_parameters.Nrefine_coarse,'f',mesh_parameters.Nrefine_fine,'r',mesh_parameters.ratio_x,'neu');
neu_mesh2d_write(filename,Coord, Elem)        %2 stands for space dimension and 3 stands for triangles

mesh_name=sprintf('%s_%d_%s_%d_%s_%d_%s_%d.%s','Nc',mesh_parameters.Nc,'c',...
    mesh_parameters.Nrefine_coarse,'f',mesh_parameters.Nrefine_fine,'r',mesh_parameters.ratio_x,'neu');
end


