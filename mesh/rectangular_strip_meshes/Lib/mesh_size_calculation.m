function [triangluar_parameters]=mesh_size_calculation(mesh_parameters)
% calculate the coarse mesh size (hc), fine mesh size(hf) and hmin and hmax

triangluar_parameters=mesh_parameters;

L=-2*mesh_parameters.a;                            % length of initial domain
H=-2*mesh_parameters.b;                             % Height of initial domain
Nf=mesh_parameters.Nc-1;                           % Number of fine domain... For our domain, this is always true
ratio_x=mesh_parameters.ratio_x;
Nc=mesh_parameters.Nc;
Nrefine_c=mesh_parameters.Nrefine_coarse;
Nrefine_f=mesh_parameters.Nrefine_fine;

triangluar_parameters.lc=ratio_x*L/(Nc-1+ratio_x*Nc);                     % length of coarse domain (before refinement)     ...... solving L=Nc*lc+(Lc-1)*lf , with lf=lc/ratio_x 
hc=H/Nc;                                                                  % height of fine mesh   (before refinement)
triangluar_parameters.lf=triangluar_parameters.lc/ratio_x;                % length of fine domain (before refinement)   
hf=hc/ratio_x;                                                             % height of fine mesh   (before refinement) 

mesh_lc=triangluar_parameters.lc/(2^Nrefine_c);
mesh_hc=hc/(2^Nrefine_c);
mesh_lf=triangluar_parameters.lf/(2^Nrefine_f);
mesh_hf=hf/(2^Nrefine_f);

triangluar_parameters.meshc_size=(mesh_hc+mesh_lc)/2;
triangluar_parameters.meshf_size=(mesh_hf+mesh_lf)/2;

triangluar_parameters.hmin=min(triangluar_parameters.meshc_size,triangluar_parameters.meshf_size);
triangluar_parameters.hmax=max(triangluar_parameters.meshc_size,triangluar_parameters.meshf_size);




end

