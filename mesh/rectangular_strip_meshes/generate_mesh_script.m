function mesh_name=generate_mesh_script(mesh_parameters)
% we create meshes with and without hanging nodes. Coarse and fine domains are arranged in alternate way, coarse,
% fine, coars, fine,....,fine, coarse. Domain starts and ends with coarse
% subdomains. 
%addpath('Lib_mesh');
tol=1e-13;
L=-2*mesh_parameters.a;                            % length of initial domain
H=-2*mesh_parameters.b;                             % Height of initial domain
Nf=mesh_parameters.Nc-1;                           % Number of fine domain... For our domain, this is always true

[Coord,Elem,Db]=Mesh_generation(mesh_parameters,L,H,Nf,tol);
% figure(11)
% triplot(Elem,Coord(:,1),Coord(:,2));
% drawnow;

%============ writing mesh into fenics or MFEM depending upon flag_package
path='mesh/Pratik_meshes/mesh_files/';
mesh_name=export_mesh(path,mesh_parameters,Coord,Elem);
end


















%================= Intializing the matrices
% K=sparse(size(Coord,1),size(Coord,1));                    % K is global stiffness matrix
% F=sparse(size(Coord,1),1);                                % global load vector
% 
% %===============   Assembly of K 
% % stima is element (local) stiffness matrices
% for i=1:size(Elem,1)
%     K(Elem(i,:),Elem(i,:))=K(Elem(i,:),Elem(i,:))+...
%                                 stima(Coord(Elem(i,:),:));                                                                
% end
% %=============== Assembly of load vector F
% for i=1:size(Elem,1)
%      F(Elem(i,:))=F(Elem(i,:))+det([1 1 1; Coord(Elem(i,:),:)'])*f(sum(Coord(Elem(i,:),:))/3)/6;
% end
% %================= Adding  Dirichlet boundary conditions
% for i=1:size(Db(:,1))
%     j=Db(i,1);
%     K(j,:)=0;
%     K(j,j)=1;
%     F(j,1)=g(Coord(j,:));
% end
% 
%   
