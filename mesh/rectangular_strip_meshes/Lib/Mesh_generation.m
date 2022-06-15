function [Coord,Elem,Db]=Mesh_generation(mesh_parameters,L,H,Nf,tol) 
% Give back cordinates, elements and Diriclet boundary
% NOTE that coarse and fine domains are arranged in alternate way, coarse,
% fine, coars, fine,....,fine, coarse. Domain starts and ends with coarse
% subdomains
% Since the fine mesh has very small length in x direction, we first divide fine mesh
% into Nc squares and then each square into ratio_x squares
% So each fine mesh is divided into ratio_x*Nc*4 triangles
% The initial domain is square [-a,a]x[-b,b]
a=mesh_parameters.a;
b=mesh_parameters.b;
ratio_x=mesh_parameters.ratio_x;
Nc=mesh_parameters.Nc;
Nrefine_c=mesh_parameters.Nrefine_coarse;
Nrefine_f=mesh_parameters.Nrefine_fine;

lc=ratio_x*L/(Nc-1+ratio_x*Nc);                     % length of coarse domain      ...... solving L=Nc*lc+(Lc-1)*lf , with lf=lc/ratio_x 
lf=lc/ratio_x;                                  % length of fine domain    
hf=(H/Nc)/ratio_x;                              % height of fine mesh 
num_c=(3*Nc+2);                             % number of cordinates in one coarse subdomain

%====================== Coarse mesh Coordinate assembly of initial mesh
count_node=0;       
for i=1:Nc                              % For cordinates of Coarse domain
    Coord(count_node+1,:)=[(i-1)*(lc+lf)+a 0+b];
    Coord(count_node+2,:)=[i*lc+(i-1)*lf+a 0+b];
    count_node=count_node+2;
    for j=1:Nc-1
        Coord(count_node+1,:)=[i*lc+(i-1)*lf+a j*H/Nc+b];
        Coord(count_node+2,:)=[(i-1)*(lc+lf)+a j*H/Nc+b];
        Coord(count_node+3,:)=[0.5*lc*(2*i-1)+(i-1)*lf+a 0.5*(2*j-1)*H/Nc+b];
        count_node=count_node+3;
    end
    j=Nc;
    Coord(count_node+1,:)=[i*lc+(i-1)*lf+a H+b];
    Coord(count_node+2,:)=[(i-1)*(lc+lf)+a H+b];
    Coord(count_node+3,:)=[0.5*lc*(2*i-1)+(i-1)*lf+a 0.5*(2*j-1)*H/Nc+b];
    count_node=count_node+3;
end
%====================== Coarse mesh Element assembly of initial mesh
count_elem=0;
for i=1:Nc                               % For elements of Coarse domain
    Elem(count_elem+1,:)=[(i-1)*num_c+1  (i-1)*num_c+2  (i-1)*num_c+5];
    Elem(count_elem+2,:)=[(i-1)*num_c+2  (i-1)*num_c+3  (i-1)*num_c+5];
    Elem(count_elem+3,:)=[(i-1)*num_c+3  (i-1)*num_c+4  (i-1)*num_c+5];
    Elem(count_elem+4,:)=[(i-1)*num_c+4  (i-1)*num_c+1  (i-1)*num_c+5];
    count_elem=count_elem+4;
    for j=2:Nc
        Elem(count_elem+1,:)=[(i-1)*num_c+3*j-2  (i-1)*num_c+3*j-3  (i-1)*num_c+3*j+2];
        Elem(count_elem+2,:)=[(i-1)*num_c+3*j-3  (i-1)*num_c+3*j    (i-1)*num_c+3*j+2];
        Elem(count_elem+3,:)=[(i-1)*num_c+3*j    (i-1)*num_c+3*j+1  (i-1)*num_c+3*j+2];
        Elem(count_elem+4,:)=[(i-1)*num_c+3*j+1  (i-1)*num_c+3*j-2  (i-1)*num_c+3*j+2];
        count_elem=count_elem+4;
    end
end

%====================== Coarse mesh Dirichlet boundary assembly of initial mesh
count_Db=0;
for j=Nc:-1:1             % on y direction
    Db(count_Db+1,:)=[3*j+1 3*j-2];                       % Left boundary
    count_Db=count_Db+1;
end
for i=1:Nc
    Db(count_Db+1,:)=[(i-1)*num_c+1     (i-1)*num_c+2];
    Db(count_Db+2,:)=[(i-1)*num_c+3*Nc  (i-1)*num_c+3*Nc+1];
    count_Db=count_Db+2;
end
Db(count_Db+1,:)=[(Nc-1)*num_c+3*j-1 (Nc-1)*num_c+3*j];                       % Right boundary
count_Db=count_Db+1;
for j=2:Nc              % in y direction
    Db(count_Db+1,:)=[(Nc-1)*num_c+3*(j-1) (Nc-1)*num_c+3*j];                       % Right boundary
    count_Db=count_Db+1;
end

%======================== First refining coarse mesh by Nrefine_c times
for p=1:Nrefine_c                                 % Number of times you want to refine the initial coarse mesh
    %================= Edge-Node-Element Connections
    [n2ed,ed2el]=edge(Elem,Coord);               % Gives nodes to edge and edge to elements connection
   
    %================== Element Redrefine
    [Coord,Elem,Db]=redrefine(Coord,Elem,n2ed,ed2el,Db); 
end
gl_count_node=size(Coord,1);
gl_count_elem=size(Elem,1);
gl_count_Db=size(Db,1);

%========== creating local elements and Db of fine mesh. it will be same
%for all fine mesh strips.  We will map them to global coordinates afterwards.

Coord_f0=zeros(3*Nc*ratio_x+2,2);        % Coordinates of each fine strip before refinement
Elem_f0=zeros(4*Nc*ratio_x,3);          % Elements of each fine strip before refinement
Db_f0=zeros(ratio_x*Nc*2,2);            % BOundary elements of each fine strip  before refinement
map_local_global=zeros(size(Coord_f0,1),1);   % map local fine coordinates to global
Db_f0(1,:)=[1 2];
Db_f0(Nc*ratio_x+1,:)=[3*ratio_x*Nc 3*ratio_x*Nc+1];
count_Db=1;
count_elem=0;                           % local counter of elements of each strip
Elem_f0(1,:)=[1  2  5]; Elem_f0(2,:)=[2  3  5];
Elem_f0(3,:)=[3  4  5]; Elem_f0(4,:)=[4  1  5];
count_elem=count_elem+4;
for j=2:Nc*ratio_x
    Elem_f0(count_elem+1,:)=[3*j-2  3*j-3  3*j+2];
    Elem_f0(count_elem+2,:)=[3*j-3  3*j    3*j+2];
    Elem_f0(count_elem+3,:)=[3*j    3*j+1  3*j+2];
    Elem_f0(count_elem+4,:)=[3*j+1  3*j-2  3*j+2];
    count_elem=count_elem+4;

    Db_f0(count_Db+1,:)=[3*j-3 3*j];
    Db_f0(Nc*ratio_x+count_Db+1,:)=[3*ratio_x*Nc+1-3*(j-1) 3*ratio_x*Nc+1-3*j];
    count_Db=count_Db+1;
end

for i=1:Nf                                  % Creating local coordinates for fine mesh
    Coord_f0=0*Coord_f0;
    Coord_f0(1:2,:)=[i*lc+(i-1)*lf+a 0+b; i*(lc+lf)+a 0+b];
    count_coord=2;
    for j=1:Nc*ratio_x
        Coord_f0(count_coord+1,:)=[i*(lc+lf)+a             j*hf+b];
        Coord_f0(count_coord+2,:)=[i*lc+(i-1)*lf+a         j*hf+b];
        Coord_f0(count_coord+3,:)=[i*lc+0.5*(2*i-1)*lf+a   0.5*(2*j-1)*hf+b];
        count_coord=count_coord+3;
    end
    Coord_f=Coord_f0;
    Elem_f=Elem_f0;                    % required for refinement
    Db_f=Db_f0;
    %======================== refining i^th fine strip mesh by Nrefine_f times
    for p=1:Nrefine_f                                 % Number of times you want to refine the initial mesh
        %================= Edge-Node-Element Connections
        [n2ed,ed2el]=edge(Elem_f,Coord_f);               % Gives nodes to edge and edge to elements connection
        %================== Element Redrefine
        [Coord_f,Elem_f,Db_f]=redrefine(Coord_f,Elem_f,n2ed,ed2el,Db_f); 
    end
    %============= Map local coordinates of fine mesh to global coordinates
    %Ncoord_c is the number of total coarse coordinates
    strip=find(abs(Coord(:,1)-(i*lc+(i-1)*lf+a))<tol | abs(Coord(:,1)-(i*(lc+lf)+a))<tol); % Finding already existing coordinates on strip
    flag_a=0;                                % flag to know if coordinate needs to be added or already present in mesh
    for j=1:size(Coord_f,1)
        for k=1:length(strip)
            if abs(Coord(strip(k),1)-Coord_f(j,1))<tol & abs(Coord(strip(k),2)-Coord_f(j,2))<tol          
                map_local_global(j)=strip(k); % if coordinates are already present, then mapping them
                flag_a=1;
            end
        end
        if flag_a==0
            map_local_global(j)=gl_count_node+1;    % creating new coordinates
            Coord(gl_count_node+1,:)=Coord_f(j,:);
            gl_count_node=gl_count_node+1;
        end
        flag_a=0;
    end

    %============= Global mapping of fine elements
    ne=size(Elem_f,1);                 % number of elements in refines fine strip
    Elem(gl_count_elem+1:gl_count_elem+ne,:)= ...
        [map_local_global(Elem_f(1:ne,1))  map_local_global(Elem_f(1:ne,2))  map_local_global(Elem_f(1:ne,3))];     

    %======== mapping of global boundary elements
    strip_Db_b=find((Coord_f(:,2)==0+b) & (Coord_f(:,1)>=i*lc+(i-1)*lf+a) & (Coord_f(:,1)<=i*(lc+lf)+a));     % Coordinates of fine mesh on bottom boundary.
    [~,Ib]=sort(Coord_f(strip_Db_b,1));    % sorting them in asscending order, with Ib giving the map
    strip_Db_t=find((Coord_f(:,2)==H+b) & (Coord_f(:,1)>=i*lc+(i-1)*lf+a) & (Coord_f(:,1)<=i*(lc+lf)+a));     % Coordinates of fine mesh on top boundary.
    [~,It]=sort(Coord_f(strip_Db_t,1));     % sorting them in asscending order, with It giving the map
    for j=1:length(Ib)-1
        Db(gl_count_Db+j,:)=[map_local_global(strip_Db_b(Ib(j))) map_local_global(strip_Db_b(Ib(j+1)))];
        Db(gl_count_Db+(2^Nrefine_f)+j,:)=[map_local_global(strip_Db_t(It(end-j+1))) map_local_global(strip_Db_t(It(end-j)))];
    end
    gl_count_elem=size(Elem,1);
    gl_count_Db=size(Db,1);
end
end
