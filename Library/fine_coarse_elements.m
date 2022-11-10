function [el_f, el_c, el_i, el_e, el_i_i, el_i_e, el_e_e, el_e_i] = fine_coarse_elements (mesh_path, h_f, tol, diam_rad, plot_mesh)

%% partition mesh into fine and coarse part
% diam_rad = 'diam', 'in_rad' -> use diameter or radius of incircle to decide if element
% is fine or coarse
%%%
% returns the global number of the fine elements (diam_rad < h_h +tol) in el_f
% and global number of the coarse elements (remaining elements) in el_c
%%%
% returns the global number of the implicit elements (fine and neighbours)
% in el_i and global number of the explicit elements (remaining elements)
% in el_e
%%%
% number im Sinne von Nummer (in der Nummerierung der Elemente im Gitter)
% und nicht Anzahl der Elemente

%% initialize
Globals2D;
N = 1;

% read in mesh
[Nv, VX, VY, K, EToV] = MeshGen_extension_mesh(mesh_path);

% construct grid and metric
StartUp2D;

%% compute diameters:
% a triangle with side lengths A, B, C has the diameter
% diam = A*B*C/sqrt( (A+B+C)*(B+C-A)*(C+A-B)*(A+B-C) )

A = sqrt( (x(1,:)-x(2,:)).^2 + (y(1,:)-y(2,:)).^2 );
B = sqrt( (x(2,:)-x(3,:)).^2 + (y(2,:)-y(3,:)).^2 );
C = sqrt( (x(3,:)-x(1,:)).^2 + (y(3,:)-y(1,:)).^2 );

diam = (A.*B.*C) ./ sqrt( (A+B+C).*(B+C-A).*(C+A-B).*(A+B-C) );

%% compute radius of in-circle
% in_rad = sqrt( (S-A)*(S-B)*(S-C)/S ) with S = (A+B+C)/2

S = (A+B+C)/2;
in_rad = sqrt( (S-A).*(S-B).*(S-C)./S );

%% get the numbers of the elements

if strcmp(diam_rad, 'diam')
    %
    dr = diam;
    %
elseif strcmp(diam_rad, 'in_rad')
    %
    dr = in_rad;
    %
else
    error('Error in choice of diam_rad. Choose either ''diam'' or ''in_rad'' as geometric factor to sort the elements in fine and coarse.');
end

% numbers of elements
el = (1:K)';

% fine elements (diam<=h_f+tol)
el_f = el(dr<=h_f+tol);

% coarse elements
el_c = setdiff((1:K)',el_f);

if isempty(el_f)
    %
    error('No fine elements for this choice of h_f.');
    %
elseif isempty(el_c)
    %
    error('No coarse elements for this choice of h_f.');
    %
end

% all (fine & coarse) neighbours of fine elements
neigh = EToE(el_f,:);
% sort out fine neighbors
not_el_f = neigh ~= el_f(1);
for k=2:length(el_f)
    %
    not_el_f = not_el_f & (neigh ~= el_f(k)); 
    %
end
% coarse neighbours of fine elements
el_neigh = neigh(not_el_f);

% implicit elements (fine elements and their neighbors)
el_i = [el_neigh; el_f];

% explicit elements
el_e = setdiff((1:K)',el_i);


%% plot
if plot_mesh
    %
    figure;
    %   
    PlotMesh2D();
    hold on
    plot(x(:, el_f), y(:, el_f), '.r');
    plot(x(:, el_c), y(:, el_c), 'ob');
    title('fine elements: red, coarse elements: blue');
    hold off
    %
    figure;
    PlotMesh2D();
    hold on
    plot(x(:, el_i), y(:, el_i), '.r');
    plot(x(:, el_e), y(:, el_e), 'ob');
    title('implicit elements: red, explicit elements: blue');
    hold off
    %
end


% all (implicit & explicit) neighbours of implicit elements
neigh = EToE(el_i,:);
% sort out implicit neighbors
not_el_i = neigh ~= el_i(1);
for k=2:length(el_i)
    %
    not_el_i = not_el_i & (neigh ~= el_i(k)); 
    %
end
% explicit neighbours of implicit elements
el_e_i = neigh(not_el_i);

% explicit elements with only explicit neighbours
el_e_e = setdiff(el_e, el_e_i);


% all (implicit & explicit) neighbours of explicit elements
neigh = EToE(el_e,:);
% sort out explicit neighbors
not_el_e = neigh ~= el_e(1);
for k=2:length(el_e)
    %
    not_el_e = not_el_e & (neigh ~= el_e(k)); 
    %
end
% implicit neighbours of explicit elements
el_i_e = neigh(not_el_e);

% implicit elements with only implicit neighbours
el_i_i = setdiff(el_i, el_i_e);


end