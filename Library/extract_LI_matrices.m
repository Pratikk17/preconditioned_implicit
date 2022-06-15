function [A_Hi,A_He,A_Ei,A_Ee] =extract_LI_matrices(triangluar_parameters, A_H, A_E,options)
% Build the stiffness matrices for LI scheme proposed by Andreas Stum 

%%======== implicit elements are in the rectangles [i*lc+(i-1)*lf,i*(lc+lf)]x[-1,1] + coarse neighbours of these fine strips
% all remaining elements are explicit
% the numbering of vertices is first vert_e then vert_i
Globals2D();

Nf=triangluar_parameters.Nc-1;  a=triangluar_parameters.a;   b=triangluar_parameters.b;
lc=triangluar_parameters.lc;    lf=triangluar_parameters.lf; meshc_size=triangluar_parameters.meshc_size; 
%%============= indices of vertices of implicit elements
vert_i=[];
epsilon=10*eps;
for j=1:Nf
    v_f=find( (sum(VX(EToV),2)/3 >= a+j*lc+(j-1)*lf) &  (sum(VX(EToV),2)/3 <= a+j*(lc+lf))  );  % elements in fine strip
    v_left=find(( (VX(EToV(:,1)) <= a+j*lc+(j-1)*lf+epsilon & VX(EToV(:,1)) >=a +j*lc+(j-1)*lf-epsilon) & ... % elements in coarse mesh to the left side of fine strip
                  (VX(EToV(:,2)) <= a+j*lc+(j-1)*lf+epsilon & VX(EToV(:,2)) >= a+j*lc+(j-1)*lf-epsilon))|  ...
                 ((VX(EToV(:,1)) <= a+j*lc+(j-1)*lf+epsilon & VX(EToV(:,1)) >= a+j*lc+(j-1)*lf-epsilon) & ...
                  (VX(EToV(:,3)) <= a+j*lc+(j-1)*lf+epsilon & VX(EToV(:,3)) >= a+j*lc+(j-1)*lf-epsilon))| ...
                 ((VX(EToV(:,3)) <= a+j*lc+(j-1)*lf+epsilon & VX(EToV(:,3)) >= a+j*lc+(j-1)*lf-epsilon) & ...
                  (VX(EToV(:,2)) <= a+j*lc+(j-1)*lf+epsilon & VX(EToV(:,2)) >= a+j*lc+(j-1)*lf-epsilon)));
    v_right=find(((VX(EToV(:,1)) <= a+j*(lc+lf)+epsilon & VX(EToV(:,1)) >= a+j*(lc+lf)-epsilon) & ...         % elements in coarse mesh to the left side of fine strip
                  (VX(EToV(:,2)) <= a+j*(lc+lf)+epsilon & VX(EToV(:,2)) >= a+j*(lc+lf)-epsilon))|  ...
                 ((VX(EToV(:,1)) <= a+j*(lc+lf)+epsilon & VX(EToV(:,1)) >= a+j*(lc+lf)-epsilon) & ...
                  (VX(EToV(:,3)) <= a+j*(lc+lf)+epsilon & VX(EToV(:,3)) >= a+j*(lc+lf)-epsilon))| ...
                 ((VX(EToV(:,3)) <= a+j*(lc+lf)+epsilon & VX(EToV(:,3)) >= a+j*(lc+lf)-epsilon) & ...
                  (VX(EToV(:,2)) <= a+j*(lc+lf)+epsilon & VX(EToV(:,2)) >= a+j*(lc+lf)-epsilon)));
    vert_i = [vert_i; v_f ;v_left'; v_right'];
end
vert_i=unique(vert_i);

%%============== indices of vertices of explicit elements
% (K = total number of vertices of all elements)
vert_e = setdiff((1:K)',vert_i);

num_vert_i = length(vert_i);         % number of implicit vertices
num_vert_e = K - num_vert_i;         % number of explicit vertices

%============= Plot implicit (red) and explicit elements (blue)
if options.plot_LI_mesh==1
    figure(41);
    PlotMesh2D();
    hold on
    plot(x(:, vert_i), y(:, vert_i), '*r');
    plot(x(:, vert_e), y(:, vert_e), 'ob');
    tit = sprintf('implicit elements: red, #=%g\nexpicit elements: blue, #=%g', ...
        num_vert_i, num_vert_e);
    title(tit);
    hold off
end

%%========== splitting: [Hx; Hy; Ez] = [Hx^e; Hx^i; Hy^e; Hy^i; Ez^e; Ez^i]
% number of dof in each element is Np = (N+1)(N+2)/2
ind_i = zeros(num_vert_i*Np,1);            % ind_i = indices of components of Hx/Hy/Ez treated implcitly
j=0;
for k=vert_i'
    j=j+1;
    ind_i((j-1)*Np+1:j*Np, 1)= (k-1)*Np+1:k*Np;
end
ind_e = setdiff((1:K*Np)',ind_i);           % ind_e = indices of components of Ez/Hx/Hy treated explicitly


%%=========== Decomposition of A
% A_Hi, A_Ei : implicitly treated parts
% A_He, A_Ee : explicitly treated parts

A_Hi = A_H;
A_Hi(:, [ind_e ind_e + K*Np]) = 0;
A_He = A_H;
A_He(:, [ind_i ind_i + K*Np]) = 0;

A_Ei = A_E;
A_Ei([ind_e ind_e + K*Np], :) = 0;
A_Ee = A_E;
A_Ee([ind_i ind_i + K*Np], :) = 0;

end

