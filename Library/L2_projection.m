function u_proj = L2_projection (u, Corder)
% L2 projection of u onto V_h = space of piecewise polynomials of degree N

Globals2D; 

% cubature nodes, weights and number of nodes
[cub.R,cub.S,cub.W, cub.Ncub] = Cubature2D(Corder);

%
va = EToV(:,1)'; vb = EToV(:,2)'; vc = EToV(:,3)';
cub.X = 0.5*(-(cub.R+cub.S)*VX(va)+(1+cub.R)*VX(vb)+(1+cub.S)*VX(vc));
cub.Y = 0.5*(-(cub.R+cub.S)*VY(va)+(1+cub.R)*VY(vb)+(1+cub.S)*VY(vc));

[a,b] = rstoab(cub.R,cub.S);

u_proj = zeros(N,K);

% compute inner product of u with the Np = 0.5*(N+1)*(N+2) orthonormal
% polynomials psi_m
for i = 0:N
    for j = 0:N-i
        m = j + (N+1)*i +1 - 0.5*i*(i-1);
        
        psi = Simplex2DP(a,b,i,j);
       
        u_proj_2 = u(cub.X,cub.Y) .* (psi*ones(K,1)');
        u_proj_2 = J(1,:) .* (cub.W' * u_proj_2);
        
        u_proj(m,:) = u_proj_2;
        
        clear u_proj_2;
    end
end

% change from psi basis to lagrange polynomial basis via Vandermonde matrix
V = Vandermonde2D(N,r,s);

u_proj =  V * (J.^(-1) .* u_proj);
