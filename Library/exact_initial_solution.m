function [Hx_ex,Hy_ex,Ez_ex,J_t,J_X,Hx0,Hy0,Ez0,u0] = exact_initial_solution(options)
Globals2D();
if strcmp(options.InitValue, 'sin')
    %
    omega = sqrt(2)*pi;
    Hx_ex = @(t,x,y) -pi/omega*sin(pi*x).*cos(pi*y) * sin(omega*t);
    Hy_ex = @(t,x,y) pi/omega*cos(pi*x).*sin(pi*y) * sin(omega*t);
    Ez_ex = @(t,x,y) sin(pi*x).*sin(pi*y) * cos(omega*t);
    J_t = @(t) 0;
    J_X = @(x,y) 0;

elseif strcmp(options.InitValue, 'exp')
    %
    Hx_ex = @(t,x,y) -2 * (x-1).*(x+1).*y * exp(t);
    Hy_ex = @(t,x,y)  2 * (y-1).*(y+1).*x * exp(t);
    Ez_ex = @(t,x,y) (x-1).*(x+1).*(y-1).*(y+1) * exp(t);
    J_t = @(t) exp(t);
    J_X = @(x,y) -(x-1).*(x+1).*(y-1).*(y+1) ...
                    + 2*(x-1).*(x+1) + 2*(y-1).*(y+1);
    %
elseif strcmp(options.InitValue, 'sinexp')
    
    Hx_ex = @(t,x,y) -pi * sin(pi*x).*cos(pi*y) * exp(t);
    Hy_ex = @(t,x,y)  pi * cos(pi*x).*sin(pi*y) * exp(t);
    Ez_ex = @(t,x,y) sin(pi*x).*sin(pi*y) * exp(t);
    J_t = @(t) exp(t);
    J_X = @(x,y) (-1-2*pi*pi) * sin(pi*x).*sin(pi*y);
    
else
    %
    error('Error in choice of initial value');
    %
end
% project source term
J_X = L2_projection(J_X, 28);
J_X = J_X(:);

% set initial conditions
Hx0 = @(x,y) Hx_ex(0,x,y);
Hy0 = @(x,y) Hy_ex(0,x,y);
Ez0 = @(x,y) Ez_ex(0,x,y);

Hx0 = L2_projection(Hx0, 28);
Hy0 = L2_projection(Hy0, 28);
Ez0 = L2_projection(Ez0, 28);

u0 = [Hx0(:); Hy0(:); Ez0(:)];


end

