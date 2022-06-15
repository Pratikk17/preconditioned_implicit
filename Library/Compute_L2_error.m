function [fullerr]=Compute_L2_error(Uout,Hx_ex,Hy_ex,Ez_ex,M)
Globals2D;
% exact solution at FinalTime
Hx_ex = @(x,y) Hx_ex(FinalTime,x,y);
Hy_ex = @(x,y) Hy_ex(FinalTime,x,y);
Ez_ex = @(x,y) Ez_ex(FinalTime,x,y);

% project exact solution
pih_Hx_ex = L2_projection(Hx_ex, 28);
pih_Hy_ex = L2_projection(Hy_ex, 28);
pih_Ez_ex = L2_projection(Ez_ex, 28);

% collect projections in u_ex
pih_u_ex = [pih_Hx_ex(:); pih_Hy_ex(:); pih_Ez_ex(:)];

% compute error
u_h = Uout(end, :)';
fullerr = sqrt( (u_h-pih_u_ex)'*M*(u_h-pih_u_ex) );


end

