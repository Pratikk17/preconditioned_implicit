function [x,w] = JacobiGQ(alpha1,beta1,N);

% function [x,w] = JacobiGQ(alpha,beta,N)
% Purpose: Compute the N'th order Gauss quadrature points, x, 
%          and weights, w, associated with the Jacobi 
%          polynomial, of type (alpha,beta) > -1 ( <> -0.5).

if (N==0) x(1)= -(alpha1-beta1)/(alpha1+beta1+2); w(1) = 2; return; end;

% Form symmetric matrix from recurrence.
J = zeros(N+1);
h1 = 2*(0:N)+alpha1+beta1;
J = diag(-1/2*(alpha1^2-beta1^2)./(h1+2)./h1) + ...
    diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha1+beta1).*...
    ((1:N)+alpha1).*((1:N)+beta1)./(h1(1:N)+1)./(h1(1:N)+3)),1);
if (alpha1+beta1<10*eps) J(1,1)=0.0;end;
J = J + J';

%Compute quadrature by eigenvalue solve
[V,D] = eig(J); x = diag(D);
w = (V(1,:)').^2*2^(alpha1+beta1+1)/(alpha1+beta1+1)*gamma(alpha1+1)*...
    gamma(beta1+1)/gamma(alpha1+beta1+1);
return;
