function [P] = JacobiP(x,alpha1,beta1,N);

% function [P] = JacobiP(x,alpha,beta,N)
% Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
%          (alpha+beta <> -1) at points x for order N and returns P[1:length(xp))]
% Note   : They are normalized to be orthonormal.

% Turn points into row if needed.
xp = x; dims = size(xp);
if (dims(2)==1) xp = xp'; end;

PL = zeros(N+1,length(xp)); 

% Initial values P_0(x) and P_1(x)
gamma0 = 2^(alpha1+beta1+1)/(alpha1+beta1+1)*gamma(alpha1+1)*...
    gamma(beta1+1)/gamma(alpha1+beta1+1);
PL(1,:) = 1.0/sqrt(gamma0);
if (N==0) P=PL'; return; end;
gamma1 = (alpha1+1)*(beta1+1)/(alpha1+beta1+3)*gamma0;
PL(2,:) = ((alpha1+beta1+2)*xp/2 + (alpha1-beta1)/2)/sqrt(gamma1);
if (N==1) P=PL(N+1,:)'; return; end;

% Repeat value in recurrence.
aold = 2/(2+alpha1+beta1)*sqrt((alpha1+1)*(beta1+1)/(alpha1+beta1+3));

% Forward recurrence using the symmetry of the recurrence.
for i=1:N-1
  h1 = 2*i+alpha1+beta1;
  anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha1+beta1)*(i+1+alpha1)*...
      (i+1+beta1)/(h1+1)/(h1+3));
  bnew = - (alpha1^2-beta1^2)/h1/(h1+2);
  PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
  aold =anew;
end;

P = PL(N+1,:)';
return