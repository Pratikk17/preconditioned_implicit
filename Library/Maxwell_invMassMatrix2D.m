function invM = Maxwell_invMassMatrix2D

% function invM = Maxwell_invMassMatrix2D
% Computes the invers of global Mass matrix of 2D Maxwell problem
% including Jacobian !! (Diference with Hesthaven, Warburton)
% order of variables: [Hx^1; Hy^1; Ez^1; ... Hx^K; Hy^K; Ez^K] = u
% where Hx^k,Hy^k Ez^k are Np x 1 vectors
% Then: du/dt =  invM S u


Globals2D;


% ----------------------------------------------------------------------
% DEFINE "MASS MATRIX" M_global ------------------------------------------
% ----------------------------------------------------------------------


%invM=zeros(3*Np*K);
invM=spalloc(3*Np*K, 3*Np*K, K*9*Np^2);

Zero=zeros(Np);
invM_k = [V*V' Zero Zero; Zero V*V' Zero; Zero Zero V*V'];


for k=1:K
    
    invM( (k-1)*3*Np+1: k*3*Np,  (k-1)*3*Np+1: k*3*Np) = (1/J(1,k))*invM_k;
 
end

%invM=sparse(invM);