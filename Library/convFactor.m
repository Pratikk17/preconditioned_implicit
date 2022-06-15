function [cf] = convFactor(rho, gamma, lambda_sq)
%Calculates Convergence Factor by enclosing FOV with a paralelogram and using Sc Toolbox
%lambda_sq is a parameter depending on RK
%rho > 0 is an arbitary value in the preconditioner (I + rho * dt^2 * CE_i CH_i)
%gamma = dt^2*sq_norm_CE_e

    alpha = real(lambda_sq);
    beta = imag(lambda_sq);
    
    a1=1;                        b1=0;
    a2=alpha/rho;                b2=beta/rho;                  
    a3=alpha/rho + alpha*gamma;  b3=beta/rho+beta*gamma; 
    a4=1+alpha*gamma;            b4=beta*gamma;

    p = polygon([a1+b1*1i a4+b4*1i a3+b3*1i a2+b2*1i]);
    %plot(p)
    f = extermap(p);
    cf = abs(evalinv(f,0,10^-12));
    

end
