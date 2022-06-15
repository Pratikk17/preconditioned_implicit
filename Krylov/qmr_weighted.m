function [x,x_iter,iter] = qmr_weighted(A, b, x,M,tol,maxit,tol_break)

% qmr for complex symmetric matrix, A=A'. Algoritm 8.1 of "An implementation
% of the QMR method based on coupled Two-term recurrences" by Roland.
% Freund and Noel Nachtigal. Note that the inner product is with respect to
% Matrix M.


    x_iter=x;
    r0 = b - A * x;      rho=MnormVector(M, r0);  %rho=sqrt(abs(r0'*M*r0)); % initial residue and its M norm
    v=r0/rho;                                                               % v_1: initial basis vector        
    p0=0*v; d0=0*v;
    c0=1; epsi0=1;  neu0=0;  eta0=-1;
    rho0=rho;
    iter=0;
    delta=1;
    d=v;                                                                    % some random d to start with ...
    while norm(d)>tol && iter<maxit                                         % d=x_n - x_{n-1}
        %============= step 1
        delta=v.'*M*v;
        if abs(epsi0)<tol_break || abs(delta)<tol_break
            disp('BREAKDOWN');
            break;
        end
        %======= step 2
        p=v-p0*rho*delta/epsi0;
        
        %==========step 3
        epsi=p.'*M*A*p;          beta=epsi/delta;
        v_tilde=A*p-v*beta;     rho=MnormVector(M, v_tilde);% rho=sqrt(abs(v_tilde.'*M*v_tilde));
        
        %========step 4
        neu=rho/(c0*abs(beta));     c=1/sqrt(1+neu^2);   eta=-eta0*rho0*(c^2)/(beta*(c0^2));
        d=p*eta +d0*((neu0*c)^2);
        x=x+d;
        x_iter=[x_iter x];       
        %=========step 5
        v=v_tilde/rho;
        
        rho0=rho;  d0=d; c0=c; epsi0=epsi;
        neu0=neu; eta0=eta; 
        p0=p;
        iter=iter+1;
    end
    
end
