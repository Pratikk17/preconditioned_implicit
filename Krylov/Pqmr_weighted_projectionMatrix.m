function [x,P, iter, errors_pre, errorBounds_const] = Pqmr_weighted_projectionMatrix(A,b,x,M,B_L,B_U,B_P,B_Q,B1_L,B1_U,B1_P,B1_Q,tol,maxit,tol_break, E_ex, sqrt_ME, inv_sqrt_ME)
% preconditioned qmr for complex symmetric matrix, A=A'. Algoritm 8.1 of "An implementation
% of the QMR method based on coupled Two-term recurrences" by Roland.
% Freund and Noel Nachtigal
% We solve Ax=b, with preconditioner B in M-norm which is induced by
% M-innerproduct. The PQMR stops when #iteration is equal to maxit or when
% difference between two consecutive x iterates is less than tol.
% Moreover, the breakdown of PQMR is controlled by tol_breakdown
% Note that [B1_L,B1_U,B1_P,B1_Q]=lu(B1) and [B_L,B_U,B_P,B_Q]=lu(B)


   
    errors_pre = [];
    errorBounds_const = [];
       
    r0 = b - A * x;                                                         % initial residue
    invB1r0=B1_Q*(B1_U\(B1_L\(B1_P*r0)));                                   % invB1r0=inv(B1)*r0; using lu decomposition
    rho= MnormVector(M, invB1r0);
    v=r0/rho;
    
    Mnorm_r0=MnormVector(M, r0);
   
    V = v;
    DIAG = [];
    Diag1=[];
    US = [];    
    LS = [];
    RHOS = [];
    PS=[];
    EPSILONS=[];
    X = [];
   
    p0=0*v; d0=0*v;
    c0=1; epsi0=1;  neu0=0;  eta0=-1;
    rho0=rho;
    iter=0;
   
    
    invBv=B_Q*(B_U\(B_L\(B_P*v)));                                      % invBv=inv(B)*v;                
    delta=v.'*M*invBv;
    DIAG = [DIAG 1/delta];
    Diag1=[Diag1 delta];
    
    d=v;                                                                    % some random d to start with ...
    while iter<maxit                                         % d=x_n - x_{n-1}
        %============= step 1
       
        if abs(epsi0)<tol_break || abs(delta)<tol_break
            disp('BREAKDOWN');
            break;
        end
        %======= step 2
        p=invBv-p0*rho*delta/epsi0;
       
        PS=[PS p];

        US = [US rho*delta/epsi0];
       
        %==========step 3
        epsi=p.'*M*A*p;          beta=epsi/delta;
        v_tilde=A*p-v*beta;      
        invB1v_tilde=B1_Q*(B1_U\(B1_L\(B1_P*v_tilde)));                     % invB1v_tilde=inv(B1)*v_tilde;

        rho=MnormVector(M, invB1v_tilde);
       
        %========step 4
        neu=rho/(c0*abs(beta));     c=1/sqrt(1+neu^2);   eta=-eta0*rho0*(c^2)/(beta*(c0^2));
        d=p*eta +d0*((neu0*c)^2);
        x=x+d;
        X = [X x];
        %=========step 5
        v=v_tilde/rho;
       
        LS = [LS beta];
       
        VM = V; %V_m
        V = [V v]; %V_m+1
       
        RHOS = [RHOS rho0];
       
        rho0=rho;  d0=d; c0=c; epsi0=epsi;
        neu0=neu; eta0=eta;
        p0=p;
        iter=iter+1;
       
        invBv=B_Q*(B_U\(B_L\(B_P*v)));                                      % invBv=inv(B)*v;                
        delta=v.'*M*invBv;
        DIAG = [DIAG 1/delta];
        Diag1=[Diag1 delta];
        m = length(LS);

        U = diag(ones(1,m)) + diag(US(2:m),1) + diag(zeros(1,m-1),-1);
        L = diag(LS) + diag(RHOS(2:m),-1) + diag(zeros(1,m-1),1);

        I = eye(m);

        IL = eye(m);
        l = transpose(rho0 * IL(:, m));
        L = [L; l];
        H = L * U; %H_tilde

        DINV = diag(DIAG);
        D=diag(Diag1);

        I = eye(length(A));

        P = I - V * H * inv(transpose(H) * H) * transpose(H) * DINV * transpose(V)*M;

        error_pre = MnormVector(M, E_ex - x);
        errors_pre = [errors_pre error_pre];

        errorBound_const = MnormMatrix(sqrt_ME, P, inv_sqrt_ME) * (1+sqrt(2)) * Mnorm_r0;
        errorBounds_const = [errorBounds_const errorBound_const];
       
    end
end
