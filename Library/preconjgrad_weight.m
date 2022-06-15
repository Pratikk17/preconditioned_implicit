function [x,iter, delta_m] = preconjgrad_weight(A, M, b, x,L,U,P,Q,tol,maxit)
    r = b - A * x;
    p = L\(P*r);   p = U\p; p = Q*p;   %  pr = M^(-1)*r;
    delta_m = r' * M * p;
    Ap = A * p;
    delta1 = p' * M * Ap;
    alpha = delta_m/delta1;
    d = alpha * p;
    iter = 0;
    while norm(d) > tol && iter<maxit
        iter = iter+1;
        Ap = A * p;
        delta1 = p' * M * Ap;
        alpha = delta_m/delta1;
        d = alpha * p;
        x = x + d;
        r = r - alpha * Ap;
        pr = L\(P*r);   pr = U\pr;    pr = Q*pr;%  pr = M^(-1)*r;
        delta_new = r' * M * pr;
        p = pr + (delta_new/delta_m) * p;
        
        delta_m = delta_new;
        
    end
    
end
