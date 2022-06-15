function [phi_q,phi_r] = convFactor_S(gamma,alpha,sq_norm_CE_i,sq_norm_CE_e)
% Calculates phi_0 (in theorem 4.5) by enclosing FOV with a quadrilateral S
% and using Sc Toolbox.
%    A = speye(size(invM_E)) + alpha * CH*CE ;                           % A    
%    B = speye(size(invM_E)) + gamma* CH_i*CE_i;     

alpha_r=real(alpha);
alpha_im=imag(alpha);

Gamma_e_g=1+gamma*sq_norm_CE_e;
Gamma_i_g=1+gamma*sq_norm_CE_i;
Gamma_e_r=1+alpha_r*sq_norm_CE_e;

% for quadrilateral Q
    q1r=1;                                             q1i=0;
    q2r=Gamma_e_r;                                     q2i=alpha_im*sq_norm_CE_e;
    q3r=(alpha_r/gamma)*(Gamma_e_g);                   q3i=(alpha_im/gamma)*(Gamma_e_g); 
    q4r=alpha_r/gamma;                                 q4i=alpha_im/gamma;

    %Q = polyshape([q1r q2r q3r q4r],[q1i q2i q3i q4i]);
    Q= polygon([q1r+q1i*1i q2r+q2i*1i q3r+q3i*1i q4r+q4i*1i]);

    %figure(1);plot(Q)
    fq = extermap(Q);
    cf_q = abs(evalinv(fq,0,10^-12)); 
    
    % for quadrilateral R

    beta= Gamma_e_g-1/Gamma_i_g;
    r1r=1;                                                 r1i=0;                  
    r2r=Gamma_e_g;                                         r2i=0;
    r3r=Gamma_e_g+(alpha_r-gamma)*beta/gamma;              r4i=alpha_im*beta/gamma;
    r4r=1+(alpha_r-gamma)*beta/gamma;   r3i=alpha_im*beta/gamma; 

    %R = polyshape([r1r r2r r3r r4r],[r1i r2i r3i r4i]);
    R = polygon([r1r+r1i*1i r2r+r2i*1i r3r+r3i*1i r4r+r4i*1i]);
    %plot(R)
    fr = extermap(R);
    cf_r = abs(evalinv(fr,0,10^-12));
    
    phi_q=1/cf_q;
    phi_r=1/cf_r;
    

end
