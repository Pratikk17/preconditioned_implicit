function [rho] = calcRoh(alpha,gamma, sq_norm_CE_e, sq_norm_CE_i, At_FOV,flag_fov)
    % calculates rho=1/|phi(0)|

    alpha_r=real(alpha);
    alpha_im=imag(alpha);
    
    Gamma_e_g=1+gamma*sq_norm_CE_e;
    Gamma_i_g=1+gamma*sq_norm_CE_i;
    Gamma_e_r=1+alpha_r*sq_norm_CE_e;

    if flag_fov==1
        %====== rho for fov_AT
        p = polygon(At_FOV);
        f = extermap(p);
        g=inv(f);
        rho_fov = abs(evalinv(f,0,10^-12));
    else    
        rho_fov=0;
    end
    %====== rho for quadrilateral Q
    q1r=1;                                             q1i=0;
    q2r=Gamma_e_r;                                     q2i=alpha_im*sq_norm_CE_e;
    q3r=(alpha_r/gamma)*(Gamma_e_g);                   q3i=(alpha_im/gamma)*(Gamma_e_g); 
    q4r=alpha_r/gamma;                                 q4i=alpha_im/gamma;

    Q= polygon([q1r+q1i*1i q2r+q2i*1i q3r+q3i*1i q4r+q4i*1i]);

    fq = extermap(Q);
%    g=inv(fq);
    rho_Q = abs(evalinv(fq,0,10^-12));

    %====== rho for quadrilateral R
    beta= Gamma_e_g-1/Gamma_i_g;
    r1r=1;                                                 r1i=0;                  
    r2r=Gamma_e_g;                                         r2i=0;
    r3r=Gamma_e_g+(alpha_r-gamma)*beta/gamma;              r4i=alpha_im*beta/gamma;
    r4r=1+(alpha_r-gamma)*beta/gamma;   r3i=alpha_im*beta/gamma; 

    R = polygon([r1r+r1i*1i r2r+r2i*1i r3r+r3i*1i r4r+r4i*1i]);
    fr = extermap(R);
   % g=inv(f);
    rho_R = abs(evalinv(fr,0,10^-12));

    rho=[rho_fov; rho_Q;rho_R];
end