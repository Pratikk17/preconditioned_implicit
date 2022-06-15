function [Ss, Sf] = Maxwell_Stif_Matrix2D_hanging(alpha_stab)
% Creates the Stiffness matrix for Maxwell's equation in TM mode.

Globals2D;

neighbors = BuildHNonCon2D(N+1, 1e-6);    % find non-conforming neighbors
savemapB = mapB; savevmapB = vmapB; mapB = mapW;  vmapB = vmapW; % 1.0 Only apply PEC boundary conditions at wall boundaries 

%% ========= DEFINE "STIFNESS MATRIX" Ss 
Ss=sparse(3*Np*K, 3*Np*K);
Zero=zeros(Np);
Dr1=MassMatrix*Dr_tilda;
Ds1=MassMatrix*Ds_tilda;

for k=1:K
    Dx_k=rx(1,k)*Dr1+sx(1,k)*Ds1;
    Dy_k=ry(1,k)*Dr1+sy(1,k)*Ds1;
    Ss_k=[Zero Zero Dy_k; Zero Zero -Dx_k; Dy_k -Dx_k, Zero];
   
    Ss( (k-1)*3*Np+1: k*3*Np,  (k-1)*3*Np+1: k*3*Np) = J(1,k)*Ss_k;
end

%% ========= DEFINE  "FLUX MATRIX" Sf 
Sf=sparse(3*Np*K, 3*Np*K);

[GL GVM GVP] = G_intermatrix();
nbd_faces=size(vmapB,1)/Nfp;
vmapM_e=reshape(vmapM, Nfp*Nfaces, K);
vmapP_e=reshape(vmapP, Nfp*Nfaces, K);
vmapb=reshape(vmapB,Nfp,nbd_faces);
for k=1:K                                    %loop over all elements
    nx_k=[nx(1,k); nx(Nfp+1,k); nx(2*Nfp+1,k)];            % outward normals on element k
    ny_k=[ny(1,k); ny(Nfp+1,k); ny(2*Nfp+1,k)];
 
    % indices of neighbouring elements to element k
    [vmapM_global_k,vmapM_local_k,vmapP_global_k,vmapP_local,nghbor]=vmap_mapping(vmapM_e,vmapP_e,vmapb,Nfp,Np,k);
    
    Jac=J(1,k);                         % jacobian indices
    for l=1:Nfaces                      %loop over all faces
        if nghbor(l)==k      %boundary          
            fluxEz_Hx = 2*(-ny_k(l)*Fscale((l-1)*Nfp+1,k))*GVM{l}; 
            fluxEz_Hy = 2*(nx_k(l)*Fscale((l-1)*Nfp+1,k))*GVM{l};
            fluxEz_Ez = 2*(-alpha_stab*Fscale((l-1)*Nfp+1,k))*GVM{l};

            rEz_Hx=Jac*MassMatrix*LIFT*GL{l}*fluxEz_Hx/2;
            rEz_Hy=Jac*MassMatrix*LIFT*GL{l}*fluxEz_Hy/2;
            rEz_Ez=Jac*MassMatrix*LIFT*GL{l}*fluxEz_Ez/2;
            
            Sf_k=[rEz_Hx rEz_Hy rEz_Ez];
            Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:k*3*Np)= Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:k*3*Np)+Sf_k; 
             
        elseif nghbor(l)==0                            %values from inside of element k--leave them unchanged. will be updated afterwards
        else                                           %values from inside of element k
            fluxEz_Hx = (-ny_k(l)*Fscale((l-1)*Nfp+1,k))*GVM{l};
            fluxEz_Hy = (nx_k(l)*Fscale((l-1)*Nfp+1,k))*GVM{l};
            fluxEz_Ez = (-alpha_stab*Fscale((l-1)*Nfp+1,k))*GVM{l};

            rEz_Hx=Jac*MassMatrix*LIFT*GL{l}*fluxEz_Hx/2;
            rEz_Hy=Jac*MassMatrix*LIFT*GL{l}*fluxEz_Hy/2;
            rEz_Ez=Jac*MassMatrix*LIFT*GL{l}*fluxEz_Ez/2;

            Sf_k=[rEz_Hx rEz_Hy rEz_Ez];
            Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:k*3*Np)= Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:k*3*Np)+Sf_k;
            %values from outside of element k
            GVP=zeros(Np);
            indx=vmapP_local((l-1)*Nfp+1:l*Nfp);
            GVP(Fmask(:,l),indx)=eye(Nfp,Nfp);                            

            fluxEzn_Hx = (-ny_k(l)*Fscale((l-1)*Nfp+1,k))*GVP; 
            fluxEzn_Hy = (nx_k(l)*Fscale((l-1)*Nfp+1,k))*GVP; 
            fluxEzn_Ez = (-alpha_stab*Fscale((l-1)*Nfp+1,k))*GVP;

            rEzn_Hx=Jac*MassMatrix*LIFT*GL{l}*fluxEzn_Hx/2;
            rEzn_Hy=Jac*MassMatrix*LIFT*GL{l}*fluxEzn_Hy/2;
            rEzn_Ez=Jac*MassMatrix*LIFT*GL{l}*fluxEzn_Ez/2;

            Sfn_k=[rEzn_Hx rEzn_Hy rEzn_Ez];
            Sf((k-1)*3*Np+2*Np+1: k*3*Np, (nghbor(l)-1)*3*Np+1:nghbor(l)*3*Np)=Sfn_k;            
        end           
    end

end

mapB = savemapB; vmapB = savevmapB;    % 1.2 Restore original boundary node lists

%%===========  2.0 Correct lifted fluxes at each non-conforming face fragment
Nnoncon = length(neighbors);
for n=1:Nnoncon
    neigh = neighbors{n};
    %2.1 Extract information about this non-conforming face fragment
    k1 = neigh.elmtM;  GVM = neigh.gVM;    k2 = neigh.elmtP;  GVP = neigh.gVP;
    lnx = neigh.nx;    lny = neigh.ny;     l1=neigh.faceM;    l2=neigh.faceP;
    lift1=neigh.lift;

    Jac=J(1,k1);
   % For Ez       %  values from inside of element k
    fluxEz_Hx = -lny*GVM;        fluxEz_Hy = lnx.*GVM;     fluxEz_Ez = -alpha_stab*GVM;

    rEz_Hx=Jac*MassMatrix*lift1*fluxEz_Hx/2;
    rEz_Hy=Jac*MassMatrix*lift1*fluxEz_Hy/2;
    rEz_Ez=Jac*MassMatrix*lift1*fluxEz_Ez/2;

    Sf_k=[rEz_Hx rEz_Hy rEz_Ez];
    Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+1:k1*3*Np)= Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+1:k1*3*Np)+Sf_k; 
  %  values from outside of element k
    fluxEzn_Hx = -lny*GVP;    fluxEzn_Hy = lnx.*GVP;     fluxEzn_Ez = -alpha_stab*GVP;
    
    rEzn_Hx=Jac*MassMatrix*lift1*fluxEzn_Hx/2;
    rEzn_Hy=Jac*MassMatrix*lift1*fluxEzn_Hy/2;
    rEzn_Ez=Jac*MassMatrix*lift1*fluxEzn_Ez/2;

    Sfn_k=[rEzn_Hx rEzn_Hy rEzn_Ez];
    Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+1:k2*3*Np)=...
                Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+1:k2*3*Np)+Sfn_k; 
end 

end


%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:(k-1)*3*Np+Np) = ...
%                 Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:(k-1)*3*Np+Np)+rEz_Hx;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+Np+1:(k-1)*3*Np+2*Np) = ...
%                 Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+Np+1:(k-1)*3*Np+2*Np)+rEz_Hy;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+2*Np+1:k*3*Np) = ...
%                  Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+2*Np+1:k*3*Np)+rEz_Ez;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:(k-1)*3*Np+Np) = ...
%                 Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+1:(k-1)*3*Np+Np)+rEz_Hx;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+Np+1:(k-1)*3*Np+2*Np) = ...
%                 Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+Np+1:(k-1)*3*Np+2*Np)+rEz_Hy;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+2*Np+1:k*3*Np) =  ...
%                 Sf((k-1)*3*Np+2*Np+1: k*3*Np, (k-1)*3*Np+2*Np+1:k*3*Np)+rEz_Ez;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (nghbor(l)-1)*3*Np+1:(nghbor(l)-1)*3*Np+Np) = rEzn_Hx;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (nghbor(l)-1)*3*Np+Np+1:(nghbor(l)-1)*3*Np+2*Np) =rEzn_Hy;
%             Sf((k-1)*3*Np+2*Np+1: k*3*Np, (nghbor(l)-1)*3*Np+2*Np+1:nghbor(l)*3*Np) =  rEzn_Ez;

%     Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+1:(k1-1)*3*Np+Np) = ...
%         Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+1:(k1-1)*3*Np+Np)+rEz_Hx;
%     Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+Np+1:(k1-1)*3*Np+2*Np) = ...
%         Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+Np+1:(k1-1)*3*Np+2*Np)+rEz_Hy;
%     Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+2*Np+1:k1*3*Np) = ...
%          Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k1-1)*3*Np+2*Np+1:k1*3*Np)+rEz_Ez;

%     Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+1:(k2-1)*3*Np+Np) = ...
%         Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+1:(k2-1)*3*Np+Np)+rEzn_Hx;
%     Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+Np+1:(k2-1)*3*Np+2*Np) = ...
%         Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+Np+1:(k2-1)*3*Np+2*Np)+rEzn_Hy;
%     Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+2*Np+1:k2*3*Np) = ...
%          Sf((k1-1)*3*Np+2*Np+1: k1*3*Np, (k2-1)*3*Np+2*Np+1:k2*3*Np)+rEzn_Ez;