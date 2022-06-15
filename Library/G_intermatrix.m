function [GL GVM GVP] = G_intermatrix()

Globals2D;
GL1 = zeros(3*Nfp, Np);
GL2 = zeros(3*Nfp, Np);
GL3 = zeros(3*Nfp, Np);

GP1 = zeros(Np, Np);  GM1 = zeros(Np, Np);
GP2 = zeros(Np, Np);  GM2 = zeros(Np, Np);
GP3 = zeros(Np, Np);  GM3 = zeros(Np, Np);

% face 1
faceR = r(Fmask(:,1));
GL1(1:Nfp,Fmask(:,1)) = eye(Nfp,Nfp);
GM1(Fmask(:,1),Fmask(:,1)) = eye(Nfp,Nfp);
GP1(Fmask(:,1),Fmask(:,1)) = eye(Nfp,Nfp);

% face 2
faceR = r(Fmask(:,2));
GL2(Nfp+1:2*Nfp,Fmask(:,2)) = eye(Nfp,Nfp);
GM2(Fmask(:,2),Fmask(:,2)) = eye(Nfp,Nfp);
GP2(Fmask(:,2),Fmask(:,2)) = eye(Nfp,Nfp);

% face 3
faceR = r(Fmask(:,3));
GL3(2*Nfp+1:3*Nfp,Fmask(:,3)) = eye(Nfp,Nfp);
GM3(Fmask(:,3),Fmask(:,3)) = eye(Nfp,Nfp);
GP3(Fmask(:,3),Fmask(:,3)) = eye(Nfp,Nfp);

GL={GL1 GL2 GL3};
GVM={GM1 GM2 GM3};
GVP={GP1 GP2 GP3};

return
