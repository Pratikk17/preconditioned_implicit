function [vmapm_g,vmapm_l,vmapp_g,vmapp_l,nghbor] = vmap_mapping(vmapm,vmapp,vmapb,Nfp,Np,k)

vmapm_g=zeros(3*Nfp,1); vmapm_l=zeros(3*Nfp,1);
vmapp_g=zeros(3*Nfp,1); vmapp_l=zeros(3*Nfp,1);
vmapm_g=vmapm(:,k);
vmapp_g=vmapp(:,k);

vmapm_l=vmapm_g(:,1);
nghbor=zeros(3,1);
for i=1:3*Nfp
    neigh_triangle(i)=ceil(vmapp_g(i)/Np);
    vmapp_l(i)=vmapp_g(i)-Np*(neigh_triangle(i)-1);
end
nghbor=neigh_triangle(1:Nfp:end);

vz=reshape(vmapp_g,Nfp,3);
for i=1:3
    flag_present=0;
    for j=1:size(vmapb,2)
        if(vz(1,i)==vmapb(1,j))
            flag_present=2;
            for l=1:size(vmapb,1)
                if (vz(l,i)==vmapb(l,j))
                    flag_present=1*flag_present;
                else
                    flag_present=0;
                end
            end
        end
        if flag_present==2
            break;
        end
    end
    fl(i)=flag_present;    
end

for i=1:3
    if nghbor(i)==k && fl(i)==0
        nghbor(i)=0;
    end
end

end

