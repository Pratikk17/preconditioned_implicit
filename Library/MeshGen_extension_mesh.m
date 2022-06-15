function [Nv, VX, VY, K, EToV] = MeshGen_extension_mesh(FileName)

% function [Nv, VX, VY, K, EToV] = MeshReaderGambit2D(FileName)
% Purpose  : Read in basic grid information to build grid
% NOTE     : .mesh format is assumed



Fid = fopen(FileName, 'rt');
line =fgetl(Fid);

while( strcmp(line,' Vertices')~=1)
    line =fgetl(Fid);
end
Nv_line=fgetl(Fid);

Nv=sscanf(Nv_line, '%d');

VX = (1:Nv); VY = (1:Nv);
for i = 1:Nv
  line = fgetl(Fid);
  tmpx = sscanf(line, '%lf');
  VX(i) = tmpx(1); VY(i) = tmpx(2);
end


while( strcmp(line,' Triangles')~=1)
line =fgetl(Fid);
end

line_K=fgetl(Fid);
K = sscanf(line_K, '%d');


EToV = zeros(K, 3);
for k = 1:K
  line   = fgetl(Fid);
  tmpcon = sscanf(line, '%lf');
  % vertices have to be oriented counter-clockwise in order
  % to keep Jacobian positive (for referent triangle vertices
  % are (-1,-1), (1,-1), (-1,1))
  A = tmpcon(1); 
  B = tmpcon(2);
  C = tmpcon(3);
  
  EToV(k,1) = A;
  EToV(k,2) = B;
  EToV(k,3) = C;
  
  if (VX(A) == VX(B) & VY(A) < VY(B) & VX(C)>VX(A)) || (VX(A) == VX(B) & VY(A) > VY(B) & VX(C)<VX(A)) 
      EToV(k,2) = C;
      EToV(k,3) = B; 
      
  else
  
  a = (VY(B)-VY(A))/(VX(B)-VX(A));
  b = VY(A) - a*VX(A);
  if VX(A) < VX(B) & VY(C) < a*VX(C) + b  % point C is below line AB
      EToV(k,2) = C;
      EToV(k,3) = B; 
  end
  
  if VX(A) > VX(B) & VY(C) > a*VX(C) + b  % point C is above line AB
      EToV(k,2) = C;
      EToV(k,3) = B; 
  end
  
  end
      
end

% Close file
st = fclose(Fid);
return