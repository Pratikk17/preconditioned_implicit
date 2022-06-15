function [n2ed,ed2el]=edge(Elem,Coord)
 % Gives nodes to edge and edge to elements connection
%  Element contains the ordered list of elements present in the mesh
%  Coord contains the ordered list of cordinates present in the mesh

n2el=sparse(size(Coord,1),size(Coord,1));         % size(Coord,1)= Number of nodes in the mesh
for j=1:size(Elem,1)                              % size(Elem,1) = Number of elements in the mesh
    n2el(Elem(j,:),Elem(j,[2 3 1]))=n2el(Elem(j,:),Elem(j,[2 3 1]))+j*eye(3,3);
end
% n2el is a matrix with (n2el)_(i,j) denotes the element to the left of
% edge i to j
B=n2el+n2el';
[I,J]=find(triu(B));                              % Gives the connecting nodes
n2ed=sparse(I,J,1:size(I,1),size(Coord,1),size(Coord,1));   
n2ed=n2ed+n2ed';                                 % if the entry is non zero, it means that i th node shares a edge with jth node

ed2el=zeros(size(I,1),4);                        % Size(I,1) gives number of edges in the mesh
for m = 1:size(Elem,1)
  for k = 1:3 
    p = n2ed(Elem(m,k),Elem(m,rem(k,3)+1)); 
     if ed2el(p,1)==0  
      ed2el(p,:)=[Elem(m,k) Elem(m,rem(k,3)+1) ...
          n2el(Elem(m,k),Elem(m,rem(k,3)+1)) n2el(Elem(m,rem(k,3)+1),Elem(m,k))];                           
    end
  end
end
  


