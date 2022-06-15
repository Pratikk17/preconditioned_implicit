function [Coord,Elem,Db]=redrefine(Coord,Elem,n2ed,ed2el,Db)

nt=size(Elem,1); 
ne=size(ed2el,1);
marker=zeros(ne,1);

for j=1:ne,
    inode=ed2el(j,1); enode=ed2el(j,2);
    coord1=Coord(inode,:); coord2=Coord(enode,:);
    nCoord=(coord1+coord2)/2;
    marker(j)=size(Coord,1)+1; 
    Coord(size(Coord,1)+1,:)=[nCoord(1) nCoord(2)];
end


for j=1:nt,
    ct=Elem(j,:);
    ce=diag(n2ed(ct([2 3 1 ]),ct([3 1 2])));
    m1=marker(ce(1));  m2=marker(ce(2)); m3=marker(ce(3));
    Elem(j,:)=[m1 m2 m3];
    nt1=size(Elem,1); 
    Elem(nt1+1,:)=[ct(1) m3 m2];  
    Elem(nt1+2,:)=[ct(2) m1 m3];  
    Elem(nt1+3,:)=[ct(3) m2 m1];  
end


%%%% Boundary Edges
if (~isempty(Db))
    nb=size(Db,1);
    for j=1:nb,
        base=n2ed(Db(j,1),Db(j,2));
        p=[Db(j,1) marker(base)  Db(j,2)];
        Db(j,:)=[p(1) p(2)];
        Db(size(Db,1)+1,:)=[p(2) p(3)];
    end
end


