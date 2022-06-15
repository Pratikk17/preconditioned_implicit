function neu_mesh2d_write (neu_filename, Coord, Elem)

% format long
%====  Open the file.
  neu = fopen (neu_filename,'wt' );
 
  %====  Write the data.
  fprintf ( neu, '\t  \t CONTROL INFO 2.2.30\n' );
  fprintf ( neu, '** GAMBIT NEUTRAL FILE \n');  
  fprintf ( neu, '1\n');
  fprintf ( neu, 'PROGRAM:                Gambit     VERSION:  2.2.30 \n');
  fprintf ( neu, '08 June 2020    10:44:37  \n');
  fprintf ( neu, '\t NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL \n');
  fprintf ( neu, '\t %d           %d       %d        %d           %d          %d  \n',size(Coord,1), size(Elem,1), 1,0, 2,2);
  fprintf ( neu, 'ENDOFSECTION \n');
  fprintf ( neu, '\t NODAL COORDINATES 2.2.30 \n');
  Coord_num=1:size(Coord,1)';
  fprintf(neu,'\t \t %d  %25.22f %25.22f   \n',[Coord_num; double(Coord(:,1))'; double(Coord(:,2))']);
  fprintf ( neu, 'ENDOFSECTION \n');
  fprintf ( neu, '      ELEMENTS/CELLS 2.2.30 \n');
  for i = 1 : size(Elem,1)
     fprintf ( neu, '\t \t %d   %d   %d \t  %d \t  %d \t %d\n', 1,3,3, Elem(i,1:3) );
  end
  fprintf ( neu, 'ENDOFSECTION \n');
  fprintf ( neu, '\t \t ELEMENT GROUP 2.2.30 \n');
  fprintf ( neu, 'GROUP:          1 ELEMENTS:          %d MATERIAL:          2 NFLAGS:          1 \n', size(Elem,1));
  fprintf ( neu, '\t \t 0\n');
  strin = num2str(size(Elem,1));
  N_t = str2double(strin(1:end-1));
  N_o = str2double(strin(end));
  for i=0:N_t-1
     fprintf(neu,'\t \t %d   ',i*10+1:i*10+10);
     fprintf(neu,'\n');
  end
  if N_o>0
      fprintf(neu,'\t \t %d   ',N_t*10+1:N_t*10+N_o);
     fprintf(neu,'\n');
  end
  fprintf ( neu, 'ENDOFSECTION \n');
  
end
