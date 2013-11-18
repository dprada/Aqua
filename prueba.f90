PROGRAM prueba
  
  INTEGER,PARAMETER::num_atoms=7
  INTEGER,PARAMETER::num_crit=3
  INTEGER,DIMENSION(num_atoms)::criterium
  INTEGER,DIMENSION(num_atoms)::orden
  INTEGER,DIMENSION(num_crit,num_atoms)::support

  INTEGER::ii
  INTEGER::base
  INTEGER,DIMENSION(num_crit)::potencias
  INTEGER,DIMENSION(num_atoms)::valores

  criterium(:)=(/1,1,1,0,1,1,1/)
  orden(:)    =(/1,2,3,4,5,6,7/)
  support(1,:)=(/6,6,4,4,2,2,1/)
  support(2,:)=(/3,2,4,4,6,6,3/)
  support(3,:)=(/4,4,4,4,3,1,3/)
  !support(1,:)=(/6,6,4,4,2,2,2/)

  !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!

  base=MAXVAL(support)+1
  
  potencias(num_crit)=1
  DO ii=(num_crit-1),1,-1
     potencias(ii)=potencias(ii+1)*base
  END DO

  DO ii=1,num_atoms
     valores(ii)=dot_product(potencias,support(:,ii))
     print*, valores(ii)
  END DO

END PROGRAM prueba
