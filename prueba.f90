PROGRAM prueba
  
  INTEGER,PARAMETER::num_atoms=7
  INTEGER,PARAMETER::num_crit=3
  INTEGER,DIMENSION(num_atoms)::criterium
  INTEGER,DIMENSION(num_atoms)::orden
  INTEGER,DIMENSION(num_crit,num_atoms)::support

  INTEGER::ii,jj,hh,gg,quedan,criterio
  LOGICAL,DIMENSION(num_atoms)::filtro,filtroaux
  INTEGER,DIMENSION(:),ALLOCATABLE::vect_aux,orden_prov,agugeros
  LOGICAL::interruptor

  criterium(:)=(/1,1,1,0,1,1,1/)
  orden(:)=(/671,672,673,674,675,676,677/)
  support(1,:)=(/6,6,4,4,2,2,1/)
  !support(1,:)=(/6,6,4,4,2,2,2/)
  support(2,:)=(/3,2,4,4,6,6,3/)
  support(3,:)=(/4,4,4,4,3,1,3/)

  !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!
 
  filtro(:)=.FALSE.
  DO ii=1,num_atoms
     IF (criterium(ii)==1) THEN
        filtro(ii)=.TRUE.
     END IF
  END DO
  filtroaux(:)=filtro(:)

  criterio=1
  quedan=COUNT(filtroaux)
  ALLOCATE(orden_prov(quedan),agugeros(quedan),vect_aux(quedan))

  orden_prov(:)=orden(:)

  jj=0
  DO ii=1,num_atoms
     IF (filtro(ii)==.TRUE.) THEN
        jj=jj+1
        agugeros(jj)=ii
     END IF
  END DO



  DO WHILE ((quedan>0) .and. (criterio<=num_crit))
     gg=0
     PRINT*,criterio
     DO WHILE (COUNT(filtroaux)>0)
        hh=MAXLOC(support(criterio,:),DIM=1,MASK=filtroaux)
        filtroaux(hh)=.FALSE.
        gg=gg+1
        vect_aux(gg)=hh
     END DO
     filtroaux(:)=filtro(:)
     interruptor=.TRUE.
     DO ii=1,(quedan-1)
        orden(agugeros(ii))=orden_prov(vect_aux(ii))
        IF (support(criterio,vect_aux(ii))>support(criterio,vect_aux(ii+1))) THEN
           IF (interruptor==.TRUE.) THEN
              filtroaux(agugeros(ii))=.FALSE.
           END IF
           interruptor=.TRUE.
        ELSE
           interruptor=.FALSE.
        END IF
     END DO
     orden(agugeros(quedan))=orden_prov(vect_aux(ii))
     IF (interruptor==.TRUE.) filtroaux(agugeros(quedan))=.FALSE.
     quedan=COUNT(filtroaux)
     criterio=criterio+1
     orden(:)=orden_prov(:)
  END DO

  print'(7I4)', orden(:)

END PROGRAM prueba
