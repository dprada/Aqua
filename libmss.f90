MODULE GLOB

INTEGER::definition_hbs

CONTAINS

SUBROUTINE ind_wat_limit_4_nosim (mss,aux,hbs,hbdists,num_wats,num_atoms,num_hbs)


  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms,num_wats,num_hbs
  INTEGER,DIMENSION(num_atoms,3),INTENT(IN)::aux
  INTEGER,DIMENSION(num_hbs,3),INTENT(IN)::hbs
  DOUBLE PRECISION,DIMENSION(num_hbs),INTENT(IN)::hbdists
  INTEGER,DIMENSION(num_wats,17),INTENT(OUT)::mss
  
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_first_shell
  INTEGER,DIMENSION(:,:),ALLOCATABLE::first_shell
  INTEGER,DIMENSION(:),ALLOCATABLE::bonds_o
  INTEGER,DIMENSION(:),ALLOCATABLE::microstate
  
  INTEGER::ii,jj,kk,bb,gg
  INTEGER::ind_oh,hi,ind_o
  
  mss=0
  ALLOCATE(dist_first_shell(num_wats,4),first_shell(num_wats,4),bonds_o(num_wats))
  ALLOCATE(microstate(17))

  first_shell=0
  dist_first_shell=0.0d0
  bonds_o=0


  IF (ANY((/2,3/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        ind_oh=aux(hbs(ii,1)+1,1)+1
        hi=aux(hbs(ii,2)+1,2)
        ind_o=aux(hbs(ii,3)+1,1)+1
        IF (first_shell(ind_oh,hi)==0) THEN
           first_shell(ind_oh,hi)=ind_o
           dist_first_shell(ind_oh,hi)=hbdists(ii)
        ELSE
           IF (dist_first_shell(ind_oh,hi)>hbdists(ii)) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           END IF
        END IF
        IF (bonds_o(ind_o)==0) THEN
           first_shell(ind_o,3)=ind_oh
           dist_first_shell(ind_o,3)=hbdists(ii)
           bonds_o(ind_o)=1
        ELSE
           IF (bonds_o(ind_o)==1) THEN
              IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                 first_shell(ind_o,4)=ind_oh
                 dist_first_shell(ind_o,4)=hbdists(ii)
              ELSE
                 first_shell(ind_o,4)=first_shell(ind_o,3)
                 dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                 first_shell(ind_o,3)=ind_oh
                 dist_first_shell(ind_o,3)=hbdists(ii)
              END IF
              bonds_o(ind_o)=2
           ELSE
              IF (bonds_o(ind_o)==2) THEN
                 IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 ELSE
                    IF (dist_first_shell(ind_o,4)>hbdists(ii)) THEN
                       first_shell(ind_o,4)=ind_oh
                       dist_first_shell(ind_o,4)=hbdists(ii)
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  ELSE IF (ANY((/1/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        ind_oh=aux(hbs(ii,1)+1,1)+1
        hi=aux(hbs(ii,2)+1,2)
        ind_o=aux(hbs(ii,3)+1,1)+1
        IF (first_shell(ind_oh,hi)==0) THEN
           first_shell(ind_oh,hi)=ind_o
           dist_first_shell(ind_oh,hi)=hbdists(ii)
        ELSE
           IF (dist_first_shell(ind_oh,hi)<hbdists(ii)) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           END IF
        END IF
        IF (bonds_o(ind_o)==0) THEN
           first_shell(ind_o,3)=ind_oh
           dist_first_shell(ind_o,3)=hbdists(ii)
           bonds_o(ind_o)=1
        ELSE
           IF (bonds_o(ind_o)==1) THEN
              IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                 first_shell(ind_o,4)=ind_oh
                 dist_first_shell(ind_o,4)=hbdists(ii)
              ELSE
                 first_shell(ind_o,4)=first_shell(ind_o,3)
                 dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                 first_shell(ind_o,3)=ind_oh
                 dist_first_shell(ind_o,3)=hbdists(ii)
              END IF
              bonds_o(ind_o)=2
           ELSE
              IF (bonds_o(ind_o)==2) THEN
                 IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 ELSE
                    IF (dist_first_shell(ind_o,4)<hbdists(ii)) THEN
                       first_shell(ind_o,4)=ind_oh
                       dist_first_shell(ind_o,4)=hbdists(ii)
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  END IF

  DO kk=1,num_wats

     microstate=0
     microstate(1)=kk
     
     microstate(2:5)=(/ first_shell(kk,:) /)

     gg=6
     DO ii=2,3
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg:gg+1)=(/ first_shell(bb,1:2) /)
           IF ((first_shell(bb,3)==kk).or.(first_shell(bb,4)==kk)) THEN
              IF (first_shell(bb,3)==kk) THEN
                 microstate(gg+2)=first_shell(bb,4)
              ELSE
                 microstate(gg+2)=first_shell(bb,3)
              END IF
              IF ((first_shell(bb,3)==kk).and.(first_shell(bb,4)==kk)) THEN
                 !print*,'problema ssn1' !!! Comprobar esto
                 microstate(gg+2)=first_shell(bb,3)
              END IF
           ELSE
              microstate(gg+2)=first_shell(bb,3)
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO

     DO ii=4,5
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg+1:gg+2)=(/ first_shell(bb,3:4) /)
           IF ((first_shell(bb,1)==kk).or.(first_shell(bb,2)==kk)) THEN
              IF (first_shell(bb,1)==kk) THEN
                 microstate(gg)=first_shell(bb,2)
              ELSE
                 microstate(gg)=first_shell(bb,1)
              END IF
              IF ((first_shell(bb,1)==kk).and.(first_shell(bb,2)==kk)) THEN
                 ! print*,'pproblema ssn2' !!! Comprobar esto
                 ! print*,a,'|',first_shell(b,:)
                 microstate(gg)=first_shell(bb,1)
                 ! print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
                 ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)
              END IF
           ELSE
              microstate(gg)=first_shell(bb,1) !!!! esto esta mal, hay que elegir entre el 1 y el 2
              !print*,'siiiiiiiiiiiiiiiiii'
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO
     
     !print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
     ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)

     !Quito indices de mols

     !!filtro=.true.
     !!aux=microstate
     !!mss_ind_wat=microstate
     !!DO i=1,17
     !!   IF (aux(i)<=0) filtro(i)=.false.
     !!END DO
     !! 
     !!shell_w(a,:)=0
     !!g=0
     !! 
     !!DO j=1,17
     !!   IF (filtro(j)==.true.) THEN
     !!      b=aux(j)
     !!      g=g+1
     !!      shell_w(a,g)=b
     !!      DO i=1,17
     !!         IF (filtro(i)==.true.) THEN
     !!            IF (aux(i)==b) THEN
     !!               microstate(i)=j
     !!               filtro(i)=.false.
     !!            END IF
     !!         END IF
     !!      END DO
     !!   END IF
     !!END DO

     mss(kk,:)=microstate

  END DO


  DEALLOCATE(dist_first_shell,first_shell,bonds_o)
  DEALLOCATE(microstate)


END SUBROUTINE ind_wat_limit_4_nosim

SUBROUTINE ind_wat_limit_4_nosim_prot (mss,aux,filt_water,hbs,hbdists,num_wats,num_atoms,num_hbs)


  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_atoms,num_wats,num_hbs
  INTEGER,DIMENSION(num_atoms,3),INTENT(IN)::aux
  LOGICAL,DIMENSION(num_atoms),INTENT(IN)::filt_water
  INTEGER,DIMENSION(num_hbs,3),INTENT(IN)::hbs
  DOUBLE PRECISION,DIMENSION(num_hbs),INTENT(IN)::hbdists
  INTEGER,DIMENSION(num_wats,17),INTENT(OUT)::mss
  
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_first_shell
  INTEGER,DIMENSION(:,:),ALLOCATABLE::first_shell
  INTEGER,DIMENSION(:),ALLOCATABLE::bonds_o
  INTEGER,DIMENSION(:),ALLOCATABLE::microstate
  
  INTEGER::ii,jj,kk,bb,gg
  INTEGER::ind_oh,hi,ind_o
  
  LOGICAL::interr_oh,interr_o

  mss=0
  ALLOCATE(dist_first_shell(num_wats,4),first_shell(num_wats,4),bonds_o(num_wats))
  ALLOCATE(microstate(17))

  first_shell=0
  dist_first_shell=0.0d0
  bonds_o=0


  IF (ANY((/2,3/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        interr_oh=.FALSE.
        interr_o=.FALSE.
        IF (filt_water(hbs(ii,1)+1)) THEN
           interr_oh=.TRUE.
           ind_oh=aux(hbs(ii,1)+1,1)+1
           hi=aux(hbs(ii,2)+1,2)
        ELSE
           ind_oh=-(hbs(ii,2)+1)
        END IF
        IF (filt_water(hbs(ii,3)+1)) THEN
           interr_o=.TRUE.
           ind_o=aux(hbs(ii,3)+1,1)+1
        ELSE
           ind_o=-(hbs(ii,3)+1)
        END IF
        IF (interr_oh==.TRUE.) THEN
           IF (first_shell(ind_oh,hi)==0) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           ELSE
              IF (dist_first_shell(ind_oh,hi)>hbdists(ii)) THEN
                 first_shell(ind_oh,hi)=ind_o
                 dist_first_shell(ind_oh,hi)=hbdists(ii)
              END IF
           END IF
        END IF
        IF (interr_o==.TRUE.) THEN
           IF (bonds_o(ind_o)==0) THEN
              first_shell(ind_o,3)=ind_oh
              dist_first_shell(ind_o,3)=hbdists(ii)
              bonds_o(ind_o)=1
           ELSE
              IF (bonds_o(ind_o)==1) THEN
                 IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                    first_shell(ind_o,4)=ind_oh
                    dist_first_shell(ind_o,4)=hbdists(ii)
                 ELSE
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 END IF
                 bonds_o(ind_o)=2
              ELSE
                 IF (bonds_o(ind_o)==2) THEN
                    IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                       first_shell(ind_o,4)=first_shell(ind_o,3)
                       dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                       first_shell(ind_o,3)=ind_oh
                       dist_first_shell(ind_o,3)=hbdists(ii)
                    ELSE
                       IF (dist_first_shell(ind_o,4)>hbdists(ii)) THEN
                          first_shell(ind_o,4)=ind_oh
                          dist_first_shell(ind_o,4)=hbdists(ii)
                       END IF
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  ELSE IF (ANY((/1/)==definition_hbs)) THEN
     DO ii=1,num_hbs
        interr_oh=.FALSE.
        interr_o=.FALSE.
        IF (filt_water(hbs(ii,1)+1)) THEN
           interr_oh=.TRUE.
           ind_oh=aux(hbs(ii,1)+1,1)+1
           hi=aux(hbs(ii,2)+1,2)
        ELSE
           ind_oh=-(hbs(ii,2)+1)
        END IF
        IF (filt_water(hbs(ii,3)+1)) THEN
           interr_o=.TRUE.
           ind_o=aux(hbs(ii,3)+1,1)+1
        ELSE
           ind_o=-(hbs(ii,3)+1)
        END IF
        IF (interr_oh==.TRUE.) THEN
           IF (first_shell(ind_oh,hi)==0) THEN
              first_shell(ind_oh,hi)=ind_o
              dist_first_shell(ind_oh,hi)=hbdists(ii)
           ELSE
              IF (dist_first_shell(ind_oh,hi)<hbdists(ii)) THEN
                 first_shell(ind_oh,hi)=ind_o
                 dist_first_shell(ind_oh,hi)=hbdists(ii)
              END IF
           END IF
        END IF
        IF (interr_o==.TRUE.) THEN
           IF (bonds_o(ind_o)==0) THEN
              first_shell(ind_o,3)=ind_oh
              dist_first_shell(ind_o,3)=hbdists(ii)
              bonds_o(ind_o)=1
           ELSE
              IF (bonds_o(ind_o)==1) THEN
                 IF (dist_first_shell(ind_o,3)>hbdists(ii)) THEN
                    first_shell(ind_o,4)=ind_oh
                    dist_first_shell(ind_o,4)=hbdists(ii)
                 ELSE
                    first_shell(ind_o,4)=first_shell(ind_o,3)
                    dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                    first_shell(ind_o,3)=ind_oh
                    dist_first_shell(ind_o,3)=hbdists(ii)
                 END IF
                 bonds_o(ind_o)=2
              ELSE
                 IF (bonds_o(ind_o)==2) THEN
                    IF (dist_first_shell(ind_o,3)<hbdists(ii)) THEN
                       first_shell(ind_o,4)=first_shell(ind_o,3)
                       dist_first_shell(ind_o,4)=dist_first_shell(ind_o,3)
                       first_shell(ind_o,3)=ind_oh
                       dist_first_shell(ind_o,3)=hbdists(ii)
                    ELSE
                       IF (dist_first_shell(ind_o,4)<hbdists(ii)) THEN
                          first_shell(ind_o,4)=ind_oh
                          dist_first_shell(ind_o,4)=hbdists(ii)
                       END IF
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
  END IF

  DO kk=1,num_wats

     microstate=0
     microstate(1)=kk
     
     microstate(2:5)=(/ first_shell(kk,:) /)

     gg=6
     DO ii=2,3
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg:gg+1)=(/ first_shell(bb,1:2) /)
           IF ((first_shell(bb,3)==kk).or.(first_shell(bb,4)==kk)) THEN
              IF (first_shell(bb,3)==kk) THEN
                 microstate(gg+2)=first_shell(bb,4)
              ELSE
                 microstate(gg+2)=first_shell(bb,3)
              END IF
              IF ((first_shell(bb,3)==kk).and.(first_shell(bb,4)==kk)) THEN
                 !print*,'problema ssn1' !!! Comprobar esto
                 microstate(gg+2)=first_shell(bb,3)
              END IF
           ELSE
              microstate(gg+2)=first_shell(bb,3)
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO

     DO ii=4,5
        bb=microstate(ii)
        IF (bb>0) THEN
           microstate(gg+1:gg+2)=(/ first_shell(bb,3:4) /)
           IF ((first_shell(bb,1)==kk).or.(first_shell(bb,2)==kk)) THEN
              IF (first_shell(bb,1)==kk) THEN
                 microstate(gg)=first_shell(bb,2)
              ELSE
                 microstate(gg)=first_shell(bb,1)
              END IF
              IF ((first_shell(bb,1)==kk).and.(first_shell(bb,2)==kk)) THEN
                 ! print*,'pproblema ssn2' !!! Comprobar esto
                 ! print*,a,'|',first_shell(b,:)
                 microstate(gg)=first_shell(bb,1)
                 ! print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
                 ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)
              END IF
           ELSE
              microstate(gg)=first_shell(bb,1) !!!! esto esta mal, hay que elegir entre el 1 y el 2
              !print*,'siiiiiiiiiiiiiiiiii'
           END IF
           gg=gg+3
        ELSE
           microstate(gg:gg+2)=0
           gg=gg+3
        END IF
     END DO
     
     !print'(I5,A,4I5,A,3I5,A,3I5,A,3I5,A,3I5)',microstate(1),'||',microstate(2:5),'||',microstate(6:8),&
     ! & '||',microstate(9:11),'||',microstate(12:14),'||',microstate(15:17)

     !Quito indices de mols

     !!filtro=.true.
     !!aux=microstate
     !!mss_ind_wat=microstate
     !!DO i=1,17
     !!   IF (aux(i)<=0) filtro(i)=.false.
     !!END DO
     !! 
     !!shell_w(a,:)=0
     !!g=0
     !! 
     !!DO j=1,17
     !!   IF (filtro(j)==.true.) THEN
     !!      b=aux(j)
     !!      g=g+1
     !!      shell_w(a,g)=b
     !!      DO i=1,17
     !!         IF (filtro(i)==.true.) THEN
     !!            IF (aux(i)==b) THEN
     !!               microstate(i)=j
     !!               filtro(i)=.false.
     !!            END IF
     !!         END IF
     !!      END DO
     !!   END IF
     !!END DO

     mss(kk,:)=microstate

  END DO


  DEALLOCATE(dist_first_shell,first_shell,bonds_o)
  DEALLOCATE(microstate)


END SUBROUTINE ind_wat_limit_4_nosim_prot


SUBROUTINE remove_index_mol(mss_ind,num_mss,mss)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::num_mss
  INTEGER,DIMENSION(num_mss,17),INTENT(IN)::mss_ind
  INTEGER,DIMENSION(num_mss,17),INTENT(OUT)::mss

  INTEGER::ii,jj,gg,ll,bb
  INTEGER,DIMENSION(17)::aux_mss,aux_mss2
  LOGICAL,DIMENSION(17)::filtro


  DO jj=1,num_mss
     filtro=.true.
     aux_mss2=mss_ind(jj,:)
     aux_mss=aux_mss2

     DO ii=1,17
        IF (aux_mss(ii)<=0) filtro(ii)=.false.
     END DO
  
     gg=0
  
     DO ii=1,17
        IF (filtro(ii)==.true.) THEN
           bb=aux_mss(ii)
           gg=gg+1
           DO ll=1,17
              IF (filtro(ll)==.true.) THEN
                 IF (aux_mss(ll)==bb) THEN
                    aux_mss2(ll)=ii
                    filtro(ll)=.false.
                 END IF
              END IF
           END DO
        END IF
     END DO
  
     mss(jj,:)=aux_mss2(:)

  END DO


END SUBROUTINE remove_index_mol

SUBROUTINE remove_permutations_limit_4_nosim(mss,mss_ind,num_mss)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::num_mss
  INTEGER,DIMENSION(num_mss,17),INTENT(INOUT)::mss,mss_ind

  INTEGER::i,j,g,h,ii,kk,vigilo
  INTEGER,DIMENSION(17)::microstate,ms_short2,bb,key,key_aux
  LOGICAL,DIMENSION(17)::filtro
  INTEGER,DIMENSION(4,2)::aux_permut
  LOGICAL::interruptor
  INTEGER::ceros4,ceros5,ceros,x_ellos
  INTEGER::x_primera4,x_primera5,x_primera
  INTEGER::x_core,x_core4,x_core5
  INTEGER::x_segunda4,x_segunda5,x_segunda

  DO kk=1,num_mss

     key=mss_ind(kk,:)
     key_aux=key
     ms_short2=mss(kk,:)
     microstate=ms_short2
     aux_permut(1,1)=6
     aux_permut(1,2)=7
     aux_permut(2,1)=9
     aux_permut(2,2)=10
     aux_permut(3,1)=13
     aux_permut(3,2)=14
     aux_permut(4,1)=16
     aux_permut(4,2)=17
     
     DO h=1,2
        i=aux_permut(h,1)
        j=aux_permut(h,2)
        IF (ms_short2(i)>0) THEN
           IF (ms_short2(i)>ms_short2(j)) THEN
              microstate(i)=ms_short2(j)
              microstate(j)=ms_short2(i)
              key_aux(i)=key(j)
              key_aux(j)=key(i)
              filtro=.true.
              DO g=1,17
                 IF ((microstate(g)==i).and.(filtro(g)==.true.)) THEN
                    microstate(g)=j
                    filtro(g)=.false.
                 END IF
              END DO
              DO g=1,17
                 IF ((microstate(g)==j).and.(filtro(g)==.true.)) THEN
                    microstate(g)=i
                    filtro(g)=.false.
                 END IF
              END DO
              ms_short2=microstate
              key=key_aux
           END IF
        END IF
     END DO
     key=key_aux
     ms_short2=microstate
     
     
     
     ceros4=0
     ceros5=0
     ceros=0
     x_core4=0
     x_core5=0
     x_core=0
     x_primera4=0
     x_primera5=0
     x_primera=0
     x_segunda4=0
     x_segunda5=0
     x_segunda=0
     x_ellos=0
     
     DO i=12,14
        IF (ms_short2(i)==0) THEN
           ceros4=ceros4+1
        END IF
        IF (ms_short2(i)==1) THEN
           x_core4=x_core4+1
        END IF
        IF ((ms_short2(i)<=5).and.(ms_short2(i)>=2)) THEN
           x_primera4=x_primera4+1
        END IF
        IF ((ms_short2(i)<=11).and.(ms_short2(i)>=6)) THEN
           x_segunda4=x_segunda4+1
        END IF
        IF ((ms_short2(i)/=i).and.(ms_short2(i)>11)) THEN
           x_ellos=x_ellos+1
        END IF
     END DO
     DO i=15,17
        IF (ms_short2(i)==0) THEN
           ceros5=ceros5+1
        ELSE
           IF (ms_short2(i)==1) THEN
              x_core5=x_core5+1
           ELSE
              IF ((ms_short2(i)<=5).and.(ms_short2(i)>=2)) THEN
                 x_primera5=x_primera5+1
              ELSE
                 IF ((ms_short2(i)<=11).and.(ms_short2(i)>=6)) THEN
                    x_segunda5=x_segunda5+1
                 ELSE
                    IF ((ms_short2(i)/=i).and.(ms_short2(i)>11)) THEN
                       x_ellos=x_ellos+1
                    END IF
                 END IF
              END IF
           END IF
        END IF
     END DO
     ceros=ceros4+ceros5
     x_core=x_core4+x_core5
     x_primera=x_primera4+x_primera5
     x_segunda=x_segunda4+x_segunda5
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     interruptor=.false.
     
     IF (ms_short2(5)==0) THEN
        interruptor=.true.
     END IF
     IF (interruptor==.false.) THEN
        IF (x_core/=0) THEN
           IF (x_core==1) THEN
              IF (x_core5==1) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 microstate=ms_short2
                 interruptor=.true.
              END IF
           ELSE
              !print*,'problema x_core',frame,mol
              !print 117, ms_short2(:)
              !STOP
           END IF
        END IF
     END IF
     
     IF (interruptor==.false.) THEN
        IF (ceros4/=ceros5) interruptor=.true.
        IF (ceros5<ceros4) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
        END IF
     END IF
     IF (interruptor==.false.) THEN
        IF (x_primera4/=x_primera5) interruptor=.true.
        IF (x_primera5<x_primera4) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
        END IF
     END IF
     IF (interruptor==.false.) THEN
        IF (x_segunda4/=x_segunda5) interruptor=.true.
        IF (x_segunda5<x_segunda4) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
        END IF
     END IF
     
     !!
     microstate=ms_short2
     key_aux=key
     DO h=3,4
        i=aux_permut(h,1)
        j=aux_permut(h,2)
        IF (ms_short2(i)/=0) THEN
           IF (ms_short2(i)>ms_short2(j)) THEN
              key_aux(i)=key(j)
              key_aux(j)=key(i)
              microstate(i)=ms_short2(j)
              microstate(j)=ms_short2(i)
              filtro=.true.
              DO g=1,17
                 IF ((microstate(g)==i).and.(filtro(g)==.true.)) THEN
                    microstate(g)=j
                    filtro(g)=.false.
                 END IF
              END DO
              DO g=1,17
                 IF ((microstate(g)==j).and.(filtro(g)==.true.)) THEN
                    microstate(g)=i
                    filtro(g)=.false.
                 END IF
              END DO
              ms_short2=microstate
              key=key_aux
           END IF
        END IF
     END DO
     ms_short2=microstate
     key=key_aux
     !!
     
     
     !sigo
     
     IF (interruptor==.false.) THEN
        IF (ceros>0) THEN
           IF ((ms_short2(12)==0).and.(ms_short2(15)/=0)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(12)/=0).and.(ms_short2(15)==0)) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
           IF (interruptor==.false.) THEN
              IF ((ms_short2(13)==0).and.(ms_short2(16)/=0)) THEN
                 interruptor=.true.
              ELSE
                 IF ((ms_short2(13)/=0).and.(ms_short2(16)==0)) THEN
                    CALL DOY_VUELTA(ms_short2)
                    CALL DOY_VUELTA_KEY (key,key_aux)
                    interruptor=.true.
                 END IF
              END IF
           END IF
           IF (interruptor==.false.) THEN
              IF ((ms_short2(14)==0).and.(ms_short2(17)/=0)) THEN
                 interruptor=.true.
              ELSE
                 IF ((ms_short2(14)/=0).and.(ms_short2(17)==0)) THEN
                    CALL DOY_VUELTA(ms_short2)
                    CALL DOY_VUELTA_KEY (key,key_aux)
                    interruptor=.true.
                 END IF
              END IF
           END IF
        END IF
     END IF
     
     
     IF (interruptor==.false.) THEN
        IF ((ms_short2(12)==12).and.(ms_short2(15)/=15)) THEN
           interruptor=.true.
        ELSE
           IF ((ms_short2(12)/=12).and.(ms_short2(15)==15)) THEN
              CALL DOY_VUELTA(ms_short2)
              CALL DOY_VUELTA_KEY (key,key_aux)
              interruptor=.true.
           END IF
        END IF
        IF (interruptor==.false.) THEN
           IF ((ms_short2(13)==13).and.(ms_short2(16)/=16)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(13)/=13).and.(ms_short2(16)==16)) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
        END IF
        IF (interruptor==.false.) THEN
           IF ((ms_short2(14)==14).and.(ms_short2(17)/=17)) THEN
              interruptor=.true.
           ELSE
              IF ((ms_short2(14)/=14).and.(ms_short2(17)==17)) THEN
                 CALL DOY_VUELTA(ms_short2)
                 CALL DOY_VUELTA_KEY (key,key_aux)
                 interruptor=.true.
              END IF
           END IF
        END IF
     END IF
     microstate=ms_short2
     key_aux=key
     
     !bb=0
     !! Cruces entre ellos:
     ! - Un cruce
     ! bb(12:17)=(/ 12,13,14,14,16,17 /) !SII Elijo este como representante
     ! bb(12:17)=(/ 12,13,14,13,16,17 /) !SII ---> 12,13,14,14,16,17
     ! bb(12:17)=(/ 12,13,14,15,16,12 /) !NO
     !bb(12:17)=(/ 12,13,14,15,12,17 /) !SII
     ! - Dos cruces
     ! bb(12:17)=(/ 12,13,14,14,16,12 /) !NO
     ! bb(12:17)=(/ 12,13,14,14,12,17 /) !SII Elijo este como representante
     ! bb(12:17)=(/ 12,13,14,13,16,12 /) !NO
     ! bb(12:17)=(/ 12,13,14,13,12,17 /) !SII ---> 12,13,14,14,12,17
     ! - Tres cruces da igual
     !! Cruces entre ellos con la primera capa tambien:
     ! - Un cruce entre ellos y uno con primera capa
     ! bb(12:17)=(/ 5,13,14,14,4,17 /)
     ! bb(12:17)=(/ 5,13,14,14,16,4 /)
     !
     ! Corrijo esto:
     
     IF (ms_short2(16)==12) THEN
        CALL DOY_VUELTA(ms_short2)
        CALL DOY_VUELTA_KEY (key,key_aux)
     END IF
     microstate=ms_short2
     key_aux=key
     IF (ms_short2(15)==13) THEN
        key(13)=key_aux(14)
        key(14)=key_aux(13)
        key_aux(13:14)=key(13:14)
        ms_short2(13)=microstate(14)
        ms_short2(14)=microstate(13)
        microstate(13:14)=ms_short2(13:14)   !!! Esto lo acabo de cambiar del de Roman
        filtro=.true.
        DO i=12,17
           j=ms_short2(i)
           IF ((j>11).and.(filtro(i)==.true.)) THEN
              ms_short2(i)=i
              filtro(i)=.false.
              DO ii=i+1,17
                 IF ((microstate(ii)==j).and.(filtro(ii)==.true.)) THEN
                    ms_short2(ii)=i
                    filtro(ii)=.false.
                 END IF
              END DO
           END IF
        END DO
     END IF
     microstate=ms_short2
     key_aux=key
     
     
     IF (interruptor==.false.) THEN
        IF ((ms_short2(13)==3).and.(ms_short2(16)==2)) THEN
           CALL DOY_VUELTA(ms_short2)
           CALL DOY_VUELTA_KEY (key,key_aux)
           microstate=ms_short2
           key_aux=key
        END IF
     END IF
     
     
     mss_ind(kk,:)=key(:)
     mss(kk,:)=microstate(:)

  END DO

!!!!!!!
     
     !IF (0) THEN
     ! interruptor=.true.
     ! DO i=12,17
     ! IF (ms_short2(i)/=bb(i)) THEN
     ! interruptor=.false.
     ! EXIT
     ! END IF
     ! END DO
     ! IF (interruptor==.true.) THEN
     ! print*,'...>',frame,mol
     ! print 117,ms_short2(:)
     ! END IF
     !END IF
     
     
     !IF (interruptor==.false.) THEN
     ! print*,'TATE',frame,mol
     ! print 117,ms_short2
     !END IF
     
     
     !117 format (1I3," |",4I3," |",3I3," |",3I3," |",3I3," |",3I3)


END SUBROUTINE remove_permutations_limit_4_nosim


SUBROUTINE DOY_VUELTA (ms_short2)

  IMPLICIT NONE

  INTEGER,DIMENSION(17),INTENT(INOUT)::ms_short2
  INTEGER::i,j,ii
  INTEGER,DIMENSION(17)::microstate
  LOGICAL,DIMENSION(17)::filtro
  
  filtro=.true.
  microstate=ms_short2
  
  ms_short2(4)=microstate(5)
  ms_short2(5)=microstate(4)
  ms_short2(12:14)=microstate(15:17)
  ms_short2(15:17)=microstate(12:14)
  microstate=ms_short2
  
  DO i=1,17
     j=ms_short2(i)
     IF ((j>1).and.(filtro(i)==.true.)) THEN
        ms_short2(i)=i
        filtro(i)=.false.
        DO ii=i+1,17
           IF ((microstate(ii)==j).and.(filtro(ii)==.true.)) THEN
              ms_short2(ii)=i
              filtro(ii)=.false.
           END IF
        END DO
     END IF
  END DO

END SUBROUTINE DOY_VUELTA

SUBROUTINE DOY_VUELTA_KEY (key,key_aux)
  
  IMPLICIT NONE
  
  INTEGER,DIMENSION(17),INTENT(INOUT)::key,key_aux

  key_aux=key
  key(4)=key_aux(5)
  key(5)=key_aux(4)
  key(12:14)=key_aux(15:17)
  key(15:17)=key_aux(12:14)
  key_aux=key

END SUBROUTINE DOY_VUELTA_KEY


SUBROUTINE addbonds (tipo,mss,mss_ind,bonds,num_wats,num_bonds)


  IMPLICIT NONE

  INTEGER,INTENT(IN)::tipo,num_wats,num_bonds
  INTEGER,DIMENSION(num_wats,17),INTENT(INOUT)::mss
  INTEGER,DIMENSION(num_wats,17),INTENT(INOUT)::mss_ind
  INTEGER,DIMENSION(num_bonds),INTENT(IN)::bonds
  
  INTEGER::ii,jj,kk

  LOGICAL,DIMENSION(num_wats)::dentro
  LOGICAL,DIMENSION(17)::corrijo

  dentro=.FALSE.
  DO ii=1,num_bonds
     dentro(bonds(ii)+1)=.TRUE.
  END DO

  SELECT CASE (tipo)
     CASE (1)
        DO jj=1,num_wats
           corrijo=.FALSE.
           DO kk=1,5
              IF (dentro(mss_ind(jj,kk))==.TRUE.) THEN
                 IF (kk==1) THEN
                    IF (mss_ind(jj,4)==0) THEN
                       corrijo(4)=.TRUE.
                    ELSE
                       IF (mss_ind(jj,5)==0) THEN
                          corrijo(5)=.TRUE.
                       END IF
                    END IF
                 ELSEIF (kk==2) THEN
                    IF (mss_ind(jj,8)==0) THEN
                       corrijo(8)=.TRUE.
                    END IF
                 ELSEIF (kk==3) THEN
                    IF (mss_ind(jj,11)==0) THEN
                       corrijo(11)=.TRUE.
                    END IF
                 ELSEIF (kk==4) THEN
                    IF (mss_ind(jj,13)==0) THEN
                       corrijo(13)=.TRUE.
                    ELSE
                       IF (mss_ind(jj,14)==0) THEN
                          corrijo(14)=.TRUE.
                       END IF
                    END IF
                 ELSEIF (kk==5) THEN
                    IF (mss_ind(jj,16)==0) THEN
                       corrijo(16)=.TRUE.
                    ELSE
                       IF (mss_ind(jj,17)==0) THEN
                          corrijo(17)=.TRUE.
                       END IF
                    END IF
                 END IF
              END IF
           END DO
           DO kk=1,17
              IF (corrijo(kk)==.TRUE.) THEN
                 mss_ind(jj,kk)=-1
                 mss(jj,kk)=-1
              END IF
           END DO
        END DO

     CASE (2)

        DO jj=1,num_wats
           DO kk=1,17
              ii=mss_ind(jj,kk)
              IF (ii>0) THEN
                 IF (dentro(ii)==.TRUE.) THEN
                    mss(jj,kk)=-mss(jj,kk)
                 END IF
              END IF
           END DO
        END DO


     CASE (3)

        DO jj=1,num_wats
           DO kk=1,5
              ii=mss_ind(jj,kk)
              IF (ii>0) THEN
                 IF (dentro(ii)==.TRUE.) THEN
                    mss(jj,kk)=-mss(jj,kk)
                 END IF
              END IF
           END DO
        END DO



     CASE DEFAULT

        PRINT*,'NO OPTION CHOOSEN'

     END SELECT

END SUBROUTINE addbonds

END MODULE GLOB
