MODULE HBONDS

  !! Indices:

  INTEGER::natoms1,natoms2
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::coor1,coor2
  DOUBLE PRECISION,DIMENSION(3,3)::Lbox1,Lbox2,Lbox_2,Lbox

  
  !! 
  INTEGER::maxHBs
  INTEGER::Ldon1,Ldon2,Lacc1,Lacc2,maxH1,maxH2
  INTEGER,DIMENSION(:),ALLOCATABLE::ind_don1,ind_don2,ind_acc1,ind_acc2,num_Hdon1,num_Hdon2
  INTEGER,DIMENSION(:,:),ALLOCATABLE::ind_Hdon1,ind_Hdon2

  !! Output:
  INTEGER::num_hbs1,num_hbs2,num_hbs3,num_hbs4
  INTEGER,DIMENSION(:,:),ALLOCATABLE::salida1,salida2,salida3,salida4


  PARAMETER (maxHBs=6)

CONTAINS

  SUBROUTINE INITIALIZE_COORS_MEMORY1()

    ALLOCATE (coor1(natoms1,3))

  END SUBROUTINE INITIALIZE_COORS_MEMORY1

  SUBROUTINE INITIALIZE_COORS_MEMORY2()

    ALLOCATE (coor2(natoms2,3))

  END SUBROUTINE INITIALIZE_COORS_MEMORY2


  SUBROUTINE INITIALIZE_MEMORY_IN()

    IF (Ldon1>0) ALLOCATE(ind_don1(Ldon1),num_Hdon1(Ldon1),ind_Hdon1(Ldon1,maxH1))
    IF (Ldon2>0) ALLOCATE(ind_don2(Ldon2),num_Hdon2(Ldon2),ind_Hdon2(Ldon2,maxH2))
    IF (Lacc1>0) ALLOCATE(ind_acc1(Lacc1))
    IF (Lacc2>0) ALLOCATE(ind_acc2(Lacc2))

  END SUBROUTINE INITIALIZE_MEMORY_IN

  SUBROUTINE FREE_MEMORY_IN()
    
    IF (ALLOCATED(ind_don1)) DEALLOCATE(ind_don1)
    IF (ALLOCATED(num_Hdon1)) DEALLOCATE(num_Hdon1)
    IF (ALLOCATED(ind_Hdon1)) DEALLOCATE(ind_Hdon1)
    IF (ALLOCATED(ind_acc1)) DEALLOCATE(ind_acc1)

    IF (ALLOCATED(ind_don2)) DEALLOCATE(ind_don2)
    IF (ALLOCATED(num_Hdon2)) DEALLOCATE(num_Hdon2)
    IF (ALLOCATED(ind_Hdon2)) DEALLOCATE(ind_Hdon2)
    IF (ALLOCATED(ind_acc2)) DEALLOCATE(ind_acc2)

  END SUBROUTINE FREE_MEMORY_IN

  SUBROUTINE FREE_MEMORY_COOR()
    
    IF (ALLOCATED(coor1)) DEALLOCATE(coor1)
    IF (ALLOCATED(coor2)) DEALLOCATE(coor2)

  END SUBROUTINE FREE_MEMORY_COOR

  SUBROUTINE FREE_MEMORY_OUT()
    
    IF (ALLOCATED(salida1)) DEALLOCATE(salida1)
    IF (ALLOCATED(salida2)) DEALLOCATE(salida2)
    IF (ALLOCATED(salida3)) DEALLOCATE(salida3)
    IF (ALLOCATED(salida4)) DEALLOCATE(salida4)

  END SUBROUTINE FREE_MEMORY_OUT


  SUBROUTINE DIFF_SET (r_param,ang_param)

    DOUBLE PRECISION,INTENT(IN)::r_param,ang_param
    INTEGER::i,j,k,l,ii,jj,g,gg,kk

    INTEGER,DIMENSION(Ldon1,maxH1)::num_HB_don1
    INTEGER,DIMENSION(Ldon2,maxH2)::num_HB_don2
    INTEGER,DIMENSION(:,:,:),ALLOCATABLE::HB_don1,HB_don2
    DOUBLE PRECISION,DIMENSION(maxH1,3)::vect_Hdon1
    DOUBLE PRECISION,DIMENSION(maxH2,3)::vect_Hdon2
    DOUBLE PRECISION,DIMENSION(3)::vect_aux,pos_d


    ALLOCATE(HB_don1(Ldon1,maxH1,maxHBs))
    ALLOCATE(HB_don2(Ldon2,maxH2,maxHBs))

    num_hbs1=0
    num_hbs2=0

    num_HB_don1=0
    num_HB_don2=0

    aux_cos=COSD(ang_param)

    Lbox=Lbox1
    Lbox_2=Lbox/2.0d0

    DO i=1,Ldon1
       ii=ind_don2(i)
       pos_d=coor1(ii,:)
       DO k=1,num_Hdon1(i)
          kk=ind_Hdon1(i,k)
          vect_aux=coor1(kk,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon1(k,:)=vect_aux/norm
       END DO

       DO j=1,Lacc2
          jj=ind_acc2(j)
          vect_aux=coor2(jj,:)-pos_d
          CALL PBC(vect_aux)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          IF (norm<r_param) THEN
             vect_aux=vect_aux/norm
             DO k=1,num_Hdon1(i)
                cos_ang=dot_product(vect_aux,vect_Hdon1(k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs1=num_hbs1+1
                   g=num_HB_don1(i,k)
                   gg=g+1
                   HB_don1(i,k,gg)=jj-1
                   num_HB_don1(i,k)=gg
                END IF
             END DO
          END IF
       END DO

    END DO

    DO i=1,Ldon2
       ii=ind_don2(i)
       pos_d=coor2(ii,:)
       DO k=1,num_Hdon2(i)
          kk=ind_Hdon2(i,k)
          vect_aux=coor2(kk,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon2(k,:)=vect_aux/norm
       END DO

       DO j=1,Lacc1
          jj=ind_acc1(j)
          vect_aux=coor1(jj,:)-pos_d
          CALL PBC(vect_aux)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          IF (norm<r_param) THEN
             vect_aux=vect_aux/norm
             DO k=1,num_Hdon2(i)
                cos_ang=dot_product(vect_aux,vect_Hdon2(k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs2=num_hbs2+1
                   g=num_HB_don2(i,k)
                   gg=g+1
                   HB_don2(i,k,gg)=jj-1
                   num_HB_don2(i,k)=gg
                END IF
             END DO
          END IF
       END DO

    END DO

    ALLOCATE(salida1(num_hbs1,3),salida2(num_hbs2,3))

    num_hbs1=0

    DO i=1,Ldon1
       DO j=1,num_Hdon1(i)
          DO k=1,num_HB_don1(i,j)
             num_hbs1=num_hbs1+1
             salida1(num_hbs1,1)=ind_don1(i)
             salida1(num_hbs1,2)=ind_Hdon1(i,j)
             salida1(num_hbs1,3)=HB_don1(i,j,k)
          END DO
       END DO
    END DO

    num_hbs2=0

    DO i=1,Ldon2
       DO j=1,num_Hdon2(i)
          DO k=1,num_HB_don2(i,j)
             num_hbs2=num_hbs2+1
             salida2(num_hbs2,1)=ind_don2(i)
             salida2(num_hbs2,2)=ind_Hdon2(i,j)
             salida2(num_hbs2,3)=HB_don2(i,j,k)
          END DO
       END DO
    END DO

    DEALLOCATE(HB_don1,HB_don2)

  END SUBROUTINE DIFF_SET

  SUBROUTINE SAME_SET (r_param,ang_param)

    DOUBLE PRECISION,INTENT(IN)::r_param,ang_param
    INTEGER::i,j,k,l,ii,jj,g,gg,kk


    INTEGER,DIMENSION(Ldon1,maxH1)::num_HB_don1
    INTEGER,DIMENSION(Ldon2,maxH2)::num_HB_don2
    INTEGER,DIMENSION(:,:,:),ALLOCATABLE::HB_don1,HB_don2
    DOUBLE PRECISION,DIMENSION(Ldon1,maxH1,3)::vect_Hdon1
    DOUBLE PRECISION,DIMENSION(Ldon2,maxH2,3)::vect_Hdon2
    DOUBLE PRECISION,DIMENSION(3)::vect_aux,pos_d

    num_hbs1=0
    num_hbs2=0
    num_hbs3=0
    num_hbs4=0
    num_HB_don1=0
    num_HB_don1=0

    aux_cos=COSD(ang_param)
    
    Lbox=Lbox1
    Lbox_2=Lbox/2.0d0

    DO i=1,Ldon1
       ii=ind_don1(i)
       pos_d=coor1(ii,:)
       DO k=1,num_Hdon1(i)
          kk=ind_Hdon1(i,k)
          vect_aux=coor1(kk,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon1(i,k,:)=vect_aux/norm
       END DO
    END DO
    DO i=1,Ldon2
       ii=ind_don2(i)
       pos_d=coor1(ii,:)
       DO k=1,num_Hdon2(i)
          kk=ind_Hdon2(i,k)
          vect_aux=coor1(kk,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon2(i,k,:)=vect_aux/norm
       END DO
    END DO


    ALLOCATE(HB_don1(Ldon1,maxH1,maxHBs))
    ALLOCATE(HB_don2(Ldon2,maxH2,maxHBs))


    num_HB_don1=0
    num_HB_don2=0

    !! Same Same

    DO i=1,Ldon1
       ii=ind_don1(i)
       pos_d=coor1(ii,:)
       DO j=i+1,Ldon1
          jj=ind_don1(j)
          vect_aux=coor1(jj,:)-pos_d
          CALL PBC(vect_aux)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          IF (norm<r_param) THEN
             vect_aux=vect_aux/norm
             DO k=1,num_Hdon1(i)
                cos_ang=dot_product(vect_aux,vect_Hdon1(i,k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs1=num_hbs1+1
                   g=num_HB_don1(i,k)
                   gg=g+1
                   HB_don1(i,k,gg)=jj-1
                   num_HB_don1(i,k)=gg
                END IF
             END DO
             DO k=1,num_Hdon1(j)
                cos_ang=dot_product(-vect_aux,vect_Hdon1(j,k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs1=num_hbs1+1
                   g=num_HB_don1(j,k)
                   gg=g+1
                   HB_don1(j,k,gg)=ii-1
                   num_HB_don1(j,k)=gg
                END IF
             END DO
          END IF
       END DO
    END DO

    ALLOCATE(salida1(num_hbs1,3))

    num_hbs1=0

    DO i=1,Ldon1
       DO j=1,num_Hdon1(i)
          DO k=1,num_HB_don1(i,j)
             num_hbs1=num_hbs1+1
             salida1(num_hbs1,1)=ind_don1(i)
             salida1(num_hbs1,2)=ind_Hdon1(i,j)
             salida1(num_hbs1,3)=HB_don1(i,j,k)
          END DO
       END DO
    END DO

    !! Same - acceptor2

    DEALLOCATE(HB_don1)
    ALLOCATE(HB_don1(Ldon1,maxH1,maxHBs))

    num_HB_don1=0

    DO i=1,Ldon1
       ii=ind_don1(i)
       pos_d=coor1(ii,:)    
       DO j=1,Lacc2
          jj=ind_acc2(j)
          vect_aux=coor1(jj,:)-pos_d
          CALL PBC(vect_aux)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          IF (norm<r_param) THEN
             vect_aux=vect_aux/norm
             DO k=1,num_Hdon1(i)
                cos_ang=dot_product(vect_aux,vect_Hdon1(i,k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs2=num_hbs2+1
                   g=num_HB_don1(i,k)
                   gg=g+1
                   HB_don1(i,k,gg)=jj-1
                   num_HB_don1(i,k)=gg
                END IF
             END DO
          END IF
       END DO
    END DO

    ALLOCATE(salida2(num_hbs2,3))

    num_hbs2=0

    DO i=1,Ldon1
       DO j=1,num_Hdon1(i)
          DO k=1,num_HB_don1(i,j)
             num_hbs2=num_hbs2+1
             salida2(num_hbs2,1)=ind_don1(i)
             salida2(num_hbs2,2)=ind_Hdon1(i,j)
             salida2(num_hbs2,3)=HB_don1(i,j,k)
          END DO
       END DO
    END DO

    !! donor2 -same

    DO i=1,Ldon2
       ii=ind_don2(i)
       pos_d=coor1(ii,:)    
       DO j=1,Lacc1
          jj=ind_acc1(j)
          vect_aux=coor1(jj,:)-pos_d
          CALL PBC(vect_aux)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          IF (norm<r_param) THEN
             vect_aux=vect_aux/norm
             DO k=1,num_Hdon2(i)
                cos_ang=dot_product(vect_aux,vect_Hdon2(i,k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs3=num_hbs3+1
                   g=num_HB_don2(i,k)
                   gg=g+1
                   HB_don2(i,k,gg)=jj-1
                   num_HB_don2(i,k)=gg
                END IF
             END DO
          END IF
       END DO
    END DO

    ALLOCATE(salida3(num_hbs3,3))

    num_hbs3=0

    DO i=1,Ldon2
       DO j=1,num_Hdon2(i)
          DO k=1,num_HB_don2(i,j)
             num_hbs3=num_hbs3+1
             salida3(num_hbs3,1)=ind_don2(i)
             salida3(num_hbs3,2)=ind_Hdon2(i,j)
             salida3(num_hbs3,3)=HB_don2(i,j,k)
          END DO
       END DO
    END DO

    !! donor2 -acceptor2

    DEALLOCATE(HB_don2)
    ALLOCATE(HB_don2(Ldon2,maxH2,maxHBs))
    num_HB_don2=0

    DO i=1,Ldon2
       ii=ind_don2(i)
       pos_d=coor1(ii,:)    
       DO j=1,Lacc2
          jj=ind_acc2(j)
          vect_aux=coor1(jj,:)-pos_d
          CALL PBC(vect_aux)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          IF (norm<r_param) THEN
             vect_aux=vect_aux/norm
             DO k=1,num_Hdon2(i)
                cos_ang=dot_product(vect_aux,vect_Hdon2(i,k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs4=num_hbs4+1
                   g=num_HB_don2(i,k)
                   gg=g+1
                   HB_don2(i,k,gg)=jj-1
                   num_HB_don2(i,k)=gg
                END IF
             END DO
          END IF
       END DO
    END DO

    ALLOCATE(salida4(num_hbs4,3))

    num_hbs4=0

    DO i=1,Ldon2
       DO j=1,num_Hdon2(i)
          DO k=1,num_HB_don2(i,j)
             num_hbs4=num_hbs4+1
             salida4(num_hbs4,1)=ind_don2(i)
             salida4(num_hbs4,2)=ind_Hdon2(i,j)
             salida4(num_hbs4,3)=HB_don2(i,j,k)
          END DO
       END DO
    END DO


    DEALLOCATE(HB_don1,HB_don2)

  END SUBROUTINE SAME_SET



  SUBROUTINE PBC(vector)
    
    implicit none
    integer::i
    DOUBLE PRECISION,dimension(3),intent(INOUT)::vector
    
    DO i=1,3
       IF (abs(vector(i))>Lbox_2(i,i)) THEN
          IF (vector(i)>Lbox_2(i,i)) THEN
             vector(i)=vector(i)-Lbox(i,i)
          ELSE
             vector(i)=vector(i)+Lbox(i,i)
          END IF
       END IF
    END DO
    
  END SUBROUTINE PBC



END MODULE HBONDS
