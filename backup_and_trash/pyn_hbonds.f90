MODULE HBONDS

  !! Indices:

  INTEGER::num_hbs1,num_hbs2,num_hbs3,num_hbs4
  INTEGER,DIMENSION(:,:),ALLOCATABLE::salida1,salida2,salida3,salida4
  INTEGER::Ldon1,Ldon2,Lacc1,Lacc2,maxH1,maxH2
  INTEGER,DIMENSION(:),ALLOCATABLE::ind_don1,ind_don2,ind_acc1,ind_acc2,num_Hdon1,num_Hdon2
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)::coor_don1,coor_don2,coor_acc1,coor_acc2
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:)::coor_Hdon1,coor_Hdon2
  DOUBLE PRECISION,DIMENSION(3)::Lbox,Lbox2
  DOUBLE PRECISION::pi
  !PARAMETER (pi=acos(-1.0d0))

CONTAINS

  SUBROUTINE INITIALIZE_PBC(box)

    DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
    INTEGER::i

    DO i=1,3
       Lbox(i)=box(i,i)
    END DO
    Lbox2=Lbox/2.0d0


  END SUBROUTINE INITIALIZE_PBC

  SUBROUTINE INITIALIZE_MEMORY()

    IF (Ldon1>0) ALLOCATE(ind_don1(Ldon1),num_Hdon1(Ldon1),coor_don1(Ldon1,3),coor_Hdon1(Ldon1,maxH1,3))
    IF (Ldon2>0) ALLOCATE(ind_don2(Ldon2),num_Hdon2(Ldon2),coor_don2(Ldon2,3),coor_Hdon2(Ldon2,maxH2,3))
    IF (Lacc1>0) ALLOCATE(ind_acc1(Lacc1),coor_acc1(Lacc1,3))
    IF (Lacc2>0) ALLOCATE(ind_acc2(Lacc2),coor_acc2(Lacc2,3))

        
    pi=acos(-1.0d0)

  END SUBROUTINE INITIALIZE_MEMORY

  SUBROUTINE FREE_MEMORY()
    
    IF (Ldon1>0) DEALLOCATE(ind_don1,num_Hdon1,coor_don1,coor_Hdon1)
    IF (Ldon2>0) DEALLOCATE(ind_don2,num_Hdon2,coor_don2,coor_Hdon2)
    IF (Lacc1>0) DEALLOCATE(ind_acc1,coor_acc1)
    IF (Lacc2>0) DEALLOCATE(ind_acc2,coor_acc2)

    IF (ALLOCATED(salida1)) DEALLOCATE(salida1)
    IF (ALLOCATED(salida2)) DEALLOCATE(salida2)
    IF (ALLOCATED(salida3)) DEALLOCATE(salida3)
    IF (ALLOCATED(salida4)) DEALLOCATE(salida4)

  END SUBROUTINE FREE_MEMORY


  SUBROUTINE DIFF_SET (r_param,ang_param)

    DOUBLE PRECISION,INTENT(IN)::r_param,ang_param
    INTEGER::i,j,k,l,ii,jj,g,gg

    TYPE array_pointer
       INTEGER,DIMENSION(:),POINTER::p1
    END TYPE array_pointer


    INTEGER,DIMENSION(Ldon1,maxH1)::num_HB_don1
    INTEGER,DIMENSION(Ldon2,maxH2)::num_HB_don2
    TYPE(array_pointer),DIMENSION(:,:),POINTER::HB_don1
    TYPE(array_pointer),DIMENSION(:,:),POINTER::HB_don2
    DOUBLE PRECISION,DIMENSION(maxH1,3)::vect_Hdon1
    DOUBLE PRECISION,DIMENSION(maxH2,3)::vect_Hdon2
    DOUBLE PRECISION,DIMENSION(3)::vect_aux,pos_d
    INTEGER,DIMENSION(:),ALLOCATABLE::aux_box


    ALLOCATE(HB_don1(Ldon1,maxH1))
    ALLOCATE(HB_don2(Ldon2,maxH2))

    num_hbs1=0
    num_hbs2=0

    num_HB_don1=0
    num_HB_don2=0

    aux_cos=COSD(ang_param)


    DO i=1,Ldon1

       pos_d=coor_don1(i,:)
       DO k=1,num_Hdon1(i)
          vect_aux=coor_Hdon1(i,k,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon1(k,:)=vect_aux/norm
       END DO

       DO j=1,Lacc2
          vect_aux=coor_acc2(j,:)-pos_d
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
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don1(i,k)%p1(h)
                   END DO
                   aux_box(gg)=j
                   IF (g>0) DEALLOCATE(HB_don1(i,k)%p1)
                   ALLOCATE(HB_don1(i,k)%p1(gg))
                   HB_don1(i,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
                   num_HB_don1(i,k)=gg
                END IF
             END DO
          END IF
       END DO

    END DO

    DO i=1,Ldon2

       pos_d=coor_don2(i,:)
       DO k=1,num_Hdon2(i)
          vect_aux=coor_Hdon2(i,k,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon2(k,:)=vect_aux/norm
       END DO

       DO j=1,Lacc1
          vect_aux=coor_acc1(j,:)-pos_d
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
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don2(i,k)%p1(h)
                   END DO
                   aux_box(gg)=j
                   IF (g>0) DEALLOCATE(HB_don2(i,k)%p1)
                   ALLOCATE(HB_don2(i,k)%p1(gg))
                   HB_don2(i,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
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
             salida1(num_hbs1,2)=j-1
             salida1(num_hbs1,3)=ind_acc2(HB_don1(i,j)%p1(k))
          END DO
       END DO
    END DO

    num_hbs2=0

    DO i=1,Ldon2
       DO j=1,num_Hdon2(i)
          DO k=1,num_HB_don2(i,j)
             num_hbs2=num_hbs2+1
             salida2(num_hbs2,1)=ind_don2(i)
             salida2(num_hbs2,2)=j-1
             salida2(num_hbs2,3)=ind_acc1(HB_don2(i,j)%p1(k))
          END DO
       END DO
    END DO

  END SUBROUTINE DIFF_SET

  SUBROUTINE SAME_SET (r_param,ang_param)


    ! I should get rid of the pointers.
    DOUBLE PRECISION,INTENT(IN)::r_param,ang_param
    INTEGER::i,j,k,l,ii,jj,g,gg

    TYPE array_pointer
       INTEGER,DIMENSION(:),POINTER::p1
    END TYPE array_pointer


    INTEGER,DIMENSION(Ldon1,maxH1)::num_HB_don1
    INTEGER,DIMENSION(Ldon2,maxH2)::num_HB_don2
    TYPE(array_pointer),DIMENSION(:,:),POINTER::HB_don1
    TYPE(array_pointer),DIMENSION(:,:),POINTER::HB_don2
    DOUBLE PRECISION,DIMENSION(Ldon1,maxH1,3)::vect_Hdon1
    DOUBLE PRECISION,DIMENSION(Ldon2,maxH2,3)::vect_Hdon2
    DOUBLE PRECISION,DIMENSION(3)::vect_aux,pos_d
    INTEGER,DIMENSION(:),ALLOCATABLE::aux_box


    num_hbs1=0
    num_hbs2=0
    num_hbs3=0
    num_hbs4=0

    aux_cos=COSD(ang_param)

    DO i=1,Ldon1
       pos_d=coor_don1(i,:)
       DO k=1,num_Hdon1(i)
          vect_aux=coor_Hdon1(i,k,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon1(i,k,:)=vect_aux/norm
       END DO
    END DO
    DO i=1,Ldon2
       pos_d=coor_don2(i,:)
       DO k=1,num_Hdon2(i)
          vect_aux=coor_Hdon2(i,k,:)-pos_d(:)
          norm=sqrt(dot_product(vect_aux,vect_aux))
          vect_Hdon2(i,k,:)=vect_aux/norm
       END DO
    END DO


    ALLOCATE(HB_don1(Ldon1,maxH1))
    ALLOCATE(HB_don2(Ldon2,maxH2))

    num_HB_don1=0
    num_HB_don2=0

    !! Same Same

    DO i=1,Ldon1
       pos_d=coor_don1(i,:)
       DO j=i+1,Ldon1
          vect_aux=coor_don1(j,:)-pos_d
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
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don1(i,k)%p1(h)
                   END DO
                   aux_box(gg)=j
                   IF (g>0) DEALLOCATE(HB_don1(i,k)%p1)
                   ALLOCATE(HB_don1(i,k)%p1(gg))
                   HB_don1(i,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
                   num_HB_don1(i,k)=gg
                END IF
             END DO
             DO k=1,num_Hdon1(j)
                cos_ang=dot_product(-vect_aux,vect_Hdon1(j,k,:))
                IF (aux_cos<cos_ang) THEN
                   num_hbs1=num_hbs1+1
                   g=num_HB_don1(j,k)
                   gg=g+1
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don1(j,k)%p1(h)
                   END DO
                   aux_box(gg)=i
                   IF (g>0) DEALLOCATE(HB_don1(j,k)%p1)
                   ALLOCATE(HB_don1(j,k)%p1(gg))
                   HB_don1(j,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
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
             salida1(num_hbs1,2)=j-1
             salida1(num_hbs1,3)=ind_acc1(HB_don1(i,j)%p1(k))
          END DO
       END DO
    END DO

    !! Same - acceptor2

    DEALLOCATE(HB_don1)
    ALLOCATE(HB_don1(Ldon1,maxH1))

    num_HB_don1=0

    DO i=1,Ldon1
       pos_d=coor_don1(i,:)    
       DO j=1,Lacc2
          vect_aux=coor_acc2(j,:)-pos_d
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
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don1(i,k)%p1(h)
                   END DO
                   aux_box(gg)=j
                   IF (g>0) DEALLOCATE(HB_don1(i,k)%p1)
                   ALLOCATE(HB_don1(i,k)%p1(gg))
                   HB_don1(i,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
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
             salida2(num_hbs2,2)=j-1
             salida2(num_hbs2,3)=ind_acc2(HB_don1(i,j)%p1(k))
          END DO
       END DO
    END DO

    !! donor2 -same

    DO i=1,Ldon2
       pos_d=coor_don2(i,:)    
       DO j=1,Lacc1
          vect_aux=coor_acc1(j,:)-pos_d
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
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don2(i,k)%p1(h)
                   END DO
                   aux_box(gg)=j
                   IF (g>0) DEALLOCATE(HB_don2(i,k)%p1)
                   ALLOCATE(HB_don2(i,k)%p1(gg))
                   HB_don2(i,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
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
             salida3(num_hbs3,2)=j-1
             salida3(num_hbs3,3)=ind_acc1(HB_don2(i,j)%p1(k))
          END DO
       END DO
    END DO

    !! donor2 -acceptor2

    DEALLOCATE(HB_don2)
    ALLOCATE(HB_don2(Ldon2,maxH2))
    num_HB_don2=0

    DO i=1,Ldon2
       pos_d=coor_don2(i,:)    
       DO j=1,Lacc2
          vect_aux=coor_acc2(j,:)-pos_d
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
                   ALLOCATE(aux_box(gg))
                   DO h=1,g
                      aux_box(h)=HB_don2(i,k)%p1(h)
                   END DO
                   aux_box(gg)=j
                   IF (g>0) DEALLOCATE(HB_don2(i,k)%p1)
                   ALLOCATE(HB_don2(i,k)%p1(gg))
                   HB_don2(i,k)%p1(:)=aux_box(:)
                   DEALLOCATE(aux_box)
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
             salida4(num_hbs4,2)=j-1
             salida4(num_hbs4,3)=ind_acc2(HB_don2(i,j)%p1(k))
          END DO
       END DO
    END DO



  END SUBROUTINE SAME_SET



  SUBROUTINE PBC(vector)
    
    implicit none
    integer::i
    DOUBLE PRECISION,dimension(3),intent(INOUT)::vector
    
    DO i=1,3
       IF (abs(vector(i))>Lbox2(i)) THEN
          IF (vector(i)>Lbox2(i)) THEN
             vector(i)=vector(i)-Lbox(i)
          ELSE
             vector(i)=vector(i)+Lbox(i)
          END IF
       END IF
    END DO
    
  END SUBROUTINE PBC



END MODULE HBONDS
