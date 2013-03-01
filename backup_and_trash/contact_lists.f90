
SUBROUTINE MAKE_CONTACT_LIST (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2

  INTEGER,DIMENSION(:),ALLOCATABLE::ilist1,ilist2
  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,ai,aj
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val
  
  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num

  IF (ALLOCATED(cl_ind)) DEALLOCATE(cl_ind)
  IF (ALLOCATED(cl_val)) DEALLOCATE(cl_val)
  IF (ALLOCATED(cl_start)) DEALLOCATE(cl_start)

  ALLOCATE(ilist1(n1),ilist2(n2))
  ilist1(:)=list1(:)+1
  ilist2(:)=list2(:)+1

  cut_off2=cut_off*cut_off

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0


  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       IF (gg==0) THEN
                          ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                       ELSE
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(:)=pcl_ind(ii)%i1(:)
                          box_val(:)=pcl_val(ii)%d1(:)
                          DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                          ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                          pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                          pcl_val(ii)%d1(1:gg)=box_val(:)
                          DEALLOCATE(box_ind,box_val)
                       END IF
                       gg=gg+1
                       pcl_ind(ii)%i1(gg)=aj
                       pcl_val(ii)%d1(gg)=val_aux
                    END IF
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       IF (gg==0) THEN
                          ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                       ELSE
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(:)=pcl_ind(ii)%i1(:)
                          box_val(:)=pcl_val(ii)%d1(:)
                          DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                          ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                          pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                          pcl_val(ii)%d1(1:gg)=box_val(:)
                          DEALLOCATE(box_ind,box_val)
                       END IF
                       gg=gg+1
                       pcl_ind(ii)%i1(gg)=aj
                       pcl_val(ii)%d1(gg)=val_aux
                    END IF
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        END IF
     ELSE
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=pcl_num(ii)
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    ll=pcl_num(jj)
                    IF (gg==0) THEN
                       ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(:)=pcl_ind(ii)%i1(:)
                       box_val(:)=pcl_val(ii)%d1(:)
                       DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                       ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                       pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                       pcl_val(ii)%d1(1:gg)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    IF (ll==0) THEN
                       ALLOCATE(pcl_ind(jj)%i1(1),pcl_val(jj)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(ll),box_val(ll))
                       box_ind(:)=pcl_ind(jj)%i1(:)
                       box_val(:)=pcl_val(jj)%d1(:)
                       DEALLOCATE(pcl_ind(jj)%i1,pcl_val(jj)%d1)
                       ALLOCATE(pcl_ind(jj)%i1(ll+1),pcl_val(jj)%d1(ll+1))
                       pcl_ind(jj)%i1(1:ll)=box_ind(:) 
                       pcl_val(jj)%d1(1:ll)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    gg=gg+1
                    ll=ll+1
                    pcl_ind(ii)%i1(gg)=aj
                    pcl_val(ii)%d1(gg)=val_aux
                    pcl_ind(jj)%i1(ll)=ai
                    pcl_val(jj)%d1(ll)=val_aux
                    pcl_num(jj)=ll
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=pcl_num(ii)
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    ll=pcl_num(jj)
                    IF (gg==0) THEN
                       ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(:)=pcl_ind(ii)%i1(:)
                       box_val(:)=pcl_val(ii)%d1(:)
                       DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                       ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                       pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                       pcl_val(ii)%d1(1:gg)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    IF (ll==0) THEN
                       ALLOCATE(pcl_ind(jj)%i1(1),pcl_val(jj)%d1(1))
                    ELSE
                       ALLOCATE(box_ind(ll),box_val(ll))
                       box_ind(:)=pcl_ind(jj)%i1(:)
                       box_val(:)=pcl_val(jj)%d1(:)
                       DEALLOCATE(pcl_ind(jj)%i1,pcl_val(jj)%d1)
                       ALLOCATE(pcl_ind(jj)%i1(ll+1),pcl_val(jj)%d1(ll+1))
                       pcl_ind(jj)%i1(1:ll)=box_ind(:) 
                       pcl_val(jj)%d1(1:ll)=box_val(:)
                       DEALLOCATE(box_ind,box_val)
                    END IF
                    gg=gg+1
                    ll=ll+1
                    pcl_ind(ii)%i1(gg)=aj
                    pcl_val(ii)%d1(gg)=val_aux
                    pcl_ind(jj)%i1(ll)=ai
                    pcl_val(jj)%d1(ll)=val_aux
                    pcl_num(jj)=ll
                 END IF
              END DO
              pcl_num(ii)=gg
           END DO
        END IF
     END IF
  END IF

  gg=SUM(pcl_num(:),DIM=1)
  ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))

  gg=0
  DO ii=1,natom1
     cl_start(ii)=gg
     DO jj=1,pcl_num(ii)
        gg=gg+1
        cl_ind(gg)=pcl_ind(ii)%i1(jj)
        cl_val(gg)=pcl_val(ii)%d1(jj)
     END DO
     DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
  END DO
  cl_start(natom1+1)=gg

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)
  DEALLOCATE(ilist1,ilist2)

  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF

END SUBROUTINE MAKE_CONTACT_LIST


SUBROUTINE MAKE_CONTACT_LIST2 (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2

  INTEGER,DIMENSION(:),ALLOCATABLE::ilist1,ilist2
  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,kk,ai,aj,lim
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind,box_ind2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val,box_val2
  
  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num,pcl_numi

  IF (ALLOCATED(cl_ind)) DEALLOCATE(cl_ind)
  IF (ALLOCATED(cl_val)) DEALLOCATE(cl_val)
  IF (ALLOCATED(cl_start)) DEALLOCATE(cl_start)

  ALLOCATE(ilist1(n1),ilist2(n2))
  ilist1(:)=list1(:)+1
  ilist2(:)=list2(:)+1

  cut_off2=cut_off*cut_off

  lim=30
  ALLOCATE(box_ind(lim),box_val(lim))

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0

  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       gg=gg+1
                       IF (gg>lim) THEN
                          ALLOCATE(box_ind2(lim),box_val2(lim))
                          box_ind2=box_ind
                          box_val2=box_val
                          DEALLOCATE(box_ind,box_val)
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(1:lim)=box_ind2
                          box_val(1:lim)=box_val2
                          DEALLOCATE(box_ind2,box_val2)
                          lim=gg
                       END IF
                       box_ind(gg)=aj
                       box_val(gg)=val_aux
                    END IF
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    IF (ai/=aj) THEN
                       gg=gg+1
                       IF (gg>lim) THEN
                          ALLOCATE(box_ind2(lim),box_val2(lim))
                          box_ind2=box_ind
                          box_val2=box_val
                          DEALLOCATE(box_ind,box_val)
                          ALLOCATE(box_ind(gg),box_val(gg))
                          box_ind(1:lim)=box_ind2
                          box_val(1:lim)=box_val2
                          DEALLOCATE(box_ind2,box_val2)
                          lim=gg
                       END IF
                       box_ind(gg)=aj
                       box_val(gg)=val_aux
                    END IF
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        END IF
     ELSE
        ALLOCATE(pcl_numi(natom1))
        pcl_numi=0
        IF (pbc_opt==1) THEN
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 CALL PBC (vect,box1,ortho1)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    gg=gg+1
                    IF (gg>lim) THEN
                       ALLOCATE(box_ind2(lim),box_val2(lim))
                       box_ind2=box_ind
                       box_val2=box_val
                       DEALLOCATE(box_ind,box_val)
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(1:lim)=box_ind2
                       box_val(1:lim)=box_val2
                       DEALLOCATE(box_ind2,box_val2)
                       lim=gg
                    END IF
                    box_ind(gg)=aj
                    box_val(gg)=val_aux
                    pcl_numi(aj)=pcl_numi(aj)+1
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        ELSE
           DO ii=1,n1
              ai=ilist1(ii)
              vect_aux=coors1(ai,:)
              gg=0
              DO jj=ii+1,n2
                 aj=ilist2(jj)
                 vect=(coors2(aj,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=cut_off2) THEN
                    gg=gg+1
                    IF (gg>lim) THEN
                       ALLOCATE(box_ind2(lim),box_val2(lim))
                       box_ind2=box_ind
                       box_val2=box_val
                       DEALLOCATE(box_ind,box_val)
                       ALLOCATE(box_ind(gg),box_val(gg))
                       box_ind(1:lim)=box_ind2
                       box_val(1:lim)=box_val2
                       DEALLOCATE(box_ind2,box_val2)
                       lim=gg
                    END IF
                    box_ind(gg)=aj
                    box_val(gg)=val_aux
                    pcl_numi(aj)=pcl_numi(aj)+1
                 END IF
              END DO
              ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
              pcl_ind(ii)%i1(:)=box_ind(1:gg)
              pcl_val(ii)%d1(:)=box_val(1:gg)
              pcl_num(ii)=gg
           END DO
        END IF
     END IF
  END IF

  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        gg=SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
        END DO
        cl_start(natom1+1)=gg
        DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
     ELSE
        gg=2*SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           gg=gg+pcl_numi(ii)
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
        END DO
        cl_start(natom1+1)=gg
        pcl_numi=0
        DO ii=1,natom1
           DO jj=1,pcl_num(ii)
              kk=pcl_ind(ii)%i1(jj)
              ll=pcl_numi(kk)+1
              pcl_numi(kk)=ll
              gg=cl_start(kk)+ll
              cl_ind(gg)=ii
              cl_val(gg)=pcl_val(ii)%d1(jj)
           END DO
           DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
        END DO
        DEALLOCATE(pcl_numi)
     END IF
  END IF

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)
  DEALLOCATE(ilist1,ilist2)

  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF


END SUBROUTINE MAKE_CONTACT_LIST2


SUBROUTINE UPDATE_CONTACT_LIST2 (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2


  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,mm,kk,ai,aj,dim_caja,amigo,amigo2,lim
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind,box_ind2,caja,aux_caja
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val,box_val2
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind,pcl_ind2   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num,pcl_numi

  cut_off2=cut_off*cut_off
  dim_caja=30
  lim=30


  ALLOCATE(box_ind(lim),box_val(lim))
  ALLOCATE(filtro(natom2),caja(dim_caja))
  filtro=.FALSE.

  cl_ind=cl_ind+1
  DEALLOCATE(cl_val)

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0


  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=0
           mm=0

           filtro(ii)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 mm=mm+1
                 IF (mm>=dim_caja) THEN
                    ALLOCATE(aux_caja(dim_caja))
                    aux_caja=caja
                    DEALLOCATE(caja)
                    ALLOCATE(caja(mm))
                    caja(1:dim_caja)=aux_caja(:)
                    DEALLOCATE(aux_caja)
                    dim_caja=mm
                 END IF
                 caja(mm)=amigo
                 filtro(amigo)=.TRUE.
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo2
                    filtro(amigo2)=.TRUE.
                 END IF
              END DO
           END DO

           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 gg=gg+1
                 IF (gg>lim) THEN
                    ALLOCATE(box_ind2(lim),box_val2(lim))
                    box_ind2=box_ind
                    box_val2=box_val
                    DEALLOCATE(box_ind,box_val)
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(1:lim)=box_ind2
                    box_val(1:lim)=box_val2
                    DEALLOCATE(box_ind2,box_val2)
                    lim=gg
                 END IF
                 box_ind(gg)=aj
                 box_val(gg)=val_aux
              END IF
           END DO
           ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
           pcl_ind(ii)%i1(:)=box_ind(1:gg)
           pcl_val(ii)%d1(:)=box_val(1:gg)
           pcl_num(ii)=gg
        END DO

     ELSE

        ALLOCATE(pcl_numi(natom1))
        pcl_numi=0

        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=0
           mm=0

           filtro(ai)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 IF (amigo>ai) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo
                    filtro(amigo)=.TRUE.
                 END IF
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    IF (amigo2>ai) THEN
                       mm=mm+1
                       IF (mm>=dim_caja) THEN
                          ALLOCATE(aux_caja(dim_caja))
                          aux_caja=caja
                          DEALLOCATE(caja)
                          ALLOCATE(caja(mm))
                          caja(1:dim_caja)=aux_caja(:)
                          DEALLOCATE(aux_caja)
                          dim_caja=mm
                       END IF
                       caja(mm)=amigo2
                       filtro(amigo2)=.TRUE.
                    END IF
                 END IF
              END DO
           END DO
           
           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 gg=gg+1
                 IF (gg>lim) THEN
                    ALLOCATE(box_ind2(lim),box_val2(lim))
                    box_ind2=box_ind
                    box_val2=box_val
                    DEALLOCATE(box_ind,box_val)
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(1:lim)=box_ind2
                    box_val(1:lim)=box_val2
                    DEALLOCATE(box_ind2,box_val2)
                    lim=gg
                 END IF
                 box_ind(gg)=aj
                 box_val(gg)=val_aux
                 pcl_numi(aj)=pcl_numi(aj)+1
              END IF
           END DO
           ALLOCATE(pcl_ind(ii)%i1(gg),pcl_val(ii)%d1(gg))
           pcl_ind(ii)%i1(:)=box_ind(1:gg)
           pcl_val(ii)%d1(:)=box_val(1:gg)
           pcl_num(ii)=gg
        END DO

     END IF
  END IF

  DEALLOCATE(filtro,caja,box_ind,box_val)
  DEALLOCATE(cl_ind,cl_start)
  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        gg=SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
        END DO
        cl_start(natom1+1)=gg
        DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
     ELSE
        gg=2*SUM(pcl_num(:),DIM=1)
        ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))
        gg=0
        DO ii=1,natom1
           cl_start(ii)=gg
           gg=gg+pcl_numi(ii)
           jj=gg+1
           gg=gg+pcl_num(ii)
           cl_val(jj:gg)=pcl_val(ii)%d1(:)
           cl_ind(jj:gg)=pcl_ind(ii)%i1(:)
        END DO
        cl_start(natom1+1)=gg
        pcl_numi=0
        DO ii=1,natom1
           DO jj=1,pcl_num(ii)
              kk=pcl_ind(ii)%i1(jj)
              ll=pcl_numi(kk)+1
              pcl_numi(kk)=ll
              gg=cl_start(kk)+ll
              cl_ind(gg)=ii
              cl_val(gg)=pcl_val(ii)%d1(jj)
           END DO
           DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
        END DO
        DEALLOCATE(pcl_numi)
     END IF
  END IF

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)


  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF


END SUBROUTINE UPDATE_CONTACT_LIST2


SUBROUTINE UPDATE_CONTACT_LIST (cut_off,sqrt_opt,diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2)

  IMPLICIT NONE

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  DOUBLE PRECISION,INTENT(IN)::cut_off
  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1,sqrt_opt
  INTEGER,INTENT(IN)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2


  DOUBLE PRECISION::cut_off2,dist2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux

  INTEGER::ii,jj,gg,ll,mm,kk,ai,aj,dim_caja,amigo,amigo2
  INTEGER,DIMENSION(:),ALLOCATABLE::box_ind,caja,aux_caja
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::box_val
  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  TYPE(iarray_pointer),DIMENSION(:),POINTER::pcl_ind,pcl_ind2   !! En indices fortran
  TYPE(darray_pointer),DIMENSION(:),POINTER::pcl_val
  INTEGER,DIMENSION(:),ALLOCATABLE::pcl_num,pcl_num2

  cut_off2=cut_off*cut_off
  dim_caja=10

  ALLOCATE(filtro(natom2),caja(dim_caja))
  filtro=.FALSE.

  cl_ind=cl_ind+1
  DEALLOCATE(cl_val)

  ALLOCATE(pcl_ind(natom1),pcl_val(natom1),pcl_num(natom1))
  pcl_num=0


  IF (diff_syst==1) THEN
     PRINT*,'NOT IMPLEMENTED YET'
  ELSE
     IF (diff_set==1) THEN
        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=0
           
           mm=0
           filtro(ii)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 mm=mm+1
                 IF (mm>=dim_caja) THEN
                    ALLOCATE(aux_caja(dim_caja))
                    aux_caja=caja
                    DEALLOCATE(caja)
                    ALLOCATE(caja(mm))
                    caja(1:dim_caja)=aux_caja(:)
                    DEALLOCATE(aux_caja)
                    dim_caja=mm
                 END IF
                 caja(mm)=amigo
                 filtro(amigo)=.TRUE.
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo2
                    filtro(amigo2)=.TRUE.
                 END IF
              END DO
           END DO

           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 IF (gg==0) THEN
                    ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                 ELSE
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(:)=pcl_ind(ii)%i1(:)
                    box_val(:)=pcl_val(ii)%d1(:)
                    DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                    ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                    pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                    pcl_val(ii)%d1(1:gg)=box_val(:)
                    DEALLOCATE(box_ind,box_val)
                 END IF
                 gg=gg+1
                 pcl_ind(ii)%i1(gg)=aj
                 pcl_val(ii)%d1(gg)=val_aux
              END IF
           END DO
           pcl_num(ii)=gg
        END DO

     ELSE

        DO ii=1,n1
           ai=list1(ii)+1
           vect_aux=coors1(ai,:)
           gg=pcl_num(ii)

           mm=0
           filtro(ai)=.TRUE.
           DO jj=cl_start(ii)+1,cl_start(ii+1)
              amigo=cl_ind(jj)
              IF (filtro(amigo)==.FALSE.) THEN
                 IF (amigo>ai) THEN
                    mm=mm+1
                    IF (mm>=dim_caja) THEN
                       ALLOCATE(aux_caja(dim_caja))
                       aux_caja=caja
                       DEALLOCATE(caja)
                       ALLOCATE(caja(mm))
                       caja(1:dim_caja)=aux_caja(:)
                       DEALLOCATE(aux_caja)
                       dim_caja=mm
                    END IF
                    caja(mm)=amigo
                    filtro(amigo)=.TRUE.
                 END IF
              END IF
              DO kk=cl_start(amigo)+1,cl_start(amigo+1)
                 amigo2=cl_ind(kk)
                 IF (filtro(amigo2)==.FALSE.) THEN
                    IF (amigo2>ai) THEN
                       mm=mm+1
                       IF (mm>=dim_caja) THEN
                          ALLOCATE(aux_caja(dim_caja))
                          aux_caja=caja
                          DEALLOCATE(caja)
                          ALLOCATE(caja(mm))
                          caja(1:dim_caja)=aux_caja(:)
                          DEALLOCATE(aux_caja)
                          dim_caja=mm
                       END IF
                       caja(mm)=amigo2
                       filtro(amigo2)=.TRUE.
                    END IF
                 END IF
              END DO
           END DO

           filtro(ai)=.FALSE.
           DO jj=1,mm
              aj=caja(jj)
              filtro(aj)=.FALSE.
              vect=(coors2(aj,:)-vect_aux)
              IF (pbc_opt) CALL PBC (vect,box1,ortho1)
              val_aux=dot_product(vect,vect)
              IF (val_aux<=cut_off2) THEN
                 ll=pcl_num(aj)
                 IF (gg==0) THEN
                    ALLOCATE(pcl_ind(ii)%i1(1),pcl_val(ii)%d1(1))
                 ELSE
                    ALLOCATE(box_ind(gg),box_val(gg))
                    box_ind(:)=pcl_ind(ii)%i1(:)
                    box_val(:)=pcl_val(ii)%d1(:)
                    DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
                    ALLOCATE(pcl_ind(ii)%i1(gg+1),pcl_val(ii)%d1(gg+1))
                    pcl_ind(ii)%i1(1:gg)=box_ind(:) 
                    pcl_val(ii)%d1(1:gg)=box_val(:)
                    DEALLOCATE(box_ind,box_val)
                 END IF
                 IF (ll==0) THEN
                    ALLOCATE(pcl_ind(aj)%i1(1),pcl_val(aj)%d1(1))
                 ELSE
                    ALLOCATE(box_ind(ll),box_val(ll))
                    box_ind(:)=pcl_ind(aj)%i1(:)
                    box_val(:)=pcl_val(aj)%d1(:)
                    DEALLOCATE(pcl_ind(aj)%i1,pcl_val(aj)%d1)
                    ALLOCATE(pcl_ind(aj)%i1(ll+1),pcl_val(aj)%d1(ll+1))
                    pcl_ind(aj)%i1(1:ll)=box_ind(:) 
                    pcl_val(aj)%d1(1:ll)=box_val(:)
                    DEALLOCATE(box_ind,box_val)
                 END IF
                 gg=gg+1
                 ll=ll+1
                 pcl_ind(ii)%i1(gg)=aj
                 pcl_val(ii)%d1(gg)=val_aux
                 pcl_ind(aj)%i1(ll)=ai
                 pcl_val(aj)%d1(ll)=val_aux
                 pcl_num(aj)=ll
              END IF
           END DO
           pcl_num(ii)=gg
        END DO

     END IF
  END IF

  gg=SUM(pcl_num(:),DIM=1)

  DEALLOCATE(filtro,caja)
  DEALLOCATE(cl_ind,cl_start)
  ALLOCATE(cl_val(gg),cl_ind(gg),cl_start(natom1+1))

  gg=0
  DO ii=1,natom1
     cl_start(ii)=gg
     DO jj=1,pcl_num(ii)
        gg=gg+1
        cl_ind(gg)=pcl_ind(ii)%i1(jj)
        cl_val(gg)=pcl_val(ii)%d1(jj)
     END DO
     DEALLOCATE(pcl_ind(ii)%i1,pcl_val(ii)%d1)
  END DO
  cl_start(natom1+1)=gg

  DEALLOCATE(pcl_ind,pcl_val,pcl_num)

  cl_ind=cl_ind-1

  IF (sqrt_opt) THEN
     cl_val=sqrt(cl_val)
  END IF

END SUBROUTINE UPDATE_CONTACT_LIST


