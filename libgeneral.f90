MODULE GLOB

  !! Syntaxis
  !! ctct:      contact
  !! ctcts :    contacts
  !! set:       set
  !! verl:      verlet
  !! ns:        neighbors
  !! gd:        grid
  !! gdns:      grid neighbors
  !! upd:       update
  !! filt:      filter
  !! cell:      cell
  !! ctoff:     cutoff
  !! ats:       atoms
  !! num:       number
  !! numat:     number of atoms
  !! glob:      global magnitude
  !! diff:      different

  INTEGER,DIMENSION(:),ALLOCATABLE::cl_ind,cl_start  !! indices fortran
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::cl_val


  INTEGER::mx_c,my_c,mz_c            !num_cells in x,y,z
  DOUBLE PRECISION::lx_c,ly_c,lz_c   !length cells in x,y,z
  INTEGER::mtot_c,mx_my_c            !num total cells, var aux
  INTEGER::nscxcell                  !number neighbor cells per cell

  LOGICAL,DIMENSION(:),ALLOCATABLE::cell_upd          !cell needs update
  INTEGER,DIMENSION(:,:),ALLOCATABLE::nsc_cell        !cell neighbors of a cell
  LOGICAL,DIMENSION(:),ALLOCATABLE::cell_pbc          !cell with none or any pbc requirement
  LOGICAL,DIMENSION(:,:),ALLOCATABLE::nsc_cell_pbc    !neighbor cell requires pbc

  INTEGER,DIMENSION(:),ALLOCATABLE::gns_head,gns_list,gns_at_cell   ! grid neighbors -atoms-

  !### Verlet neighbours lists should be this way, but f2py does not support derived types yet. 
  !TYPE I_NO_RECT_MAT
  !   INTEGER,DIMENSION(:),ALLOCATABLE::column
  !   INTEGER::dim
  !END TYPE I_NO_RECT_MAT
  !TYPE(I_NO_RECT_MAT),DIMENSION(:),ALLOCATABLE::ver_ic_ind,ver_oc_ind  !! indices fortran

  INTEGER,DIMENSION(:,:),ALLOCATABLE::ver_ic_ind,ver_oc_ind  !! indices fortran
  INTEGER,DIMENSION(:),ALLOCATABLE::ver_ic_dim,ver_oc_dim
  INTEGER::dim_ic,dim_oc

  LOGICAL,DIMENSION(:,:),ALLOCATABLE::filt_sets_ns_ind


  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::pos_ant

  INTEGER,DIMENSION(:,:),ALLOCATABLE::hbs_out
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::hbs_vals_out

CONTAINS

SUBROUTINE PERPENDICULAR_NORMED_VECT (vect1,vect2,vect_out)

  DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::vect1,vect2
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::vect_out

  CALL PRODUCT_VECT(vect1,vect2,vect_out)
  CALL NORMALIZE_VECT(vect_out)

END SUBROUTINE PERPENDICULAR_NORMED_VECT

SUBROUTINE PRODUCT_VECT(a,b,normal)

  DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::a,b
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::normal
  
  normal(1)=a(2)*b(3)-a(3)*b(2)
  normal(2)=-a(1)*b(3)+a(3)*b(1)
  normal(3)=a(1)*b(2)-a(2)*b(1)

END SUBROUTINE PRODUCT_VECT

SUBROUTINE NORMALIZE_VECT (a)

  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::a
  DOUBLE PRECISION::norm

  norm=sqrt(dot_product(a,a))
  a=a/norm

END SUBROUTINE NORMALIZE_VECT


SUBROUTINE PBC(vector,box,ortho)
 
  IMPLICIT NONE
 
  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::vector
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  INTEGER,INTENT(IN)::ortho
  INTEGER::i
  DOUBLE PRECISION::x,L,Lhalf
 
  IF (ortho==1) THEN
     DO i=1,3
        L=box(i,i)
        Lhalf=0.50d0*L
        x=vector(i)
        IF (abs(x)>Lhalf) THEN
           IF (x>Lhalf) THEN
              x=x-L
           ELSE
              x=x+L
           END IF
           vector(i)=x
        END IF
     END DO
  ELSE
 
     print*, 'Not implemented'
 
  END IF
  
END SUBROUTINE PBC


SUBROUTINE CENTER (pbc_opt,list_com,list_mov,coors,box,ortho,numat_com,numat_mov,numat_glob)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::pbc_opt,numat_com,numat_mov,numat_glob,ortho
  DOUBLE PRECISION,DIMENSION(numat_glob,3),INTENT(INOUT)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  INTEGER,DIMENSION(numat_com),INTENT(IN)::list_com
  INTEGER,DIMENSION(numat_mov),INTENT(IN)::list_mov

  INTEGER::ii,jj
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::com

  ALLOCATE(com(3))
  com=0.0d0

  DO ii=1,numat_com
     jj=list_com(ii)+1
     com(:)=com(:)+coors(jj,:)
  END DO
  com(:)=com(:)/(numat_com*1.0d0)

  com(1)=box(1,1)/2.0d0-com(1)
  com(2)=box(2,2)/2.0d0-com(2)
  com(3)=box(3,3)/2.0d0-com(3)


  DO ii=1,numat_glob
     jj=list_mov(ii)+1
     coors(jj,:)=coors(jj,:)+com(:)
  END DO

  DEALLOCATE(com)

END SUBROUTINE CENTER


INTEGER FUNCTION CELL_INDEX(ix,iy,iz)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::ix,iy,iz

  cell_index=1+MOD(ix+mx_c,mx_c)+MOD(iy+my_c,my_c)*mx_c+MOD(iz+mz_c,mz_c)*mx_my_c

END FUNCTION CELL_INDEX

LOGICAL FUNCTION CELLS_PBC(icell,ncell)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::icell,ncell
  INTEGER::ix,iy,iz
  INTEGER::nx,ny,nz
  INTEGER::ll,lx,lxy

  ll=icell-1
  ix=MOD(ll,mx_c)
  lx=INT(ll/mx_c)
  iy=MOD(lx,my_c)
  lxy=INT(ll/mx_my_c)
  iz=MOD(lxy,mz_c)

  ll=ncell-1
  nx=MOD(ll,mx_c)
  lx=INT(ll/mx_c)
  ny=MOD(lx,my_c)
  lxy=INT(ll/mx_my_c)
  nz=MOD(lxy,mz_c)


  cells_pbc=.FALSE.
  IF ((1.0d0*(abs(nx-ix)+1))>((1.0d0*mx_c)/2.0d0)) cells_pbc=.TRUE.
  IF ((1.0d0*(abs(ny-iy)+1))>((1.0d0*my_c)/2.0d0)) cells_pbc=.TRUE.
  IF ((1.0d0*(abs(nz-iz)+1))>((1.0d0*mz_c)/2.0d0)) cells_pbc=.TRUE.

END FUNCTION CELLS_PBC

LOGICAL FUNCTION CHECK_CELL(rcut,ix,iy,iz)

  IMPLICIT NONE

  DOUBLE PRECISION,INTENT(IN)::rcut
  INTEGER,INTENT(IN)::ix,iy,iz
  INTEGER::aix,aiy,aiz
  INTEGER::ii,jj,kk,pp,qq,rr
  DOUBLE PRECISION,DIMENSION(3)::o1,o2,vect
  DOUBLE PRECISION::val_aux,rcut2
  LOGICAL::filter

  check_cell=.FALSE.
  rcut2=rcut*rcut
  aix=abs(ix)
  aiy=abs(iy)
  aiz=abs(iz)

  filter=.FALSE.
  DO ii=0,1
     DO jj=0,1
        DO kk=0,1
           o1(:)=(/ii*lx_c,jj*ly_c,kk*lz_c/)
           DO pp=0,1
              DO qq=0,1
                 DO rr=0,1
                    o2(:)=(/pp*lx_c,qq*ly_c,rr*lz_c/)
                    vect=o2-o1
                    val_aux=dot_product(vect,vect)
                    IF (val_aux<rcut2) THEN
                       filter=.TRUE.
                       EXIT
                    END IF
                 END DO
                 IF (filter.eqv..TRUE.) EXIT
              END DO
              IF (filter.eqv..TRUE.) EXIT
           END DO
           IF (filter.eqv..TRUE.) EXIT
        END DO
        IF (filter.eqv..TRUE.) EXIT
     END DO
     IF (filter.eqv..TRUE.) EXIT
  END DO

  check_cell=filter

END FUNCTION CHECK_CELL

SUBROUTINE MAKE_CELL_NS (rcell,rcut,box,natom)

  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,INTENT(IN)::rcell,rcut
  double precision,DIMENSION(3,3),INTENT(IN)::box

  INTEGER::ii,jj,kk,gg,ix,iy,iz,deltx,delty,deltz
  INTEGER::icell,ncell
  LOGICAL::filter,filter2
  LOGICAL,DIMENSION(:,:,:),ALLOCATABLE::aux_mask

  mx_c=FLOOR(box(1,1)/rcell)
  lx_c=box(1,1)/mx_c
  my_c=FLOOR(box(2,2)/rcell)
  ly_c=box(2,2)/my_c
  mz_c=FLOOR(box(3,3)/rcell)
  lz_c=box(3,3)/mz_c
  mx_my_c=mx_c*my_c
  mtot_c=mx_my_c*mz_c

  !print*,mx_c,my_c,mz_c
  !print*,box(1,1)/mx_c,box(2,2)/my_c,box(3,3)/mz_c
  ! 
  !print*,((natom*1.0d0)/(mtot_c*1.0d0))

  deltx=CEILING(rcut/lx_c)
  delty=CEILING(rcut/ly_c)
  deltz=CEILING(rcut/lz_c)

  !print*,mtot_c
  !print*,mx_c,my_c,mz_c
  !print*,lx_c,ly_c,lz_c
  !print*,'================='
  !print*,'( -',deltx,',',deltx,' )' 
  !print*,'( -',delty,',',delty,' )' 
  !print*,'( -',deltz,',',deltz,' )' 
  
  ALLOCATE(aux_mask(-deltx:deltx,-delty:delty,-deltz:deltz))
  aux_mask=.FALSE.

  gg=0
  DO ii=-deltx,deltx
     DO jj=-delty,delty
        DO kk=-deltz,deltz
           IF (check_cell(rcut,ii,jj,kk)) THEN
              aux_mask(ii,jj,kk)=.TRUE.
              gg=gg+1
           END IF
        END DO
     END DO
  END DO

  nscxcell=gg

  !print*,nscxcell,'out of',(2*deltx+1)*(2*delty+1)*(2*deltz+1)

  IF (ALLOCATED(nsc_cell)) DEALLOCATE(nsc_cell)
  IF (ALLOCATED(cell_upd)) DEALLOCATE(cell_upd)
  IF (ALLOCATED(cell_pbc)) DEALLOCATE(cell_pbc)
  IF (ALLOCATED(nsc_cell_pbc)) DEALLOCATE(nsc_cell_pbc)

  ALLOCATE(nsc_cell(mtot_c,nscxcell),nsc_cell_pbc(mtot_c,nscxcell))
  ALLOCATE(cell_upd(mtot_c),cell_pbc(mtot_c))
  cell_pbc=.FALSE.

  DO ix=0,mx_c-1
     DO iy=0,my_c-1
        DO iz=0,mz_c-1
           icell=cell_index(ix,iy,iz)
           filter=.false.
           gg=0
           DO ii=-deltx,deltx
              DO jj=-delty,delty
                 DO kk=-deltz,deltz
                    IF (aux_mask(ii,jj,kk).eqv..TRUE.) THEN
                       gg=gg+1
                       ncell=cell_index(ix+ii,iy+jj,iz+kk)
                       nsc_cell(icell,gg)=ncell
                       filter2=cells_pbc(icell,ncell)
                       IF (filter2.eqv..TRUE.) filter=.TRUE.
                       nsc_cell_pbc(icell,gg)=filter2
                    END IF
                 END DO
              END DO
           END DO
           cell_pbc(icell)=filter
        END DO
     END DO
  END DO
  
  !PRINT*,COUNT(cell_pbc),'cells with pbc'

  DEALLOCATE(aux_mask)
  IF (ALLOCATED(gns_head))    DEALLOCATE(gns_head)
  IF (ALLOCATED(gns_list))    DEALLOCATE(gns_list)
  IF (ALLOCATED(gns_at_cell)) DEALLOCATE(gns_at_cell)
  ALLOCATE(gns_head(mtot_c))
  ALLOCATE(gns_list(natom))
  ALLOCATE(gns_at_cell(natom))

END SUBROUTINE MAKE_CELL_NS


SUBROUTINE GRID_NS_LIST (coors,box,natom)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,DIMENSION(natom,3),INTENT(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box

  DOUBLE PRECISION::mx_c_L,my_c_L,mz_c_L
  INTEGER::ii,icell

  mx_c_L=(mx_c*1.0d0)/box(1,1)
  my_c_L=(my_c*1.0d0)/box(2,2)
  mz_c_L=(mz_c*1.0d0)/box(3,3)

  gns_head=0
  DO ii=1,natom
     icell=1+int(coors(ii,1)*mx_c_L)+int(coors(ii,2)*my_c_L)*mx_c+int(coors(ii,3)*mz_c_L)*mx_my_c
     gns_at_cell(ii)=icell
     gns_list(ii)=gns_head(icell)
     gns_head(icell)=ii
  END DO

END SUBROUTINE GRID_NS_LIST

SUBROUTINE MAKE_VERLET_LIST_GRID_NS (r_ic,r_oc,pbc_opt,coors,box,vol,ortho,natom)

  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN)::r_ic,r_oc
  INTEGER,INTENT(IN)::pbc_opt,ortho
  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,INTENT(IN)::vol
  DOUBLE PRECISION,DIMENSION(natom,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box

  DOUBLE PRECISION::r_ic2,r_oc2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  DOUBLE PRECISION::val_aux,pi,fact

  !!! Hasta aqui 7.2 seg 1000 frames (3 min 25000)
  INTEGER::ii,jj,gg,ll,kk,icell,ncell
  INTEGER::gg_ic,gg_oc
  INTEGER::gg2_ic,gg2_oc
  LOGICAL::filter


  CALL GRID_NS_LIST(coors,box,natom)

  r_ic2=r_ic*r_ic
  r_oc2=r_oc*r_oc

  pi=acos(-1.0d0)
  fact=(natom/vol)*(4.0d0*pi/3.0d0)*2.0d0   !! factor 1.50 to be safe, since f2py does not support derived types
  dim_ic=CEILING(fact*(r_ic2*r_ic))+1
  dim_oc=CEILING(fact*(r_oc2*r_oc))+1

  IF (ALLOCATED(ver_ic_ind))   DEALLOCATE(ver_ic_ind)
  IF (ALLOCATED(ver_ic_dim))   DEALLOCATE(ver_ic_dim)
  IF (ALLOCATED(ver_oc_ind))   DEALLOCATE(ver_oc_ind)
  IF (ALLOCATED(ver_oc_dim))   DEALLOCATE(ver_oc_dim)
  IF (ALLOCATED(pos_ant))      DEALLOCATE(pos_ant)

  ALLOCATE(ver_ic_ind(natom,dim_ic),ver_oc_ind(natom,dim_oc))
  ALLOCATE(ver_ic_dim(natom),ver_oc_dim(natom))
  ALLOCATE(pos_ant(natom,3))
  pos_ant=coors
  ver_ic_dim=0
  ver_oc_dim=0

  DO ii=1,natom
     vect_aux=coors(ii,:)
     icell=gns_at_cell(ii)
     gg_oc=ver_oc_dim(ii)
     gg_ic=ver_ic_dim(ii)
     IF (cell_pbc(icell)) THEN
        DO jj=1,nscxcell
           ncell=nsc_cell(icell,jj)
           filter=nsc_cell_pbc(icell,jj)
           kk=gns_head(ncell)
           DO
              IF (kk==0) EXIT
              IF (kk>ii) THEN
                 vect=(coors(kk,:)-vect_aux)
                 IF (filter) CALL PBC (vect,box,ortho)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=r_oc2) THEN
                    gg_oc=gg_oc+1
                    ver_oc_ind(ii,gg_oc)=kk
                    gg2_oc=ver_oc_dim(kk)+1
                    ver_oc_ind(kk,gg2_oc)=ii
                    ver_oc_dim(kk)=gg2_oc
                    IF (val_aux<=r_ic2) THEN
                       gg_ic=gg_ic+1
                       ver_ic_ind(ii,gg_ic)=kk
                       gg2_ic=ver_ic_dim(kk)+1
                       ver_ic_ind(kk,gg2_ic)=ii
                       ver_ic_dim(kk)=gg2_ic
                    END IF
                 END IF
              END IF
              kk=gns_list(kk)
           END DO
        END DO
     ELSE
        DO jj=1,nscxcell
           ncell=nsc_cell(icell,jj)
           kk=gns_head(ncell)
           DO
              IF (kk==0) EXIT
              IF (kk>ii) THEN
                 vect=(coors(kk,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=r_oc2) THEN
                    gg_oc=gg_oc+1
                    ver_oc_ind(ii,gg_oc)=kk
                    gg2_oc=ver_oc_dim(kk)+1
                    ver_oc_ind(kk,gg2_oc)=ii
                    ver_oc_dim(kk)=gg2_oc
                    IF (val_aux<=r_ic2) THEN
                       gg_ic=gg_ic+1
                       ver_ic_ind(ii,gg_ic)=kk
                       gg2_ic=ver_ic_dim(kk)+1
                       ver_ic_ind(kk,gg2_ic)=ii
                       ver_ic_dim(kk)=gg2_ic
                    END IF
                 END IF
              END IF
              kk=gns_list(kk)
           END DO
        END DO
     END IF
     ver_oc_dim(ii)=gg_oc
     ver_ic_dim(ii)=gg_ic

  END DO


END SUBROUTINE MAKE_VERLET_LIST_GRID_NS



SUBROUTINE UPDATE_VERLET_LIST_GRID_NS (r_ic,r_oc,pbc_opt,coors,box,vol,ortho,natom)

  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN)::r_ic,r_oc
  INTEGER,INTENT(IN)::pbc_opt,ortho
  INTEGER,INTENT(IN)::natom
  DOUBLE PRECISION,INTENT(IN)::vol
  DOUBLE PRECISION,DIMENSION(natom,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box

  LOGICAL::update,filter
  INTEGER::ii,jj,kk,icell,ncell
  DOUBLE PRECISION::val_aux
  DOUBLE PRECISION::r_diff,drneimax,drneimax2,r_ic2,r_oc2
  DOUBLE PRECISION,DIMENSION(3)::vect,vect_aux
  INTEGER::gg_ic,gg_oc
  INTEGER::gg2_ic,gg2_oc

  update=.FALSE.
  drneimax=0.0 
  drneimax2=0.0
  r_oc2=r_oc*r_oc
  r_ic2=r_ic*r_ic
  r_diff=r_oc-r_ic


  DO ii=1,natom
     vect=pos_ant(ii,:)-coors(ii,:)
     CALL PBC (vect,box,ortho)
     val_aux=dot_product(vect,vect)
     IF (val_aux > drneimax) THEN
        drneimax2=drneimax
        drneimax=val_aux
     ELSE
        IF (val_aux > drneimax2) THEN
           drneimax2=val_aux
        END IF
     END IF
  END DO

  IF ((sqrt(drneimax)+sqrt(drneimax2))>r_diff) THEN
     update=.TRUE.
  END IF

  IF (update) THEN

     CALL GRID_NS_LIST(coors,box,natom)
     pos_ant=coors
     ver_ic_dim=0
     ver_oc_dim=0
     
     DO ii=1,natom
        vect_aux=coors(ii,:)
        icell=gns_at_cell(ii)
        gg_oc=ver_oc_dim(ii)
        gg_ic=ver_ic_dim(ii)
        IF (cell_pbc(icell)) THEN
           DO jj=1,nscxcell
              ncell=nsc_cell(icell,jj)
              filter=nsc_cell_pbc(icell,jj)
              kk=gns_head(ncell)
              DO
                 IF (kk==0) EXIT
                 IF (kk>ii) THEN
                    vect=(coors(kk,:)-vect_aux)
                    IF (filter) CALL PBC (vect,box,ortho)
                    val_aux=dot_product(vect,vect)
                    IF (val_aux<=r_oc2) THEN
                       gg_oc=gg_oc+1
                       ver_oc_ind(ii,gg_oc)=kk
                       gg2_oc=ver_oc_dim(kk)+1
                       ver_oc_ind(kk,gg2_oc)=ii
                       ver_oc_dim(kk)=gg2_oc
                       IF (val_aux<=r_ic2) THEN
                          gg_ic=gg_ic+1
                          ver_ic_ind(ii,gg_ic)=kk
                          gg2_ic=ver_ic_dim(kk)+1
                          ver_ic_ind(kk,gg2_ic)=ii
                          ver_ic_dim(kk)=gg2_ic
                       END IF
                    END IF
                 END IF
                 kk=gns_list(kk)
              END DO
           END DO
        ELSE
           DO jj=1,nscxcell
              ncell=nsc_cell(icell,jj)
              kk=gns_head(ncell)
              DO
                 IF (kk==0) EXIT
                 IF (kk>ii) THEN
                    vect=(coors(kk,:)-vect_aux)
                    val_aux=dot_product(vect,vect)
                    IF (val_aux<=r_oc2) THEN
                       gg_oc=gg_oc+1
                       ver_oc_ind(ii,gg_oc)=kk
                       gg2_oc=ver_oc_dim(kk)+1
                       ver_oc_ind(kk,gg2_oc)=ii
                       ver_oc_dim(kk)=gg2_oc
                       IF (val_aux<=r_ic2) THEN
                          gg_ic=gg_ic+1
                          ver_ic_ind(ii,gg_ic)=kk
                          gg2_ic=ver_ic_dim(kk)+1
                          ver_ic_ind(kk,gg2_ic)=ii
                          ver_ic_dim(kk)=gg2_ic
                       END IF
                    END IF
                 END IF
                 kk=gns_list(kk)
              END DO
           END DO
        END IF
        ver_oc_dim(ii)=gg_oc
        ver_ic_dim(ii)=gg_ic
        
     END DO

  ELSE

     ver_ic_dim=0
     DO ii=1,natom
        vect_aux=coors(ii,:)
        gg_ic=ver_ic_dim(ii)
        icell=gns_at_cell(ii)
        IF (cell_pbc(icell)) THEN
           DO jj=1,ver_oc_dim(ii)
              kk=ver_oc_ind(ii,jj)
              IF (kk>ii) THEN
                 vect=(coors(kk,:)-vect_aux)
                 CALL PBC (vect,box,ortho)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=r_ic2) THEN
                    gg_ic=gg_ic+1
                    ver_ic_ind(ii,gg_ic)=kk
                    gg2_ic=ver_ic_dim(kk)+1
                    ver_ic_ind(kk,gg2_ic)=ii
                    ver_ic_dim(kk)=gg2_ic
                 END IF
              END IF
           END DO
           ver_ic_dim(ii)=gg_ic
        ELSE
           DO jj=1,ver_oc_dim(ii)
              kk=ver_oc_ind(ii,jj)
              IF (kk>ii) THEN
                 vect=(coors(kk,:)-vect_aux)
                 val_aux=dot_product(vect,vect)
                 IF (val_aux<=r_ic2) THEN
                    gg_ic=gg_ic+1
                    ver_ic_ind(ii,gg_ic)=kk
                    gg2_ic=ver_ic_dim(kk)+1
                    ver_ic_ind(kk,gg2_ic)=ii
                    ver_ic_dim(kk)=gg2_ic
                 END IF
              END IF
           END DO
           ver_ic_dim(ii)=gg_ic
        END IF

     END DO

  END IF

END SUBROUTINE UPDATE_VERLET_LIST_GRID_NS

SUBROUTINE EXTRACT_NS_LIST_SETS(diff_sets,list1,list2,n1,n2,numat_glob)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_sets,n1,n2,numat_glob
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2

  INTEGER::ii,jj,ai
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter


  IF (ALLOCATED(filt_sets_ns_ind)) DEALLOCATE(filt_sets_ns_ind)
  ALLOCATE(filt_sets_ns_ind(numat_glob,dim_ic))
  filt_sets_ns_ind=.FALSE.

  IF (diff_sets==1) THEN
     
     ALLOCATE(filter(numat_glob))
     filter=.FALSE.
     DO ii=1,n2
        filter(list2(ii)+1)=.TRUE.
     END DO

     DO ii=1,n1
        ai=list1(ii)+1
        DO jj=1,ver_ic_dim(ai)
           IF (filter(ver_ic_ind(ai,jj))) THEN
              filt_sets_ns_ind(ai,jj)=.TRUE.
           END IF
        END DO
     END DO

     filter=.FALSE.
     DO ii=1,n1
        filter(list1(ii)+1)=.TRUE.
     END DO

     DO ii=1,n2
        ai=list2(ii)+1
        DO jj=1,ver_ic_dim(ai)
           IF (filter(ver_ic_ind(ai,jj))) THEN
              filt_sets_ns_ind(ai,jj)=.TRUE.
           END IF
        END DO
     END DO

     DEALLOCATE(filter)

  ELSE

     ALLOCATE(filter(numat_glob))
     filter=.FALSE.
     DO ii=1,n2
        filter(list2(ii)+1)=.TRUE.
     END DO

     DO ii=1,n1
        ai=list1(ii)+1
        DO jj=1,ver_ic_dim(ai)
           IF (filter(ver_ic_ind(ai,jj))) THEN
              filt_sets_ns_ind(ai,jj)=.TRUE.
           END IF
        END DO
     END DO

     DEALLOCATE(filter)
     
  END IF

END SUBROUTINE EXTRACT_NS_LIST_SETS

SUBROUTINE DISTANCE (diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,matrix)

IMPLICIT NONE

INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1
integer,intent(in)::n1,n2,natom1,natom2
INTEGER,DIMENSION(n1),INTENT(IN)::list1
INTEGER,DIMENSION(n2),INTENT(IN)::list2
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,dimension(natom2,3),intent(in)::coors2
double precision,dimension(n1,n2),intent(out)::matrix
integer::i,j,ai,aj
double precision,dimension(:),allocatable::vect,vect_aux
integer,dimension(:),allocatable::llist1,llist2
double precision::val_aux

ALLOCATE(vect(3),vect_aux(3))
ALLOCATE(llist1(n1),llist2(n2))
llist1=list1+1
llist2=list2+1

matrix=0.0d0

IF ((diff_syst==1) .or. (diff_set==1)) THEN
   IF (pbc_opt==1) THEN
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            vect=(coors2(aj,:)-vect_aux)
            CALL PBC (vect,box1,ortho1)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
         end do
      end do
   ELSE
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=1,n2
            aj=llist2(j)
            vect=(coors2(aj,:)-vect_aux)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
         end do
      end do
   END IF
ELSE
   IF (pbc_opt==1) THEN
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=i+1,n2
            aj=llist2(j)
            vect=(coors1(aj,:)-vect_aux)
            CALL PBC (vect,box1,ortho1)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
            matrix(j,i)=val_aux
         end do
      end do
   ELSE
      do i=1,n1
         ai=llist1(i)
         vect_aux=coors1(ai,:)
         do j=i+1,n2
            aj=llist2(j)
            vect=(coors1(aj,:)-vect_aux)
            val_aux=sqrt(dot_product(vect,vect))
            matrix(i,j)=val_aux
            matrix(j,i)=val_aux
         end do
      end do
   END IF
END IF

DEALLOCATE(vect,vect_aux)
DEALLOCATE(llist1,llist2)

END SUBROUTINE DISTANCE


SUBROUTINE DISTANCE_IMAGES (diff_syst,diff_set,list1,coors1,box1,ortho1,list2,coors2,&
                            n1,n2,natom1,natom2,min_dists,ind_atoms_min,min_image)

IMPLICIT NONE
  
INTEGER,INTENT(IN)::diff_syst,diff_set,ortho1
integer,intent(in)::n1,n2,natom1,natom2
INTEGER,DIMENSION(n1),INTENT(IN)::list1
INTEGER,DIMENSION(n2),INTENT(IN)::list2
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,dimension(natom2,3),intent(in)::coors2
double precision,dimension(n1),intent(out)::min_dists
integer,dimension(n1),intent(out):: ind_atoms_min
integer,dimension(n1,3),intent(out):: min_image

integer::ii,jj,gg,kk,ai,aj,ind_aux_min,val_imin,prov_imin
double precision,dimension(:),allocatable::vect,vect_aux,vect_aux2
integer,dimension(:),allocatable::llist1,llist2
double precision::val_aux,val_aux_min,val_ref_min,prov_val_min
INTEGER,DIMENSION(:,:),ALLOCATABLE::imag
DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::trans_imag

val_ref_min=10000.0d0
DO ii=1,3
   val_aux=dot_product(box1(ii,:),box1(ii,:))
   IF (val_ref_min>val_aux) THEN
      val_ref_min=val_aux
   END IF
END DO

ALLOCATE(vect(3),vect_aux(3),vect_aux2(3),imag(26,3),trans_imag(26,3))
gg=0
DO ii=-1,1
   DO jj=-1,1
      DO kk=-1,1
         IF ((ii/=0).or.(jj/=0).or.(kk/=0)) THEN
            gg=gg+1
            imag(gg,:)=(/ii,jj,kk/)
            trans_imag(gg,:)=ii*box1(1,:)+jj*box1(2,:)+kk*box1(3,:)
         END IF
      END DO
   END DO
END DO

ALLOCATE(llist1(n1),llist2(n2))
llist1=list1+1
llist2=list2+1

IF (diff_set==1) THEN
   DO ii=1,n1
      ai=llist1(ii)
      vect_aux=coors1(ai,:)
      val_aux_min=val_ref_min
      DO jj=1,n2
         aj=llist2(jj)
         vect_aux2=coors2(aj,:)-vect_aux
         DO gg=1,26
            vect=vect_aux2+trans_imag(gg,:)
            val_aux=dot_product(vect,vect)
            IF (val_aux<val_aux_min) THEN
               val_aux_min=val_aux
               val_imin=gg
               ind_aux_min=jj
            END IF
         END DO
      END DO
      min_dists(ii)=val_aux_min
      ind_atoms_min(ii)=ind_aux_min
      min_image(ii,1)=val_imin
   END DO
ELSE
   min_dists(:)=val_ref_min
   DO ii=1,n1
      ai=llist1(ii)
      vect_aux=coors1(ai,:)
      val_aux_min=min_dists(ii)
      DO jj=ii+1,n2
         aj=llist2(jj)
         vect_aux2=coors2(aj,:)-vect_aux
         prov_val_min=val_ref_min
         DO gg=1,26
            vect=vect_aux2+trans_imag(gg,:)
            val_aux=dot_product(vect,vect)
            IF (val_aux<prov_val_min) THEN
               prov_val_min=val_aux
               prov_imin=gg
            END IF
         END DO
         IF (prov_val_min<val_aux_min) THEN
            val_aux_min=prov_val_min
            val_imin=prov_imin
            ind_aux_min=jj
         END IF
         IF (prov_val_min<min_dists(jj)) THEN
            min_dists(jj)=prov_val_min
            ind_atoms_min(jj)=ii
            min_image(jj,1)=prov_imin
         END IF
      END DO
      min_dists(ii)=val_aux_min
      ind_atoms_min(ii)=ind_aux_min
      min_image(ii,1)=val_imin
   END DO
END IF

DO ii=1,n1
   min_dists(ii)=sqrt(min_dists(ii))
   ind_atoms_min(ii)=list2(ind_atoms_min(ii))
   min_image(ii,:)=imag(min_image(ii,1),:)
END DO

DEALLOCATE(vect,vect_aux,vect_aux2,imag,trans_imag)
DEALLOCATE(llist1,llist2)

END SUBROUTINE DISTANCE_IMAGES


SUBROUTINE RADIUS_GYRATION (list1,coors1,box1,ortho1,n1,natom1,val_Rg)

IMPLICIT NONE

INTEGER,INTENT(IN)::ortho1
integer,intent(in)::n1,natom1
INTEGER,DIMENSION(n1),INTENT(IN)::list1
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,INTENT(OUT)::val_Rg

integer::ii,ai
double precision,dimension(:),allocatable:: cdm,vect_aux
integer,dimension(:),allocatable::llist1

ALLOCATE(llist1(n1),cdm(3),vect_aux(3))
llist1=list1+1

val_Rg=0.0d0
cdm=0.0d0

DO ii=1,n1
   ai=llist1(ii)
   cdm=cdm+coors1(ai,:)
END DO
cdm=cdm/(n1*1.0d0)

DO ii=1,n1
   ai=llist1(ii)
   vect_aux=(coors1(ai,:)-cdm)
   val_Rg=val_Rg+dot_product(vect_aux,vect_aux)
END DO

DEALLOCATE(llist1,cdm,vect_aux)
val_Rg=sqrt(val_Rg/(1.0d0*n1))

END SUBROUTINE RADIUS_GYRATION


SUBROUTINE PRINCIPAL_INERTIA_AXIS (list1,coors1,box1,ortho1,n1,natom1,axis)

IMPLICIT NONE

INTEGER,INTENT(IN)::ortho1
integer,intent(in)::n1,natom1
INTEGER,DIMENSION(n1),INTENT(IN)::list1
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,DIMENSION(3,3),INTENT(OUT)::axis


integer::ii,ai
double precision:: dd
double precision,dimension(:),allocatable:: cdm,vect_aux,values
double precision,dimension(:,:),allocatable:: matrix
integer,dimension(:),allocatable::llist1
integer::num_val,info
INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work

ALLOCATE(work(8*3),iwork(5*3),ifail(3))
ALLOCATE(llist1(n1),cdm(3),vect_aux(3),matrix(3,3),values(3))

llist1=list1+1

cdm=0.0d0
matrix=0.0d0

DO ii=1,n1
   ai=llist1(ii)
   cdm=cdm+coors1(ai,:)
END DO
cdm=cdm/(n1*1.0d0)

DO ii=1,n1
   ai=llist1(ii)
   vect_aux=(coors1(ai,:)-cdm)
   dd=dot_product(vect_aux,vect_aux)
   matrix(1,1)=matrix(1,1)+dd-vect_aux(1)**2
   matrix(2,2)=matrix(2,2)+dd-vect_aux(2)**2
   matrix(3,3)=matrix(3,3)+dd-vect_aux(3)**2
   matrix(1,2)=matrix(1,2)-vect_aux(1)*vect_aux(2)
   matrix(1,3)=matrix(1,3)-vect_aux(1)*vect_aux(3)
   matrix(2,3)=matrix(2,3)-vect_aux(2)*vect_aux(3)
END DO
matrix(2,1)=matrix(1,2)
matrix(3,1)=matrix(1,3)
matrix(3,2)=matrix(2,3)


matrix=matrix/(n1*1.0d0)

DEALLOCATE(llist1)

values=0.0d0

CALL dsyevx ('V','I','U',3,matrix,3,0,0,1,3,0.0d0,num_val&
       &,values,axis,3,work,8*3,iwork,ifail,info)

IF (info/=0) THEN
   print*,"Error with the diagonalization."
   print*,"The array 'work' should has the dimension:",work(1)
END IF

DEALLOCATE(work,iwork,ifail)
DEALLOCATE(cdm,vect_aux,matrix,values)


END SUBROUTINE PRINCIPAL_INERTIA_AXIS


SUBROUTINE PRINCIPAL_GEOMETRIC_AXIS (list1,coors1,box1,ortho1,n1,natom1,axis)

IMPLICIT NONE

INTEGER,INTENT(IN)::ortho1
integer,intent(in)::n1,natom1
INTEGER,DIMENSION(n1),INTENT(IN)::list1
double precision,dimension(natom1,3),intent(in)::coors1
double precision,DIMENSION(3,3),INTENT(IN)::box1
double precision,DIMENSION(3,3),INTENT(OUT)::axis


integer::ii,ai
double precision:: dd
double precision,dimension(:),allocatable:: cdm,vect_aux,values
double precision,dimension(:,:),allocatable:: matrix
integer,dimension(:),allocatable::llist1
integer::num_val,info
INTEGER, DIMENSION(:), ALLOCATABLE::iwork,ifail
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::work

ALLOCATE(work(8*3),iwork(5*3),ifail(3))
ALLOCATE(llist1(n1),cdm(3),vect_aux(3),matrix(3,3),values(3))

llist1=list1+1

cdm=0.0d0
matrix=0.0d0

DO ii=1,n1
   ai=llist1(ii)
   cdm=cdm+coors1(ai,:)
END DO
cdm=cdm/(n1*1.0d0)

DO ii=1,n1
   ai=llist1(ii)
   vect_aux=(coors1(ai,:)-cdm)
   dd=dot_product(vect_aux,vect_aux)
   matrix(1,1)=matrix(1,1)+vect_aux(1)**2
   matrix(2,2)=matrix(2,2)+vect_aux(2)**2
   matrix(3,3)=matrix(3,3)+vect_aux(3)**2
   matrix(1,2)=matrix(1,2)+vect_aux(1)*vect_aux(2)
   matrix(1,3)=matrix(1,3)+vect_aux(1)*vect_aux(3)
   matrix(2,3)=matrix(2,3)+vect_aux(2)*vect_aux(3)
END DO
matrix(2,1)=matrix(1,2)
matrix(3,1)=matrix(1,3)
matrix(3,2)=matrix(2,3)

matrix=matrix/(n1*1.0d0)

DEALLOCATE(llist1)

values=0.0d0

CALL dsyevx ('V','I','U',3,matrix,3,0,0,1,3,0.0d0,num_val&
       &,values,axis,3,work,8*3,iwork,ifail,info)

IF (info/=0) THEN
   print*,"Error with the diagonalization."
   print*,"The array 'work' should has the dimension:",work(1)
END IF

DEALLOCATE(work,iwork,ifail)
DEALLOCATE(cdm,vect_aux,matrix,values)


END SUBROUTINE PRINCIPAL_GEOMETRIC_AXIS


SUBROUTINE DIHEDRAL_ANGLES (dih_angs,coors,box,ortho,list_angs,num_dih_angs,natom)

integer,intent(in)::ortho,natom,num_dih_angs
INTEGER,DIMENSION(num_dih_angs,4),INTENT(IN)::list_angs
double precision,dimension(natom,3),intent(in)::coors
double precision,DIMENSION(3,3),INTENT(IN)::box
DOUBLE PRECISION,DIMENSION(num_dih_angs),INTENT(OUT)::dih_angs

INTEGER::ii,jj
DOUBLE PRECISION,DIMENSION(3)::vect1,vect2,vect3
DOUBLE PRECISION::ang

DO ii=1,num_dih_angs
   vect1=coors(list_angs(ii,2)+1,:)-coors(list_angs(ii,1)+1,:)
   vect2=coors(list_angs(ii,3)+1,:)-coors(list_angs(ii,2)+1,:)
   vect3=coors(list_angs(ii,4)+1,:)-coors(list_angs(ii,3)+1,:)
   CALL PBC (vect1,box,ortho)
   CALL PBC (vect2,box,ortho)
   CALL PBC (vect3,box,ortho)
   CALL calculo_dihed (vect1,vect2,vect3,ang)
   dih_angs(ii)=ang
END DO

END SUBROUTINE DIHEDRAL_ANGLES



SUBROUTINE calculo_dihed (vec1,vec2,vec3,angulo)

  IMPLICIT NONE
  
  DOUBLE PRECISION,intent(out)::angulo
  DOUBLE PRECISION,dimension(3),intent(in)::vec1,vec2,vec3
  DOUBLE PRECISION,dimension(3)::aux1,aux2,aux3
  DOUBLE PRECISION::cosa
  integer::signo
  double precision::pi

  cosa=0.0d0
  signo=0
  angulo=0.0d0
  pi=3.14159265358979
  
  !vec1(:)=atom2(:)-atom1(:)
  !vec2(:)=atom3(:)-atom2(:)
  !vec3(:)=atom4(:)-atom3(:)
  
  aux1=0.0d0
  aux2=0.0d0
  
  CALL PRODUCT_VECT(vec1,vec2,aux1)
  CALL PRODUCT_VECT(vec2,vec3,aux2)
  
  cosa=(aux1(1)*aux2(1)+aux1(2)*aux2(2)+aux1(3)*aux2(3))/(sqrt(dot_product(aux1,aux1))*sqrt(dot_product(aux2,aux2)))

  IF (cosa>=1.0d0) THEN 
     cosa=1.0d0
  END IF
  IF (cosa<=-1.0d0) THEN 
     cosa=-1.0d0
  END IF


  cosa=(acos(cosa))!*(180.0d0/pi)

 ! IF (cosa>180.0d0) THEN
 !    print*,'ERROR EN ANGULOS'
 !    STOP
 ! END IF
  
  aux3=0.0d0
  
  CALL PRODUCT_VECT(aux1,aux2,aux3)
  
  IF ( dot_product(aux3,vec2) <= 0.0d0 ) THEN
     signo=-1
  ELSE
     signo=+1
  END IF

  angulo=cosa*signo

END SUBROUTINE calculo_dihed



SUBROUTINE neighbs_ranking (diff_syst,diff_set,pbc_opt,limit,list1,coors1,box1,ortho1,list2,coors2,&
                            n1,n2,natom1,natom2,neighb_list) !before: neighb_dist,neighb_uvect

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1
  INTEGER,INTENT(IN)::limit
  integer,intent(in)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2
  INTEGER,dimension(n1,limit),intent(out)::neighb_list

  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::dist_matrix
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dist_aux
  INTEGER,DIMENSION(:),ALLOCATABLE::neight_aux
  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
  INTEGER::ii,jj,gg

  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_dist
  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit,3),INTENT(OUT)::neighb_uvect

  ! The indexes in list1 and list2 not corrected yet (because of function dist)

  ALLOCATE(dist_matrix(n1,n2),dist_aux(n2),filter(n2),neight_aux(limit))
  filter=.TRUE.
  CALL distance (diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,dist_matrix)

  IF ((diff_syst==1) .or. (diff_set==1)) THEN
     DO ii=1,n1
        dist_aux=dist_matrix(ii,:)
        DO jj=1,limit
           gg=MINLOC(dist_aux(:),DIM=1,MASK=filter(:))
           neight_aux(jj)=gg
           filter(gg)=.FALSE.
        END DO
        DO jj=1,limit
           gg=neight_aux(jj)
           filter(gg)=.TRUE.
           neighb_list(ii,jj)=list2(gg) ! No correction was done on indexes
        END DO
     END DO
  ELSE
     DO ii=1,n1
        dist_aux=dist_matrix(ii,:)
        filter(ii)=.FALSE.
        DO jj=1,limit
           gg=MINLOC(dist_aux(:),DIM=1,MASK=filter(:))
           neight_aux(jj)=gg
           filter(gg)=.FALSE.
        END DO
        DO jj=1,limit
           gg=neight_aux(jj)
           filter(gg)=.TRUE.
           neighb_list(ii,jj)=list2(gg) ! No correction was done on indexes
        END DO
        filter(ii)=.TRUE.
     END DO
  END IF

  DEALLOCATE(dist_matrix,dist_aux,filter,neight_aux)

END SUBROUTINE neighbs_ranking


SUBROUTINE neighbs_dist (diff_syst,diff_set,pbc_opt,limit,list1,coors1,box1,ortho1,list2,&
                         coors2,n1,n2,natom1,natom2,contact_map,num_neighbs,dist_matrix) !before: neighb_dist,neighb_uvect

  IMPLICIT NONE

  INTEGER,INTENT(IN)::diff_syst,diff_set,pbc_opt,ortho1
  DOUBLE PRECISION,INTENT(IN)::limit
  integer,intent(in)::n1,n2,natom1,natom2
  INTEGER,DIMENSION(n1),INTENT(IN)::list1
  INTEGER,DIMENSION(n2),INTENT(IN)::list2
  double precision,dimension(natom1,3),intent(in)::coors1
  double precision,DIMENSION(3,3),INTENT(IN)::box1
  double precision,dimension(natom2,3),intent(in)::coors2
  INTEGER,dimension(n1,n2),intent(out)::contact_map
  INTEGER,DIMENSION(n1),INTENT(OUT)::num_neighbs
  DOUBLE PRECISION,DIMENSION(n1,n2),intent(out)::dist_matrix

  INTEGER::ii,jj,gg

  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit),INTENT(OUT)::neighb_dist
  !DOUBLE PRECISION,DIMENSION(n_atoms1,limit,3),INTENT(OUT)::neighb_uvect

  ! The indexes in list1 and list2 not corrected yet (because of function dist)
  CALL DISTANCE (diff_syst,diff_set,pbc_opt,list1,coors1,box1,ortho1,list2,coors2,n1,n2,natom1,natom2,dist_matrix)

  contact_map=0
  num_neighbs=0

  IF ((diff_syst==1) .or. (diff_set==1)) THEN
     DO ii=1,n1
        gg=0
        DO jj=1,n2
           IF (dist_matrix(ii,jj)<=limit) THEN
              contact_map(ii,jj)=1
              gg=gg+1
           END IF
        END DO
        num_neighbs(ii)=gg
     END DO
  ELSE
     DO ii=1,n1
        gg=num_neighbs(ii)
        DO jj=ii+1,n2
           IF (dist_matrix(ii,jj)<=limit) THEN
              contact_map(ii,jj)=1
              contact_map(jj,ii)=1
              gg=gg+1
              num_neighbs(jj)=num_neighbs(jj)+1
           END IF
        END DO
        num_neighbs(ii)=gg
     END DO
  END IF
  
END SUBROUTINE neighbs_dist

SUBROUTINE translate_list (sort,list,filter,distances,dim_out,n_list,trans_inds)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::n_list,dim_out,sort
  INTEGER,DIMENSION(n_list),INTENT(IN)::list,filter
  DOUBLE PRECISION,DIMENSION(n_list),INTENT(IN)::distances
  INTEGER,DIMENSION(dim_out),INTENT(OUT)::trans_inds

  LOGICAL,DIMENSION(:),ALLOCATABLE::ifilter
  INTEGER::ii,gg

  IF (sort==1) THEN

     ALLOCATE(ifilter(n_list))
     ifilter=filter

     DO ii=1,dim_out
        gg=MINLOC(distances(:),DIM=1,MASK=ifilter(:))
        ifilter(gg)=.FALSE.
        trans_inds(ii)=list(gg) ! The indexes were not corrected
     END DO

     DEALLOCATE(ifilter)

  ELSE
     gg=0
     DO ii=1,n_list
        IF (filter(ii)==1) THEN
           gg=gg+1
           trans_inds(gg)=list(ii)
        END IF
     END DO

  END IF

END SUBROUTINE translate_list


SUBROUTINE within (list_dists,cutoff,dim_list,ISIN)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::dim_list
  DOUBLE PRECISION,INTENT(IN)::cutoff
  DOUBLE PRECISION,DIMENSION(dim_list),INTENT(IN)::list_dists
  INTEGER,INTENT(OUT)::ISIN
  INTEGER::ii

  ISIN=0

  DO ii=1,dim_list
     IF (list_dists(ii)<=cutoff) THEN
        ISIN=1
        EXIT
     END IF
  END DO

END SUBROUTINE within





SUBROUTINE water_bisector (opt_pbc,list_atoms,coors,box,ortho,nlist,natom,bisectors)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::ortho,opt_pbc
  integer,intent(in)::nlist,natom
  INTEGER,DIMENSION(nlist,3),INTENT(IN)::list_atoms
  double precision,dimension(natom,3),intent(in)::coors
  double precision,DIMENSION(3,3),INTENT(IN)::box
  double precision,DIMENSION(nlist,3),INTENT(OUT)::bisectors

  INTEGER::ii,jj,ind_o,ind_h1,ind_h2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_oh1,vect_oh2,vect_aux
  DOUBLE PRECISION::val_aux

  ALLOCATE(vect_oh1(3),vect_oh2(3),vect_aux(3))

  bisectors=0.0d0

  IF (opt_pbc==1) THEN

     DO ii=1,nlist
        ind_o=list_atoms(ii,1)+1
        ind_h1=list_atoms(ii,2)+1
        ind_h2=list_atoms(ii,3)+1
        vect_oh1=coors(ind_h1,:)-coors(ind_o,:)
        vect_oh2=coors(ind_h2,:)-coors(ind_o,:)
        CALL PBC (vect_oh1,box,ortho)
        CALL PBC (vect_oh2,box,ortho)
        vect_aux=vect_oh1+vect_oh2
        val_aux=sqrt(dot_product(vect_aux,vect_aux))
        vect_aux=vect_aux/val_aux
        bisectors(ii,:)=vect_aux
     END DO

  ELSE

     DO ii=1,nlist
        ind_o=list_atoms(ii,1)+1
        ind_h1=list_atoms(ii,2)+1
        ind_h2=list_atoms(ii,3)+1
        vect_oh1=coors(ind_h1,:)-coors(ind_o,:)
        vect_oh2=coors(ind_h2,:)-coors(ind_o,:)
        vect_aux=vect_oh1+vect_oh2
        val_aux=sqrt(dot_product(vect_aux,vect_aux))
        vect_aux=vect_aux/val_aux
        bisectors(ii,:)=vect_aux
     END DO

  END IF

  DEALLOCATE(vect_oh1,vect_oh2,vect_aux)

END SUBROUTINE WATER_BISECTOR

SUBROUTINE WATER_ANGLE_BISECTOR_ATOM (opt_pbc,list_atoms,list_wats,coors,box,ortho,nwats,nlist,natoms,angles)
 
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::ortho,opt_pbc
  integer,intent(in)::natoms,nwats,nlist
  INTEGER,DIMENSION(nlist),INTENT(IN)::list_atoms
  INTEGER,DIMENSION(nwats,3),INTENT(IN)::list_wats
  double precision,dimension(natoms,3),intent(in)::coors
  double precision,DIMENSION(3,3),INTENT(IN)::box
  double precision,DIMENSION(nwats,nlist),INTENT(OUT)::angles
 
  INTEGER::ii,jj,kk,ind_o,ind_h1,ind_h2
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::vect_oh1,vect_oh2,vect_aux,vect_atomo,pos_atom
  DOUBLE PRECISION::val_aux,cosang
 
  ALLOCATE(vect_oh1(3),vect_oh2(3),vect_aux(3),vect_atomo(3),pos_atom(3))
 
  angles=0.0d0
 
  IF (opt_pbc==1) THEN

     DO kk=1,nlist

        pos_atom=coors(list_atoms(kk)+1,:)
 
        DO ii=1,nwats
           ind_o=list_wats(ii,1)+1
           ind_h1=list_wats(ii,2)+1
           ind_h2=list_wats(ii,3)+1
           vect_oh1=coors(ind_h1,:)-coors(ind_o,:)
           vect_oh2=coors(ind_h2,:)-coors(ind_o,:)
           CALL PBC (vect_oh1,box,ortho)
           CALL PBC (vect_oh2,box,ortho)
           vect_aux=vect_oh1+vect_oh2
           val_aux=sqrt(dot_product(vect_aux,vect_aux))
           vect_aux=vect_aux/val_aux
           vect_atomo=coors(ind_o,:)-pos_atom(:)
           CALL PBC (vect_atomo,box,ortho)
           val_aux=sqrt(dot_product(vect_atomo,vect_atomo))
           vect_atomo=vect_atomo/val_aux
           cosang=dot_product(vect_atomo,vect_aux)
           cosang=acos(cosang)
           angles(ii,kk)=cosang
        END DO

     END DO
 
  ELSE
 
     DO kk=1,nlist

        pos_atom=coors(list_atoms(kk)+1,:)
 
        DO ii=1,nwats
           ind_o=list_wats(ii,1)+1
           ind_h1=list_wats(ii,2)+1
           ind_h2=list_wats(ii,3)+1
           vect_oh1=coors(ind_h1,:)-coors(ind_o,:)
           vect_oh2=coors(ind_h2,:)-coors(ind_o,:)
           vect_aux=vect_oh1+vect_oh2
           val_aux=sqrt(dot_product(vect_aux,vect_aux))
           vect_aux=vect_aux/val_aux
           vect_atomo=coors(ind_o,:)-pos_atom(:)
           val_aux=sqrt(dot_product(vect_atomo,vect_atomo))
           vect_atomo=vect_atomo/val_aux
           cosang=dot_product(vect_atomo,vect_aux)
           cosang=acos(cosang)
           angles(ii,1)=cosang
        END DO
 
     END DO
 
  END IF
 
  DEALLOCATE(vect_oh1,vect_oh2,vect_aux,vect_atomo,pos_atom)

END SUBROUTINE WATER_ANGLE_BISECTOR_ATOM



!!$SUBROUTINE min_dist_atoms (pbc_opt,eq_opt,coors,box,ortho,list_a,list_b,N_tot,N_a,N_b,ind_a,ind_b,min_dist)
!!$
!!$  IMPLICIT NONE
!!$  integer,intent(in)::N_tot,N_a,N_b,pbc_opt,eq_opt,ortho
!!$  real,dimension(N_tot,3),intent(in)::coors
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box
!!$  INTEGER,DIMENSION(N_a),INTENT(IN)::list_a
!!$  INTEGER,DIMENSION(N_b),INTENT(IN)::list_b
!!$  INTEGER,DIMENSION(N_a)::auxlist_a
!!$  INTEGER,DIMENSION(N_b)::auxlist_b
!!$  INTEGER,INTENT(OUT)::ind_a,ind_b
!!$  REAL,INTENT(OUT)::min_dist
!!$
!!$  REAL,DIMENSION(3)::vect,vect_a
!!$  REAL::aux_dist
!!$  INTEGER::i,j,ia,jb
!!$
!!$  auxlist_a=list_a+1
!!$  auxlist_b=list_b+1
!!$
!!$  vect=0.0d0
!!$  min_dist=1.0d0/0.0d0
!!$  DO i=1,N_a
!!$     ia=auxlist_a(i)
!!$     vect_a=coors(ia,:)
!!$     DO j=1,N_b
!!$        jb=auxlist_b(j)
!!$        IF ((eq_opt==0).or.(jb>ia)) THEN
!!$           vect=(coors(jb,:)-vect_a(:))
!!$           IF (pbc_opt==1) CALL PBC (vect,box,ortho)
!!$           aux_dist=sqrt(dot_product(vect,vect))
!!$           IF (aux_dist<min_dist) THEN
!!$              min_dist=aux_dist
!!$              ind_a=ia
!!$              ind_b=jb
!!$           END IF
!!$        END IF
!!$     END DO
!!$  END DO
!!$  
!!$  ind_a=ind_a-1
!!$  ind_b=ind_b-1
!!$
!!$END SUBROUTINE min_dist_atoms
!!$
!!$SUBROUTINE min_dist_atoms_ref (pbc_opt,coors,box,ortho,list_a,list_coors_b,N_tot,N_a,N_b,ind_a,ind_b,min_dist)
!!$
!!$  IMPLICIT NONE
!!$  integer,intent(in)::N_tot,N_a,N_b,pbc_opt
!!$  real,dimension(N_tot,3),intent(in)::coors
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box
!!$  INTEGER,DIMENSION(N_a),INTENT(IN)::list_a
!!$  REAL,DIMENSION(N_b,3),INTENT(IN)::list_coors_b
!!$  INTEGER,DIMENSION(N_a)::auxlist_a
!!$  INTEGER,INTENT(OUT)::ind_a,ind_b
!!$  REAL,INTENT(OUT)::min_dist
!!$
!!$  REAL,DIMENSION(3)::vect,vect_a
!!$  REAL::aux_dist
!!$  INTEGER::i,j,ia,jb
!!$
!!$  auxlist_a=list_a+1
!!$
!!$  vect=0.0d0
!!$  min_dist=1.0d0/0.0d0
!!$  do i=1,N_a
!!$     ia=auxlist_a(i)
!!$     vect_a=coors(ia,:)
!!$     do j=1,N_b
!!$        vect=(list_coors_b(j,:)-vect_a(:))
!!$        IF (pbc_opt==1) CALL PBC (vect,box,ortho)
!!$        aux_dist=sqrt(dot_product(vect,vect))
!!$        IF (aux_dist<min_dist) THEN
!!$           min_dist=aux_dist
!!$           ind_a=ia
!!$           ind_b=j
!!$        END IF
!!$     end do
!!$  end do
!!$
!!$  ind_a=ind_a-1
!!$  ind_b=ind_b-1
!!$
!!$END SUBROUTINE min_dist_atoms_ref
!!$


!!$SUBROUTINE neighbs_dist2(pbc_opt,ident,ii,dist,coors1,box1,ortho1,coors2,n_atoms2,neighb_list,neighb_dist,neighb_uvect)
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,INTENT(IN)::pbc_opt,ident,ortho1
!!$  REAL,INTENT(IN)::dist
!!$  INTEGER,INTENT(IN)::n_atoms2,ii
!!$  REAL,DIMENSION(3),INTENT(IN)::coors1
!!$  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box1
!!$
!!$  INTEGER,DIMENSION(n_atoms2),INTENT(OUT)::neighb_list
!!$  REAL,DIMENSION(n_atoms2),INTENT(OUT)::neighb_dist
!!$  REAL,DIMENSION(n_atoms2,3),INTENT(OUT)::neighb_uvect
!!$
!!$  LOGICAL::lpbc,lident
!!$  INTEGER::j,g,limit
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::list
!!$  REAL,DIMENSION(:),ALLOCATABLE::list_dists
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
!!$  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
!!$  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
!!$  REAL::norm
!!$
!!$  lident=.FALSE.
!!$  lpbc=.FALSE.
!!$  IF (ident>0) lident=.TRUE.
!!$  IF (pbc_opt>0) lpbc=.TRUE.
!!$
!!$  
!!$
!!$  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))
!!$
!!$  list=0
!!$  list_dists=0.0d0
!!$  list_vects=0.0d0
!!$  aux=0.0d0
!!$  aux2=0.0d0
!!$  filter=.false.
!!$
!!$
!!$  aux=coors1(:)
!!$  limit=0
!!$
!!$  DO j=1,n_atoms2
!!$     aux2=coors2(j,:)-aux
!!$     IF (lpbc.eqv..true.) CALL PBC (aux2,box1,ortho1)
!!$     norm=sqrt(dot_product(aux2,aux2))
!!$     IF (norm<=dist) THEN
!!$        limit=limit+1
!!$        filter(j)=.true.
!!$        list_dists(j)=norm
!!$        list_vects(j,:)=aux2
!!$     END IF
!!$  END DO
!!$  
!!$  IF (lident.eqv..true.) THEN 
!!$     filter(ii)=.false.
!!$     limit=limit-1
!!$  END IF
!!$  
!!$  
!!$  DO j=1,limit
!!$     g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
!!$     list(j)=g
!!$     norm=list_dists(g)
!!$     neighb_dist(j)=norm
!!$     neighb_uvect(j,:)=list_vects(g,:)/norm
!!$     neighb_list(j)=g
!!$     filter(g)=.false.
!!$  END DO
!!$  
!!$  DO j=1,limit
!!$     g=list(j)
!!$     filter(g)=.false.
!!$  END DO
!!$  filter(ii)=.false.
!!$
!!$
!!$  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)
!!$
!!$  neighb_list=neighb_list-1
!!$
!!$
!!$END SUBROUTINE NEIGHBS_DIST2
!!$
!!$
!!$SUBROUTINE neighbs_dist1(pbc_opt,ident,dist,coors1,box1,ortho1,coors2,n_atoms1,n_atoms2,neighb_list,neighb_dist,neighb_uvect)
!!$
!!$  IMPLICIT NONE
!!$  
!!$  INTEGER,INTENT(IN)::pbc_opt,ident,ortho1
!!$  REAL,INTENT(IN)::dist
!!$  INTEGER,INTENT(IN)::n_atoms1,n_atoms2
!!$  REAL,DIMENSION(n_atoms1,3),INTENT(IN)::coors1
!!$  REAL,DIMENSION(n_atoms2,3),INTENT(IN)::coors2
!!$  REAL,DIMENSION(3,3),INTENT(IN)::box1
!!$  INTEGER,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_list
!!$  REAL,DIMENSION(n_atoms1,n_atoms2),INTENT(OUT)::neighb_dist
!!$  REAL,DIMENSION(n_atoms1,n_atoms2,3),INTENT(OUT)::neighb_uvect
!!$
!!$  LOGICAL::lpbc,lident
!!$  INTEGER::i,j,g,limit
!!$  INTEGER,DIMENSION(:),ALLOCATABLE::list
!!$  REAL,DIMENSION(:),ALLOCATABLE::list_dists
!!$  LOGICAL,DIMENSION(:),ALLOCATABLE::filter
!!$  REAL,DIMENSION(:,:),ALLOCATABLE::list_vects
!!$  REAL,DIMENSION(:),ALLOCATABLE::aux,aux2
!!$  REAL::norm
!!$
!!$  lident=.FALSE.
!!$  lpbc=.FALSE.
!!$  IF (ident>0) lident=.TRUE.
!!$  IF (pbc_opt>0) lpbc=.TRUE.
!!$
!!$
!!$  ALLOCATE(list(n_atoms2),list_vects(n_atoms2,3),list_dists(n_atoms2),filter(n_atoms2),aux(3),aux2(3))
!!$
!!$  list=0
!!$  list_dists=0.0d0
!!$  list_vects=0.0d0
!!$  aux=0.0d0
!!$  aux2=0.0d0
!!$  filter=.false.
!!$
!!$  DO i=1,n_atoms1
!!$     aux=coors1(i,:)
!!$     limit=0
!!$     DO j=1,n_atoms2
!!$        aux2=coors2(j,:)-aux
!!$        IF (lpbc.eqv..true.) CALL PBC (aux2,box1,ortho1)
!!$        norm=sqrt(dot_product(aux2,aux2))
!!$        IF (norm<=dist) THEN
!!$           limit=limit+1
!!$           filter(j)=.true.
!!$           list_dists(j)=norm
!!$           list_vects(j,:)=aux2
!!$        END IF
!!$     END DO
!!$
!!$     IF (lident.eqv..true.) THEN 
!!$        filter(i)=.false.
!!$        limit=limit-1
!!$     END IF
!!$
!!$     DO j=1,limit
!!$        g=MINLOC(list_dists(:),DIM=1,MASK=filter(:))
!!$        list(j)=g
!!$        norm=list_dists(g)
!!$        neighb_dist(i,j)=norm
!!$        neighb_uvect(i,j,:)=list_vects(g,:)/norm
!!$        neighb_list(i,j)=g
!!$        filter(g)=.false.
!!$     END DO
!!$
!!$     DO j=1,limit
!!$        g=list(j)
!!$        filter(g)=.false.
!!$     END DO
!!$     filter(i)=.false.
!!$  END DO
!!$
!!$  DEALLOCATE(list,list_vects,list_dists,filter,aux,aux2)
!!$
!!$  neighb_list=neighb_list-1
!!$
!!$END SUBROUTINE NEIGHBS_DIST1


!!$subroutine min_rmsd(struct_ref,struct_2,N,U,center_ref,center_2,rmsd,g)
!!$
!!$INTEGER,INTENT(IN)::N
!!$REAL,DIMENSION(N,3),INTENT(IN)::struct_ref,struct_2
!!$
!!$DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::U
!!$DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::center_ref,center_2
!!$DOUBLE PRECISION,INTENT(OUT)::rmsd
!!$DOUBLE PRECISION,DIMENSION(N,3),INTENT(OUT)::g
!!$
!!$
!!$INTEGER::i,j
!!$DOUBLE PRECISION,DIMENSION(N,3)::x,y
!!$DOUBLE PRECISION,DIMENSION(N)::w
!!$DOUBLE PRECISION::sw,msd,x_norm,y_norm
!!$DOUBLE PRECISION,DIMENSION(3,3)::R
!!$DOUBLE PRECISION,DIMENSION(4,4)::F
!!$DOUBLE PRECISION,DIMENSION(3)::tmp
!!$
!!$!To diagonalise:
!!$DOUBLE PRECISION,DIMENSION(4,4)::CC
!!$INTEGER::num_val,info
!!$INTEGER, DIMENSION(:),ALLOCATABLE::iwork,ifail
!!$DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE::values,work
!!$DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::vectors
!!$
!!$
!!$ALLOCATE(values(4),vectors(4,4),work(8*4),iwork(5*4),ifail(4))
!!$
!!$w=0.0d0
!!$x=0.0d0
!!$y=0.0d0
!!$CC=0.0d0
!!$rmsd=0.0d0
!!$msd=0.0d0
!!$sw=0.0d0
!!$values=0.0d0
!!$vectors=0.0d0
!!$center_ref=0.0d0
!!$center_2=0.0d0
!!$U=0.0d0
!!$g=0.0d0
!!$x_norm=0.0d0
!!$y_norm=0.0d0
!!$R=0.0d0
!!$F=0.0d0
!!$tmp=0.0d0
!!$
!!$
!!$!!! copio y peso las coordenadas:
!!$w=1.0d0
!!$DO i=1,N
!!$   sw=w(i)
!!$   x(i,:)=sw*dble(struct_ref(i,:))
!!$   y(i,:)=sw*dble(struct_2(i,:))
!!$END DO
!!$
!!$!!! calculo baricentros, centroides y normas:
!!$
!!$DO i=1,3
!!$   center_ref(i)=sum(x(:,i))/dble(N)
!!$   center_2(i)=sum(y(:,i))/dble(N)
!!$   x(:,i)=x(:,i)-center_ref(i)
!!$   y(:,i)=y(:,i)-center_2(i)
!!$   x_norm=x_norm+dot_product(x(:,i),x(:,i))
!!$   y_norm=y_norm+dot_product(y(:,i),y(:,i))
!!$END DO
!!$
!!$!!! calculo la matriz R
!!$DO i=1,3
!!$   DO j=1,3
!!$      R(i,j)=dot_product(x(:,i),y(:,j))
!!$   END DO
!!$END DO
!!$
!!$!!! construimos la matriz F:
!!$
!!$F(1,1)=R(1,1)+R(2,2)+R(3,3)
!!$F(2,1)=R(2,3)-R(3,2)
!!$F(3,1)=R(3,1)-R(1,3)
!!$F(4,1)=R(1,2)-R(2,1)
!!$F(1,2)=F(2,1)
!!$F(2,2)=R(1,1)-R(2,2)-R(3,3)
!!$F(3,2)=R(1,2)+R(2,1)
!!$F(4,2)=R(1,3)+R(3,1)
!!$F(1,3)=F(3,1)
!!$F(2,3)=F(3,2)
!!$F(3,3)=-R(1,1)+R(2,2)-R(3,3)
!!$F(4,3)=R(2,3)+R(3,2)
!!$F(1,4)=F(4,1)
!!$F(2,4)=F(4,2)
!!$F(3,4)=F(4,3)
!!$F(4,4)=-R(1,1)-R(2,2)+R(3,3) 
!!$
!!$!!! calculos los autovalores y autovectores:
!!$CC=F
!!$call dsyevx ('V','I','U',4,CC,4,0,0,1,4,0.0d0,num_val&
!!$     &,values,vectors,4,work,8*4,iwork,ifail,info)
!!$
!!$!!! computo el rmsd, la matriz de rotacion y g
!!$
!!$msd=max(0.0d0,((x_norm+y_norm)-2.0d0*values(4)))/dble(N)
!!$rmsd=sqrt(msd)
!!$
!!$
!!$call rotation_matrix(vectors(:,4),U)
!!$
!!$DO i=1,N
!!$   DO j=1,3
!!$      tmp(:)=matmul(transpose(U(:,:)),y(i,:))
!!$      g(i,j)=(x(i,j)-tmp(j))/(rmsd*dble(N))
!!$   END DO
!!$END DO
!!$
!!$!!! calculo las nuevas posiciones con la traslacion y rotacion
!!$
!!$!DO i=1,N
!!$!   pos_new(i,:)=matmul(transpose(U(:,:)),struct_2(i,:)-center_2)+center_ref
!!$!END DO
!!$
!!$
!!$END subroutine min_rmsd
!!$
!!$subroutine rotation_matrix(q, U)
!!$
!!$DOUBLE PRECISION,DIMENSION(4),INTENT(in)::q
!!$DOUBLE PRECISION,DIMENSION(3,3),INTENT(out)::U
!!$DOUBLE PRECISION::q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33
!!$
!!$q0=q(1)
!!$q1=q(2)
!!$q2=q(3)
!!$q3=q(4)
!!$
!!$b0=2.0d0*q0
!!$b1=2.0d0*q1
!!$b2=2.0d0*q2
!!$b3=2.0d0*q3
!!$
!!$q00=b0*q0-1.0d0
!!$q01=b0*q1
!!$q02=b0*q2
!!$q03=b0*q3
!!$
!!$q11=b1*q1
!!$q12=b1*q2
!!$q13=b1*q3  
!!$
!!$q22=b2*q2
!!$q23=b2*q3
!!$
!!$q33=b3*q3 
!!$
!!$U(1,1)=q00+q11
!!$U(1,2)=q12-q03
!!$U(1,3)=q13+q02
!!$
!!$U(2,1)=q12+q03
!!$U(2,2)=q00+q22
!!$U(2,3)=q23-q01
!!$
!!$U(3,1)=q13-q02
!!$U(3,2)=q23+q01
!!$U(3,3)=q00+q33
!!$
!!$end subroutine rotation_matrix
!!$
!!$
!!$subroutine rot_trans(struct,rot,center_2,center_ref,N,new_struct)
!!$
!!$IMPLICIT NONE
!!$INTEGER,INTENT(IN)::N
!!$REAL,DIMENSION(N,3),INTENT(IN)::struct
!!$DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::rot
!!$DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::center_2,center_ref
!!$REAL,DIMENSION(N,3),INTENT(OUT)::new_struct
!!$INTEGER::i
!!$
!!$DO i=1,N
!!$   new_struct(i,:)=matmul(transpose(rot(:,:)),(struct(i,:)-center_2(:)))+center_ref(:)
!!$END DO
!!$
!!$END subroutine rot_trans
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$subroutine proj3d(vect1,vect2,N,val)
!!$
!!$IMPLICIT NONE
!!$INTEGER,INTENT(IN)::N
!!$REAL,DIMENSION(N,3),INTENT(IN)::vect1,vect2
!!$REAL,DIMENSION(N,3)::vect_norm
!!$REAL::norm
!!$REAL,INTENT(OUT)::val
!!$INTEGER::i
!!$
!!$vect_norm=0.0d0
!!$norm=0.0d0
!!$DO i=1,N
!!$   norm=norm+dot_product(vect2(i,:),vect2(i,:))
!!$END DO
!!$vect_norm=vect2/(sqrt(norm))
!!$
!!$val=0.0d0
!!$DO i=1,N
!!$   val=val+dot_product(vect1(i,:),vect_norm(i,:))
!!$END DO
!!$
!!$end subroutine proj3d
!!$
!!$
!!$
END MODULE GLOB



MODULE RDF

CONTAINS

SUBROUTINE rdf_frame(distances,box,segment_min,segment_max,bins,n_A,n_B,rdf)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::n_A,n_B,bins
  DOUBLE PRECISION,INTENT(IN)::segment_min,segment_max
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box
  DOUBLE PRECISION,DIMENSION(n_A,n_B),INTENT(IN)::distances
  DOUBLE PRECISION,DIMENSION(bins),INTENT(OUT)::rdf

  INTEGER::ii,jj,kk,tt
  DOUBLE PRECISION::pi,factor,density,aux
  DOUBLE PRECISION::ro,delta_x

  pi=acos(-1.0d0)
  delta_x=(segment_max-segment_min)/(1.0d0*bins)
  factor=0.0d0
  factor=(4.0d0*pi*n_B*delta_x)
  density=(box(1,1)*box(2,2)*box(3,3))/(1.0d0*n_A)
  factor=density/factor
  rdf=0.0d0

  DO jj=1,n_A
     DO ii=1,n_B
        
        aux=(distances(jj,ii)-segment_min)/delta_x
        tt=FLOOR(aux)+1
        rdf(tt)=rdf(tt)+factor

     END DO
  END DO

  ro=segment_min+delta_x/2.0d0
  DO ii=1,bins
     rdf(ii)=rdf(ii)/(ro**2)
     ro=ro+delta_x
  END DO
  
END SUBROUTINE rdf_frame


END MODULE RDF




MODULE HBONDS

USE GLOB

!! Parameters
INTEGER::definition
DOUBLE PRECISION::sk_param,roh2_param,roo2_param
DOUBLE PRECISION::cos_angooh_param  ! the cosine


!!Output
INTEGER,DIMENSION(:),ALLOCATABLE::hbs_s_A,hbs_inds_A,hbs_s_B,hbs_inds_B
DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::hbs_vals_A,hbs_vals_B

CONTAINS

SUBROUTINE FREE_MEMORY ()

  IF (ALLOCATED(hbs_s_A)) DEALLOCATE(hbs_s_A)
  IF (ALLOCATED(hbs_s_B)) DEALLOCATE(hbs_s_B)
  IF (ALLOCATED(hbs_inds_A)) DEALLOCATE(hbs_inds_A)
  IF (ALLOCATED(hbs_inds_B)) DEALLOCATE(hbs_inds_B)
  IF (ALLOCATED(hbs_vals_A)) DEALLOCATE(hbs_vals_A)
  IF (ALLOCATED(hbs_vals_B)) DEALLOCATE(hbs_vals_B)

END SUBROUTINE FREE_MEMORY


SUBROUTINE GET_HBONDS (effic,diff_syst,diff_set,pbc_opt,acc_A,acc_sH_A,acc_H_A,don_A,don_sH_A,don_H_A,coors1,box1,ortho1, &
     acc_B,acc_sH_B,acc_H_B,don_B,don_sH_B,don_H_B,coors2,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, &
     nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H,natomA,natomB)

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  INTEGER,INTENT(IN)::effic,diff_syst,diff_set,pbc_opt,ortho1
  INTEGER,INTENT(IN)::nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H
  INTEGER,INTENT(IN)::nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H
  INTEGER,DIMENSION(nA_acc),INTENT(IN)    ::acc_A
  INTEGER,DIMENSION(nB_acc),INTENT(IN)    ::acc_B
  INTEGER,DIMENSION(nA_don),INTENT(IN)    ::don_A
  INTEGER,DIMENSION(nB_don),INTENT(IN)    ::don_B
  INTEGER,DIMENSION(nA_acc_sH),INTENT(IN) ::acc_sH_A
  INTEGER,DIMENSION(nB_acc_sH),INTENT(IN) ::acc_sH_B
  INTEGER,DIMENSION(nA_acc_H),INTENT(IN)  ::acc_H_A
  INTEGER,DIMENSION(nB_acc_H),INTENT(IN)  ::acc_H_B
  INTEGER,DIMENSION(nA_don_sH),INTENT(IN) ::don_sH_A
  INTEGER,DIMENSION(nB_don_sH),INTENT(IN) ::don_sH_B
  INTEGER,DIMENSION(nA_don_H),INTENT(IN)  ::don_H_A
  INTEGER,DIMENSION(nB_don_H),INTENT(IN)  ::don_H_B
  DOUBLE PRECISION,DIMENSION(natomA,3),intent(in)::coors1
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box1
  DOUBLE PRECISION,DIMENSION(natomB,3),intent(in)::coors2

  TYPE(iarray_pointer),DIMENSION(:),POINTER::hbs_a_ind,hbs_b_ind 
  TYPE(darray_pointer),DIMENSION(:),POINTER::hbs_a_val,hbs_b_val 
  INTEGER,DIMENSION(:),ALLOCATABLE::num_hbs_a,num_hbs_b

  INTEGER::ii,jj,gg
  INTEGER::don,don_h,acc
  INTEGER::lim_hbs

  DOUBLE PRECISION,DIMENSION(3)::pos_acc,pos_don,pos_h
  DOUBLE PRECISION,DIMENSION(3)::vect_don_acc,vect_h_acc,vect_don_h
  DOUBLE PRECISION,DIMENSION(3)::aux_vect_1,aux_vect_2,aux_vect_3
  DOUBLE PRECISION::dist_h_acc,dist_don_acc,dist_don_h,aux_cos,sk_val
  DOUBLE PRECISION::dist2_h_acc,dist2_don_acc,dist2_don_h

  INTEGER::acc_H1,acc_H2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_perp

  INTEGER,DIMENSION(:),ALLOCATABLE::aux_box_ind,aux2_box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_box_val,aux2_box_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  lim_hbs=3

  ALLOCATE(aux_box_ind(lim_hbs),aux_box_val(lim_hbs),filtro(lim_hbs))
  filtro=.FALSE.

  ALLOCATE(hbs_a_ind(nA_don_H),hbs_a_val(nA_don_H))
  ALLOCATE(hbs_b_ind(nB_don_H),hbs_b_val(nB_don_H))
  ALLOCATE(num_hbs_a(nA_don_H),num_hbs_b(nB_don_H))

  num_hbs_a=0
  num_hbs_b=0


  IF (diff_syst==0) THEN

     SELECT CASE (definition)

     CASE (1) ! Skinner                                                                                  !SK 
        !!!! Source: R. Kumar, J. R. Schmidt and J. L. Skinner. J. Chem. Phys. 126, 204107 (2007)        !SK 
        ALLOCATE(vect_perp(nB_acc,3))                                                                    !SK 
        DO ii=1,nB_acc                                                                                   !SK
           acc=acc_B(jj)+1                                                                               !SK
           jj=acc_sH_B(ii)+1                                                                             !SK
           acc_H1=acc_H_B(jj)                                                                            !SK 
           acc_H2=acc_H_B(jj+1)                                                                          !SK 
           pos_acc=coors2(acc,:)                                                                         !SK 
           aux_vect_1=coors2(acc_H1,:)-pos_acc(:)                                                        !SK 
           aux_vect_2=coors2(acc_H2,:)-pos_acc(:)                                                        !SK 
           CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)                             !SK 
           vect_perp(ii,:)=aux_vect_3                                                                    !SK 
        END DO                                                                                           !SK 
                                                                                                         !SK 
        DO ii=1,nA_don                                                                                   !SK 
           DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                                           !SK 
              don_h=don_H_A(hh)+1                                                                        !SK 
              pos_h=coors1(don_h,:)                                                                      !SK 
              gg=0                                                                                       !SK 
              DO jj=1,nB_acc                                                                             !SK
                 acc=acc_B(jj)+1                                                                         !SK 
                 vect_h_acc=pos_h(:)-coors2(acc,:)                                                       !SK 
                 IF (pbc_opt==1) CALL PBC(vect_h_acc,box1,ortho1)                                           !SK 
                 dist_h_acc=sqrt(dot_product(vect_h_acc,vect_h_acc))                                     !SK 
                 aux_cos=dot_product(vect_perp(jj,:),vect_h_acc(:))/dist_h_acc                           !SK 
                 if (aux_cos>=1.0d0) aux_cos=1.0d0                                                       !SK 
                 if (aux_cos<=-1.0d0) aux_cos=-1.0d0                                                     !SK 
                 aux_cos=acos(aux_cos)                                                                   !SK 
                 aux_cos=aux_cos*(90/pi)                                                                 !SK 
                 IF (aux_cos>90) THEN                                                                    !SK 
                    print*,'aquiii error 3.14',aux_cos                                                   !SK 
                    STOP                                                                                 !SK 
                 END IF                                                                                  !SK 
                 IF (aux_cos<0) THEN                                                                     !SK 
                    print*,'aquiii error 3.14',aux_cos                                                   !SK 
                    STOP                                                                                 !SK 
                 END IF                                                                                  !SK 
                 sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)         !SK 
                 IF (sk_val>sk_param) THEN                                                               !SK 
                    gg=gg+1                                                                              !SK 
                    IF (gg>lim_hbs) THEN                                                                 !SK 
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                             !SK 
                       aux2_box_ind=aux_box_ind                                                          !SK 
                       aux2_box_val=aux_box_val                                                          !SK 
                       DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                        !SK 
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                              !SK 
                       aux_box_ind(1:lim_hbs)=aux2_box_ind                                               !SK 
                       aux_box_val(1:lim_hbs)=aux2_box_val                                               !SK 
                       filtro=.FALSE.
                       DEALLOCATE(aux2_box_ind,aux2_box_val)                                             !SK 
                       lim_hbs=gg                                                                        !SK 
                    END IF                                                                               !SK 
                    aux_box_ind(gg)=acc_B(jj)                                                            !SK 
                    aux_box_val(gg)=sk_val                                                               !SK 
                 END IF                                                                                  !SK 
              END DO                                                                                     !SK 
                                                                                                         !SK 
              IF (gg>0) THEN                                                                             !SK 
                 ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                                     !SK 
                 filtro(1:gg)=.TRUE.                                                                     !SK 
                 DO jj=1,gg                                                                              !SK 
                    ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                          !SK 
                    filtro(ll)=.FALSE.                                                                   !SK 
                    hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                                 !SK 
                    hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                                 !SK 
                 END DO                                                                                  !SK 
              END IF                                                                                     !SK 
                                                                                                         !SK 
              num_hbs_a(hh)=gg                                                                           !SK 
                                                                                                         !SK 
           END DO                                                                                        !SK 
        END DO                                                                                           !SK 
        DEALLOCATE(vect_perp)                                                                            !SK 
                                                                                                         !SK 
        IF (diff_set==1) THEN                                                                               !SK 
                                                                                                         !SK 
           DO ii=1,nA_acc                                                                                !SK 
              acc=acc_A(jj)+1                                                                            !SK
              jj=acc_sH_A(ii)+1                                                                          !SK
              acc_H1=acc_H_A(jj)                                                                         !SK 
              acc_H2=acc_H_A(jj+1)                                                                       !SK 
              pos_acc=coors1(acc,:)                                                                      !SK 
              aux_vect_1=coors1(acc_H1,:)-pos_acc(:)                                                     !SK 
              aux_vect_2=coors1(acc_H2,:)-pos_acc(:)                                                     !SK 
              CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)                          !SK 
              vect_perp(ii,:)=aux_vect_3                                                                 !SK 
           END DO                                                                                        !SK 
                                                                                                         !SK 
           DO ii=1,nB_don                                                                                !SK 
              DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)                                                        !SK 
                 don_h=don_H_B(hh)+1                                                                     !SK 
                 pos_h=coors2(don_h,:)                                                                   !SK 
                 gg=0                                                                                    !SK 
                 DO jj=1,nA_acc                                                                          !SK 
                    acc=acc_A(jj)+1                                                                      !SK 
                    vect_h_acc=pos_h(:)-coors1(acc,:)                                                    !SK 
                    IF (pbc_opt==1) CALL PBC(vect_h_acc,box1,ortho1)                                        !SK 
                    dist_h_acc=sqrt(dot_product(vect_h_acc,vect_h_acc))                                  !SK 
                    aux_cos=dot_product(vect_perp(jj,:),vect_h_acc(:))/dist_h_acc                        !SK 
                    if (aux_cos>=1.0d0) aux_cos=1.0d0                                                    !SK 
                    if (aux_cos<=-1.0d0) aux_cos=-1.0d0                                                  !SK 
                    aux_cos=acos(aux_cos)                                                                !SK 
                    aux_cos=aux_cos*(90/pi)                                                              !SK 
                    IF (aux_cos>90) THEN                                                                 !SK 
                       print*,'aquiii error 3.14',aux_cos                                                !SK 
                       STOP                                                                              !SK 
                    END IF                                                                               !SK 
                    IF (aux_cos<0) THEN                                                                  !SK 
                       print*,'aquiii error 3.14',aux_cos                                                !SK 
                       STOP                                                                              !SK 
                    END IF                                                                               !SK 
                    sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)      !SK 
                    IF (sk_val>sk_param) THEN                                                            !SK 
                       gg=gg+1                                                                           !SK 
                       IF (gg>lim_hbs) THEN                                                              !SK 
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                          !SK 
                          aux2_box_ind=aux_box_ind                                                       !SK 
                          aux2_box_val=aux_box_val                                                       !SK 
                          DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                     !SK 
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                           !SK 
                          aux_box_ind(1:lim_hbs)=aux2_box_ind                                            !SK 
                          aux_box_val(1:lim_hbs)=aux2_box_val                                            !SK 
                          filtro=.FALSE.
                          DEALLOCATE(aux2_box_ind,aux2_box_val)                                          !SK 
                          lim_hbs=gg                                                                     !SK 
                       END IF                                                                            !SK 
                       aux_box_ind(gg)=acc_A(jj)                                                         !SK 
                       aux_box_val(gg)=sk_val                                                            !SK 
                    END IF                                                                               !SK 
                 END DO                                                                                  !SK 
                                                                                                         !SK 
                 IF (gg>0) THEN                                                                          !SK 
                    ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                                  !SK 
                    filtro(1:gg)=.TRUE.                                                                  !SK 
                    DO jj=1,gg                                                                           !SK 
                       ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                       !SK 
                       filtro(ll)=.FALSE.                                                                !SK 
                       hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                                              !SK 
                       hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                                              !SK 
                    END DO                                                                               !SK 
                 END IF                                                                                  !SK 
                                                                                                         !SK 
                 num_hbs_b(hh)=gg                                                                        !SK 
                                                                                                         !SK 
              END DO                                                                                     !SK 
           END DO                                                                                        !SK 
           DEALLOCATE(vect_perp)                                                                         !SK 
                                                                                                         !SK 
        END IF                                                                                           !SK 

     CASE (2) ! R(o,h)                                                                   !ROH 
        !!!! Source: V. J. Buch. J. Chem. Phys. 96, 3814-3823 (1992)                     !ROH 
        DO ii=1,nA_don                                                                   !ROH
           DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                           !ROH 
              don_h=don_H_A(hh)+1                                                        !ROH 
              pos_h=coors1(don_h,:)                                                      !ROH 
              gg=0                                                                       !ROH 
              DO jj=1,nB_acc                                                             !ROH 
                 acc=acc_B(jj)+1                                                         !ROH 
                 vect_h_acc=coors2(acc,:)-pos_h(:)                                       !ROH 
                 IF (pbc_opt==1) CALL PBC(vect_h_acc,box1,ortho1)                           !ROH 
                 dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)                          !ROH 
                 IF (dist2_h_acc<roh2_param) THEN                                          !ROH 
                    gg=gg+1                                                              !ROH 
                    IF (gg>lim_hbs) THEN                                                 !ROH 
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))             !ROH 
                       aux2_box_ind=aux_box_ind                                          !ROH 
                       aux2_box_val=aux_box_val                                          !ROH 
                       DEALLOCATE(aux_box_ind,aux_box_val,filtro)                        !ROH 
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))              !ROH 
                       aux_box_ind(1:lim_hbs)=aux2_box_ind                               !ROH 
                       aux_box_val(1:lim_hbs)=aux2_box_val                               !ROH 
                       filtro=.FALSE.
                       DEALLOCATE(aux2_box_ind,aux2_box_val)                             !ROH 
                       lim_hbs=gg                                                        !ROH 
                    END IF                                                               !ROH 
                    aux_box_ind(gg)=acc_B(jj)                                            !ROH 
                    aux_box_val(gg)=dist2_h_acc                                           !ROH 
                 END IF                                                                  !ROH 
              END DO                                                                     !ROH 
                                                                                         !ROH 
              IF (gg>0) THEN                                                             !ROH 
                 ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                     !ROH 
                 filtro(1:gg)=.TRUE.                                                     !ROH 
                 DO jj=1,gg                                                              !ROH 
                    ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                          !ROH 
                    filtro(ll)=.FALSE.                                                   !ROH 
                    hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                 !ROH 
                    hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                 !ROH 
                 END DO                                                                  !ROH 
              END IF                                                                     !ROH 
                                                                                         !ROH 
              num_hbs_a(hh)=gg                                                           !ROH 
                                                                                         !ROH 
           END DO                                                                        !ROH 
        END DO                                                                           !ROH 
                                                                                         !ROH 
                                                                                         !ROH 
        IF (diff_set==1) THEN                                                               !ROH 
                                                                                         !ROH 
           DO ii=1,num_don_B                                                             !ROH 
              DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)                                        !ROH 
                 don_h=don_H_B(hh)+1                                                     !ROH 
                 pos_h=coors2(don_h,:)                                                   !ROH 
                 gg=0                                                                    !ROH 
                 DO jj=1,num_acc_A                                                       !ROH 
                    acc=acc_A(jj)+1                                                      !ROH 
                    vect_h_acc=coors1(acc,:)-pos_h(:)                                    !ROH 
                    IF (pbc_opt==1) CALL PBC(vect_h_acc,box1,ortho1)                        !ROH 
                    dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)                  !ROH 
                    IF (dist2_h_acc<roh2_param) THEN                                       !ROH 
                       gg=gg+1                                                           !ROH 
                       IF (gg>lim_hbs) THEN                                              !ROH 
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))          !ROH 
                          aux2_box_ind=aux_box_ind                                       !ROH 
                          aux2_box_val=aux_box_val                                       !ROH 
                          DEALLOCATE(aux_box_ind,aux_box_val,filtro)                     !ROH 
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))           !ROH 
                          aux_box_ind(1:lim_hbs)=aux2_box_ind                            !ROH 
                          aux_box_val(1:lim_hbs)=aux2_box_val                            !ROH 
                          filtro=.FALSE.
                          DEALLOCATE(aux2_box_ind,aux2_box_val)                          !ROH 
                          lim_hbs=gg                                                     !ROH 
                       END IF                                                            !ROH 
                       aux_box_ind(gg)=acc_A(jj)                                         !ROH 
                       aux_box_val(gg)=dist2_h_acc                                        !ROH 
                    END IF                                                               !ROH 
                 END DO                                                                  !ROH 
                                                                                         !ROH 
                 IF (gg>0) THEN                                                          !ROH 
                    ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                  !ROH 
                    filtro(1:gg)=.TRUE.                                                  !ROH 
                    DO jj=1,gg                                                           !ROH 
                       ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                       !ROH 
                       filtro(ll)=.FALSE.                                                !ROH 
                       hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                              !ROH 
                       hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                              !ROH 
                    END DO                                                               !ROH 
                 END IF                                                                  !ROH 
                                                                                         !ROH 
                 num_hbs_b(hh)=gg                                                        !ROH 
                                                                                         !ROH 
              END DO                                                                     !ROH 
           END DO                                                                        !ROH 
                                                                                         !ROH 
        END IF


     CASE (3) ! R(o,o)-Ang(o,o,h)                                                                        !ROO_ANG 
        !!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)                          !ROO_ANG 
!!$        DO ii=1,nA_don                                                                                   !ROO_ANG 
!!$           don=don_A(ii)+1                                                                               !ROO_ANG 
!!$           pos_don=coors1(don,:)                                                                         !ROO_ANG 
!!$           DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                                           !ROO_ANG 
!!$              don_h=don_H_A(hh)+1                                                                        !ROO_ANG 
!!$              vect_don_h=coors1(don_h,:)-pos_don                                                         !ROO_ANG
!!$              dist2_don_h=dot_product(vect_don_h,vect_don_h)                                             !ROO_ANG
!!$              gg=0                                                                                       !ROO_ANG 
!!$              DO jj=1,nB_acc                                                                             !ROO_ANG 
!!$                 acc=acc_B(jj)+1                                                                         !ROO_ANG 
!!$                 vect_don_acc=coors2(acc,:)-pos_don(:)                                                   !ROO_ANG 
!!$                 IF (pbc_opt==1) CALL PBC(vect_don_acc,box1,ortho1)                                         !ROO_ANG 
!!$                 dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)                                    !ROO_ANG
!!$                 IF (dist2_don_acc<roo2_param) THEN                                                      !ROO_ANG 
!!$                    aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))       !ROO_ANG
!!$                    IF (aux_cos>cos_angooh_param) THEN                                                   !ROO_ANG
!!$                       gg=gg+1                                                                           !ROO_ANG 
!!$                       IF (gg>lim_hbs) THEN                                                              !ROO_ANG 
!!$                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                          !ROO_ANG 
!!$                          aux2_box_ind=aux_box_ind                                                       !ROO_ANG 
!!$                          aux2_box_val=aux_box_val                                                       !ROO_ANG 
!!$                          DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                     !ROO_ANG 
!!$                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                           !ROO_ANG 
!!$                          aux_box_ind(1:lim_hbs)=aux2_box_ind                                            !ROO_ANG 
!!$                          aux_box_val(1:lim_hbs)=aux2_box_val                                            !ROO_ANG 
!!$                          filtro=.FALSE.                                                                 !ROO_ANG
!!$                          DEALLOCATE(aux2_box_ind,aux2_box_val)                                          !ROO_ANG 
!!$                          lim_hbs=gg                                                                     !ROO_ANG 
!!$                       END IF                                                                            !ROO_ANG 
!!$                       aux_box_ind(gg)=acc_B(jj)                                                         !ROO_ANG 
!!$                       aux_box_val(gg)=dist2_don_acc                                                     !ROO_ANG 
!!$                    END IF                                                                               !ROO_ANG 
!!$                 END IF                                                                                  !ROO_ANG
!!$              END DO                                                                                     !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$              IF (gg>0) THEN                                                                             !ROO_ANG 
!!$                 ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                                     !ROO_ANG 
!!$                 filtro(1:gg)=.TRUE.                                                                     !ROO_ANG 
!!$                 DO jj=1,gg                                                                              !ROO_ANG 
!!$                    ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                          !ROO_ANG 
!!$                    filtro(ll)=.FALSE.                                                                   !ROO_ANG 
!!$                    hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                                 !ROO_ANG 
!!$                    hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                                 !ROO_ANG 
!!$                 END DO                                                                                  !ROO_ANG 
!!$              END IF                                                                                     !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$              num_hbs_a(hh)=gg                                                                           !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$           END DO                                                                                        !ROO_ANG 
!!$        END DO                                                                                           !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$        IF (diff_set==1) THEN                                                                               !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$           DO ii=1,nB_don                                                                                !ROO_ANG 
!!$              don=don_B(ii)+1                                                                            !ROO_ANG 
!!$              pos_don=coors2(don,:)                                                                      !ROO_ANG 
!!$              DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)                                                        !ROO_ANG 
!!$                 don_h=don_H_B(hh)+1                                                                     !ROO_ANG 
!!$                 vect_don_h=coors2(don_h,:)-pos_don                                                      !ROO_ANG
!!$                 dist2_don_h=dot_product(vect_don_h,vect_don_h)                                          !ROO_ANG
!!$                 gg=0                                                                                    !ROO_ANG 
!!$                 DO jj=1,nA_acc                                                                          !ROO_ANG 
!!$                    acc=acc_A(jj)+1                                                                      !ROO_ANG 
!!$                    vect_don_acc=coors1(acc,:)-pos_don(:)                                                !ROO_ANG 
!!$                    IF (pbc_opt==1) CALL PBC(vect_don_acc,box1,ortho1)                                      !ROO_ANG 
!!$                    dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)                                 !ROO_ANG
!!$                    IF (dist2_don_acc<roo2_param) THEN                                                     !ROO_ANG
!!$                       aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))    !ROO_ANG
!!$                       IF (aux_cos>cos_angooh_param) THEN                                                !ROO_ANG
!!$                          gg=gg+1                                                                        !ROO_ANG 
!!$                          IF (gg>lim_hbs) THEN                                                           !ROO_ANG 
!!$                             ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                       !ROO_ANG 
!!$                             aux2_box_ind=aux_box_ind                                                    !ROO_ANG 
!!$                             aux2_box_val=aux_box_val                                                    !ROO_ANG 
!!$                             DEALLOCATE(aux_box_ind,aux_box_val,filtro)                                  !ROO_ANG 
!!$                             ALLOCATE(aux_box_ind(gg),aux_box_val(gg),filtro(gg))                        !ROO_ANG 
!!$                             aux_box_ind(1:lim_hbs)=aux2_box_ind                                         !ROO_ANG 
!!$                             aux_box_val(1:lim_hbs)=aux2_box_val                                         !ROO_ANG 
!!$                             filtro=.FALSE.                                                              !ROO_ANG
!!$                             DEALLOCATE(aux2_box_ind,aux2_box_val)                                       !ROO_ANG 
!!$                             lim_hbs=gg                                                                  !ROO_ANG 
!!$                          END IF                                                                         !ROO_ANG 
!!$                          aux_box_ind(gg)=acc_A(jj)                                                      !ROO_ANG 
!!$                          aux_box_val(gg)=dist2_don_acc                                                  !ROO_ANG
!!$                       END IF                                                                            !ROO_ANG 
!!$                    END IF                                                                               !ROO_ANG
!!$                 END DO                                                                                  !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$                 IF (gg>0) THEN                                                                          !ROO_ANG 
!!$                    ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                                  !ROO_ANG 
!!$                    filtro(1:gg)=.TRUE.                                                                  !ROO_ANG 
!!$                    DO jj=1,gg                                                                           !ROO_ANG 
!!$                       ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                       !ROO_ANG 
!!$                       filtro(ll)=.FALSE.                                                                !ROO_ANG 
!!$                       hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                                              !ROO_ANG 
!!$                       hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                                              !ROO_ANG 
!!$                    END DO                                                                               !ROO_ANG 
!!$                 END IF                                                                                  !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$                 num_hbs_b(hh)=gg                                                                        !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$              END DO                                                                                     !ROO_ANG 
!!$           END DO                                                                                        !ROO_ANG 
!!$                                                                                                         !ROO_ANG 
!!$        END IF                                                                                           !ROO_ANG


     !!!CASE (4) ! Donor-Acceptor-Number
     !!!   !!!! Source: A. D. Hammerich, V. J. Buch. J. Chem. Phys. 128, 111101 (2008)
     !!!   print*, "NOT IMPLEMENTED YET"
     !!! 
     !!!CASE (5) ! Topological
     !!!   !!!! Source: R. H. Henchman and S. J. Irudayam. J. Phys. Chem. B. 114, 16792-16810 (2010)
     !!! 
     !!!   DO ii=1,nA_don                                                                                 
     !!!      DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)                                                                  
     !!!         don_h=don_H_A(hh)+1                                                                             
     !!!         pos_h=coors1(don_h,:)                                                                       
     !!!         gg=0
     !!!         candidato=0
     !!!         val_out=100.0d0 !! Should be NaN
     !!!         DO jj=1,nB_acc                                                                           
     !!!            acc=acc_B(jj)+1                                                                          
     !!!            vect_h_acc=coors2(acc,:)-pos_don(:)                                                    
     !!!            IF (pbc_opt==1) CALL PBC(vect_h_acc,box1,ortho1)                                          
     !!!            dist_h_acc=dot_product(vect_h_acc,vect_h_acc)
     !!!            IF (dist_don_acc<val_out) THEN
     !!!               candidato=jj
     !!!            END IF
     !!!         END DO
     !!!         
     !!! 
     !!! 
     !!! 
     !!!                  gg=gg+1                                                                            
     !!!                  IF (gg>lim_hbs) THEN                                                               
     !!!                     ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                           
     !!!                     aux2_box_ind=aux_box_ind                                                        
     !!!                     aux2_box_val=aux_box_val                                                        
     !!!                     DEALLOCATE(aux_box_ind,aux_box_val)                                             
     !!!                     ALLOCATE(aux_box_ind(gg),aux_box_val(gg))                                       
     !!!                     aux_box_ind(1:lim_hbs)=aux2_box_ind                                             
     !!!                     aux_box_val(1:lim_hbs)=aux2_box_val                                             
     !!!                     DEALLOCATE(aux2_box_ind,aux2_box_val)                                           
     !!!                     lim_hbs=gg                                                                      
     !!!                  END IF                                                                             
     !!!                  aux_box_ind(gg)=acc_B(jj)                                                          
     !!!                  aux_box_val(gg)=val_out                                                            
     !!!               END IF                                                                                
     !!!            END IF                                                                                  
     !!!         END DO                                                                                      
     !!!                                                                                                     
     !!!         IF (gg>0) THEN                                                                              
     !!!            ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))                                      
     !!!            ALLOCATE(filtro(gg))                                                                     
     !!!            filtro=.TRUE.                                                                            
     !!!            DO jj=1,gg                                                                               
     !!!               ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                           
     !!!               filtro(ll)=.FALSE.                                                                    
     !!!               hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)                                                  
     !!!               hbs_a_val(hh)%d1(jj)=aux_box_val(ll)                                                  
     !!!            END DO                                                                                   
     !!!            DEALLOCATE(filtro)                                                                       
     !!!         END IF                                                                                      
     !!!                                                                                                     
     !!!         num_hbs_a(hh)=gg                                                                            
     !!!                                                                                                     
     !!!      END DO                                                                                         
     !!!   END DO                                                                                            
     !!!                                                                                                     
     !!!                                                                                                     
     !!!   IF (diff_set==1) THEN                                                                                
     !!!                                                                                                     
     !!!      DO ii=1,num_don_B                                                                              
     !!!         don=don_B(ii,1)+1                                                                           
     !!!         pos_don=coors2(don,:)                                                                       
     !!!         DO hh=H_s_B(ii)+1,H_s_B(ii+1)                                                               
     !!!            don_h=H_B(hh)+1                                                                          
     !!!            pos_h=coors2(don_h,:)                                                                    
     !!!            vect_don_h=coors2(don_h,:)-pos_don_h                                                    
     !!!            gg=0                                                                                     
     !!!            DO jj=1,num_acc_A                                                                        
     !!!               acc=acc_A(jj)+1                                                                       
     !!!               vect_don_acc=coors1(acc,:)-pos_don(:)                                                 
     !!!               IF (pbc_opt==1) CALL PBC(vect_don_acc,box1,ortho1)                                       
     !!!               dist_don_acc=sqrt(dot_product(vect_don_acc,vect_don_acc))                            
     !!!               IF (dist_don_acc<roo_param) THEN                                                      
     !!!                  dist_don_h_acc=sqrt(dot_product(vect_don_h_acc,vect_don_h_acc))                   
     !!!                  aux_cos=dot_product(vect_don_h_acc,vect_don_acc)/(dist_don_acc*dist_don_h_acc)    
     !!!                  IF (aux_cos>cos_angooh_param) THEN                                                
     !!!                     gg=gg+1                                                                         
     !!!                     IF (gg>lim_hbs) THEN                                                            
     !!!                        ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))                        
     !!!                        aux2_box_ind=aux_box_ind                                                     
     !!!                        aux2_box_val=aux_box_val                                                     
     !!!                        DEALLOCATE(aux_box_ind,aux_box_val)                                          
     !!!                        ALLOCATE(aux_box_ind(gg),aux_box_val(gg))                                    
     !!!                        aux_box_ind(1:lim_hbs)=aux2_box_ind                                          
     !!!                        aux_box_val(1:lim_hbs)=aux2_box_val                                          
     !!!                        DEALLOCATE(aux2_box_ind,aux2_box_val)                                        
     !!!                        lim_hbs=gg                                                                   
     !!!                     END IF                                                                          
     !!!                     aux_box_ind(gg)=acc_A(jj)                                                       
     !!!                     aux_box_val(gg)=val_out                                                         
     !!!                  END IF                                                                             
     !!!               END IF                                                                               
     !!!            END DO                                                                                   
     !!!                                                                                                     
     !!!            IF (gg>0) THEN                                                                           
     !!!               ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))                                   
     !!!               ALLOCATE(filtro(gg))                                                                  
     !!!               filtro=.TRUE.                                                                         
     !!!               DO jj=1,gg                                                                            
     !!!                  ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)                                        
     !!!                  filtro(ll)=.FALSE.                                                                 
     !!!                  hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)                                               
     !!!                  hbs_b_val(hh)%d1(jj)=aux_box_val(ll)                                               
     !!!               END DO                                                                                
     !!!               DEALLOCATE(filtro)                                                                    
     !!!            END IF                                                                                   
     !!!                                                                                                     
     !!!            num_hbs_b(hh)=gg                                                                         
     !!!                                                                                                     
     !!!         END DO                                                                                      
     !!!      END DO                                                                                         
     !!!                                                                                                     
     !!!   END IF                                                                                           
     !!! 
     !!!CASE (6) ! Donor-Number-Ang(o,o,h)
     !!!   print*, "NOT IMPLEMENTED YET"
     !!! 
     !!!CASE (7) ! Nearest-Neighbour
     !!!   print*, "NOT IMPLEMENTED YET"

     CASE DEFAULT
        PRINT*, 'Error: Hbond definition unknown'
     
     END SELECT

     print*,SUM(num_hbs_a(:),DIM=1),SUM(num_hbs_b(:),DIM=1)



  ELSE
     PRINT*,'NOT IMPLEMENTED YET'
  END IF


  DEALLOCATE(aux_box_ind,aux_box_val)
  DEALLOCATE(hbs_a_ind,hbs_a_val)
  DEALLOCATE(hbs_b_ind,hbs_b_val)
  DEALLOCATE(num_hbs_a,num_hbs_b)
  DEALLOCATE(filtro)


END SUBROUTINE GET_HBONDS


SUBROUTINE GET_HBONDS_ROO_ANG_NS_LIST (diff_set,pbc_opt,acc_A,acc_sH_A,acc_H_A,don_A,don_sH_A,don_H_A,coors,box,ortho, &
     acc_B,acc_sH_B,acc_H_B,don_B,don_sH_B,don_H_B,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, &
     nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H,numat_glob)

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  INTEGER,INTENT(IN)::diff_set,pbc_opt,ortho,numat_glob
  INTEGER,INTENT(IN)::nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H
  INTEGER,INTENT(IN)::nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H
  INTEGER,DIMENSION(nA_acc),INTENT(IN)    ::acc_A
  INTEGER,DIMENSION(nB_acc),INTENT(IN)    ::acc_B
  INTEGER,DIMENSION(nA_don),INTENT(IN)    ::don_A
  INTEGER,DIMENSION(nB_don),INTENT(IN)    ::don_B
  INTEGER,DIMENSION(nA_acc_sH),INTENT(IN) ::acc_sH_A
  INTEGER,DIMENSION(nB_acc_sH),INTENT(IN) ::acc_sH_B
  INTEGER,DIMENSION(nA_acc_H),INTENT(IN)  ::acc_H_A
  INTEGER,DIMENSION(nB_acc_H),INTENT(IN)  ::acc_H_B
  INTEGER,DIMENSION(nA_don_sH),INTENT(IN) ::don_sH_A
  INTEGER,DIMENSION(nB_don_sH),INTENT(IN) ::don_sH_B
  INTEGER,DIMENSION(nA_don_H),INTENT(IN)  ::don_H_A
  INTEGER,DIMENSION(nB_don_H),INTENT(IN)  ::don_H_B
  DOUBLE PRECISION,DIMENSION(numat_glob,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box


  TYPE(iarray_pointer),DIMENSION(:),POINTER::hbs_a_ind,hbs_b_ind 
  TYPE(darray_pointer),DIMENSION(:),POINTER::hbs_a_val,hbs_b_val 
  INTEGER,DIMENSION(:),ALLOCATABLE::num_hbs_a,num_hbs_b

  INTEGER::ii,jj,gg
  INTEGER::don,don_h,acc
  INTEGER::lim_hbs

  DOUBLE PRECISION,DIMENSION(3)::pos_acc,pos_don,pos_h
  DOUBLE PRECISION,DIMENSION(3)::vect_don_acc,vect_h_acc,vect_don_h
  DOUBLE PRECISION,DIMENSION(3)::aux_vect_1,aux_vect_2,aux_vect_3
  DOUBLE PRECISION::dist_h_acc,dist_don_acc,dist_don_h,aux_cos,sk_val
  DOUBLE PRECISION::dist2_h_acc,dist2_don_acc,dist2_don_h

  INTEGER::acc_H1,acc_H2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_perp

  INTEGER,DIMENSION(:),ALLOCATABLE::aux_box_ind,aux2_box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_box_val,aux2_box_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  lim_hbs=3

  ALLOCATE(aux_box_ind(lim_hbs),aux_box_val(lim_hbs),filtro(lim_hbs))
  filtro=.FALSE.

  ALLOCATE(hbs_a_ind(nA_don_H),hbs_a_val(nA_don_H))
  ALLOCATE(hbs_b_ind(nB_don_H),hbs_b_val(nB_don_H))
  ALLOCATE(num_hbs_a(nA_don_H),num_hbs_b(nB_don_H))

  num_hbs_a=0
  num_hbs_b=0

  !!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)


  CALL EXTRACT_NS_LIST_SETS(diff_set,don_A,acc_B,nA_don,nB_acc,numat_glob)

  DO ii=1,nA_don
     don=don_A(ii)+1
     pos_don=coors(don,:)
     DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
        don_h=don_H_A(hh)+1
        vect_don_h=coors(don_h,:)-pos_don
        IF (pbc_opt==1) CALL PBC(vect_don_h,box,ortho)
        dist2_don_h=dot_product(vect_don_h,vect_don_h)
        gg=0
        DO jj=1,ver_ic_dim(don)
           IF (filt_sets_ns_ind(don,jj)) THEN
              acc=ver_ic_ind(don,jj)
              vect_don_acc=coors(acc,:)-pos_don(:)
              IF (pbc_opt==1) CALL PBC(vect_don_acc,box,ortho)
              dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)
              aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))
              IF (aux_cos>cos_angooh_param) THEN
                 gg=gg+1
                 IF (gg>lim_hbs) THEN
                    ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                    aux2_box_ind=aux_box_ind
                    aux2_box_val=aux_box_val
                    DEALLOCATE(aux_box_ind,aux_box_val)
                    ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                    aux_box_ind(1:lim_hbs)=aux2_box_ind
                    aux_box_val(1:lim_hbs)=aux2_box_val
                    DEALLOCATE(aux2_box_ind,aux2_box_val)
                    lim_hbs=gg
                 END IF
                 aux_box_ind(gg)=acc-1
                 aux_box_val(gg)=dist2_don_acc
              END IF
           END IF
        END DO
        IF (gg>0) THEN
           ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))
           filtro(1:gg)=.TRUE.
           DO jj=1,gg
              ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
              filtro(ll)=.FALSE.
              hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)
              hbs_a_val(hh)%d1(jj)=aux_box_val(ll)
           END DO
        END IF

        num_hbs_a(hh)=gg

     END DO
  END DO

  IF (diff_set==1) THEN

     CALL EXTRACT_NS_LIST_SETS(diff_set,don_B,acc_A,nB_don,nA_acc,numat_glob)

     DO ii=1,nB_don
        don=don_B(ii)+1
        pos_don=coors(don,:)
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           don_h=don_H_B(hh)+1
           vect_don_h=coors(don_h,:)-pos_don
           IF (pbc_opt==1) CALL PBC(vect_don_h,box,ortho)
           dist2_don_h=dot_product(vect_don_h,vect_don_h)
           gg=0
           DO jj=1,ver_ic_dim(don)
              IF (filt_sets_ns_ind(don,jj)) THEN
                 acc=ver_ic_ind(don,jj)
                 vect_don_acc=coors(acc,:)-pos_don(:)
                 IF (pbc_opt==1) CALL PBC(vect_don_acc,box,ortho)
                 dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)
                 aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))
                 IF (aux_cos>cos_angooh_param) THEN
                    gg=gg+1
                    IF (gg>lim_hbs) THEN
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                       aux2_box_ind=aux_box_ind
                       aux2_box_val=aux_box_val
                       DEALLOCATE(aux_box_ind,aux_box_val)
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                       aux_box_ind(1:lim_hbs)=aux2_box_ind
                       aux_box_val(1:lim_hbs)=aux2_box_val
                       DEALLOCATE(aux2_box_ind,aux2_box_val)
                       lim_hbs=gg
                    END IF
                    aux_box_ind(gg)=acc-1
                    aux_box_val(gg)=dist2_don_acc
                 END IF
              END IF
           END DO
           IF (gg>0) THEN
              ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))
              filtro(1:gg)=.TRUE.
              DO jj=1,gg
                 ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
                 filtro(ll)=.FALSE.
                 hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)
                 hbs_b_val(hh)%d1(jj)=aux_box_val(ll)
              END DO
           END IF
           
           num_hbs_b(hh)=gg
           
        END DO
     END DO
  END IF

  IF (ALLOCATED(hbs_out)) DEALLOCATE(hbs_out)
  IF (ALLOCATED(hbs_vals_out)) DEALLOCATE(hbs_vals_out)


  IF (diff_set==1) THEN
     ii=SUM(num_hbs_a)+SUM(num_hbs_b)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=sqrt(hbs_a_val(hh)%d1(jj))
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
     DO ii=1,nB_don
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           IF (num_hbs_b(hh)>0) THEN
              DO jj=1,num_hbs_b(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_B(ii)
                 hbs_out(gg,2)=don_H_B(hh)
                 hbs_out(gg,3)=hbs_b_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=sqrt(hbs_b_val(hh)%d1(jj))
              END DO
              DEALLOCATE(hbs_b_ind(hh)%i1,hbs_b_val(hh)%d1)
           END IF
        END DO
     END DO
  ELSE
     ii=SUM(num_hbs_a)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=sqrt(hbs_a_val(hh)%d1(jj))
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
  END IF

  DEALLOCATE(aux_box_ind,aux_box_val,filtro)
  DEALLOCATE(hbs_a_ind,hbs_a_val)
  DEALLOCATE(hbs_b_ind,hbs_b_val)
  DEALLOCATE(num_hbs_a,num_hbs_b)

END SUBROUTINE GET_HBONDS_ROO_ANG_NS_LIST


SUBROUTINE GET_HBONDS_ROO_ANG (diff_set,pbc_opt,acc_A,acc_sH_A,acc_H_A,don_A,don_sH_A,don_H_A,coors,box,ortho, &
     acc_B,acc_sH_B,acc_H_B,don_B,don_sH_B,don_H_B,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, &
     nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H,numat_glob)

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  INTEGER,INTENT(IN)::diff_set,pbc_opt,ortho,numat_glob
  INTEGER,INTENT(IN)::nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H
  INTEGER,INTENT(IN)::nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H
  INTEGER,DIMENSION(nA_acc),INTENT(IN)    ::acc_A
  INTEGER,DIMENSION(nB_acc),INTENT(IN)    ::acc_B
  INTEGER,DIMENSION(nA_don),INTENT(IN)    ::don_A
  INTEGER,DIMENSION(nB_don),INTENT(IN)    ::don_B
  INTEGER,DIMENSION(nA_acc_sH),INTENT(IN) ::acc_sH_A
  INTEGER,DIMENSION(nB_acc_sH),INTENT(IN) ::acc_sH_B
  INTEGER,DIMENSION(nA_acc_H),INTENT(IN)  ::acc_H_A
  INTEGER,DIMENSION(nB_acc_H),INTENT(IN)  ::acc_H_B
  INTEGER,DIMENSION(nA_don_sH),INTENT(IN) ::don_sH_A
  INTEGER,DIMENSION(nB_don_sH),INTENT(IN) ::don_sH_B
  INTEGER,DIMENSION(nA_don_H),INTENT(IN)  ::don_H_A
  INTEGER,DIMENSION(nB_don_H),INTENT(IN)  ::don_H_B
  DOUBLE PRECISION,DIMENSION(numat_glob,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box


  TYPE(iarray_pointer),DIMENSION(:),POINTER::hbs_a_ind,hbs_b_ind 
  TYPE(darray_pointer),DIMENSION(:),POINTER::hbs_a_val,hbs_b_val 
  INTEGER,DIMENSION(:),ALLOCATABLE::num_hbs_a,num_hbs_b

  INTEGER::ii,jj,gg
  INTEGER::don,don_h,acc
  INTEGER::lim_hbs

  DOUBLE PRECISION,DIMENSION(3)::pos_acc,pos_don,pos_h
  DOUBLE PRECISION,DIMENSION(3)::vect_don_acc,vect_h_acc,vect_don_h
  DOUBLE PRECISION,DIMENSION(3)::aux_vect_1,aux_vect_2,aux_vect_3
  DOUBLE PRECISION::dist_h_acc,dist_don_acc,dist_don_h,aux_cos,sk_val
  DOUBLE PRECISION::dist2_h_acc,dist2_don_acc,dist2_don_h

  INTEGER::acc_H1,acc_H2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_perp

  INTEGER,DIMENSION(:),ALLOCATABLE::aux_box_ind,aux2_box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_box_val,aux2_box_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  lim_hbs=3

  ALLOCATE(aux_box_ind(lim_hbs),aux_box_val(lim_hbs),filtro(lim_hbs))
  filtro=.FALSE.

  ALLOCATE(hbs_a_ind(nA_don_H),hbs_a_val(nA_don_H))
  ALLOCATE(hbs_b_ind(nB_don_H),hbs_b_val(nB_don_H))
  ALLOCATE(num_hbs_a(nA_don_H),num_hbs_b(nB_don_H))

  num_hbs_a=0
  num_hbs_b=0

  !!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)

  DO ii=1,nA_don
     don=don_A(ii)+1
     pos_don=coors(don,:)
     DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
        don_h=don_H_A(hh)+1
        vect_don_h=coors(don_h,:)-pos_don
        IF (pbc_opt==1) CALL PBC(vect_don_h,box,ortho)
        dist2_don_h=dot_product(vect_don_h,vect_don_h)
        gg=0
        DO jj=1,nB_acc
           acc=acc_B(jj)+1
           IF (acc/=don) THEN
              vect_don_acc=coors(acc,:)-pos_don(:)
              IF (pbc_opt==1) CALL PBC(vect_don_acc,box,ortho)
              dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)
              IF (dist2_don_acc<roo2_param) THEN
                 aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))
                 IF (aux_cos>cos_angooh_param) THEN
                    gg=gg+1
                    IF (gg>lim_hbs) THEN
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                       aux2_box_ind=aux_box_ind
                       aux2_box_val=aux_box_val
                       DEALLOCATE(aux_box_ind,aux_box_val)
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                       aux_box_ind(1:lim_hbs)=aux2_box_ind
                       aux_box_val(1:lim_hbs)=aux2_box_val
                       DEALLOCATE(aux2_box_ind,aux2_box_val)
                       lim_hbs=gg
                    END IF
                    aux_box_ind(gg)=acc-1
                    aux_box_val(gg)=dist2_don_acc
                 END IF
              END IF
           END IF
        END DO
        IF (gg>0) THEN
           ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))
           filtro(1:gg)=.TRUE.
           DO jj=1,gg
              ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
              filtro(ll)=.FALSE.
              hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)
              hbs_a_val(hh)%d1(jj)=aux_box_val(ll)
           END DO
        END IF

        num_hbs_a(hh)=gg

     END DO
  END DO

  IF (diff_set==1) THEN
     DO ii=1,nB_don
        don=don_B(ii)+1
        pos_don=coors(don,:)
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           don_h=don_H_B(hh)+1
           vect_don_h=coors(don_h,:)-pos_don
           IF (pbc_opt==1) CALL PBC(vect_don_h,box,ortho)
           dist2_don_h=dot_product(vect_don_h,vect_don_h)
           gg=0
           DO jj=1,nA_acc
              acc=acc_A(jj)+1
              vect_don_acc=coors(acc,:)-pos_don(:)
              IF (pbc_opt==1) CALL PBC(vect_don_acc,box,ortho)
              dist2_don_acc=dot_product(vect_don_acc,vect_don_acc)
              IF (dist2_don_acc<roo2_param) THEN
                 aux_cos=dot_product(vect_don_h,vect_don_acc)/(sqrt(dist2_don_acc*dist2_don_h))
                 IF (aux_cos>cos_angooh_param) THEN
                    gg=gg+1
                    IF (gg>lim_hbs) THEN
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                       aux2_box_ind=aux_box_ind
                       aux2_box_val=aux_box_val
                       DEALLOCATE(aux_box_ind,aux_box_val)
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                       aux_box_ind(1:lim_hbs)=aux2_box_ind
                       aux_box_val(1:lim_hbs)=aux2_box_val
                       DEALLOCATE(aux2_box_ind,aux2_box_val)
                       lim_hbs=gg
                    END IF
                    aux_box_ind(gg)=acc-1
                    aux_box_val(gg)=dist2_don_acc
                 END IF
              END IF
           END DO
           IF (gg>0) THEN
              ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))
              filtro(1:gg)=.TRUE.
              DO jj=1,gg
                 ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
                 filtro(ll)=.FALSE.
                 hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)
                 hbs_b_val(hh)%d1(jj)=aux_box_val(ll)
              END DO
           END IF
           
           num_hbs_b(hh)=gg
           
        END DO
     END DO
  END IF

  IF (ALLOCATED(hbs_out)) DEALLOCATE(hbs_out)
  IF (ALLOCATED(hbs_vals_out)) DEALLOCATE(hbs_vals_out)

  IF (diff_set==1) THEN
     ii=SUM(num_hbs_a)+SUM(num_hbs_b)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=sqrt(hbs_a_val(hh)%d1(jj))
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
     DO ii=1,nB_don
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           IF (num_hbs_b(hh)>0) THEN
              DO jj=1,num_hbs_b(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_B(ii)
                 hbs_out(gg,2)=don_H_B(hh)
                 hbs_out(gg,3)=hbs_b_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=sqrt(hbs_b_val(hh)%d1(jj))
              END DO
              DEALLOCATE(hbs_b_ind(hh)%i1,hbs_b_val(hh)%d1)
           END IF
        END DO
     END DO
  ELSE
     ii=SUM(num_hbs_a)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=sqrt(hbs_a_val(hh)%d1(jj))
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
  END IF

  DEALLOCATE(aux_box_ind,aux_box_val,filtro)
  DEALLOCATE(hbs_a_ind,hbs_a_val)
  DEALLOCATE(hbs_b_ind,hbs_b_val)
  DEALLOCATE(num_hbs_a,num_hbs_b)

END SUBROUTINE GET_HBONDS_ROO_ANG



SUBROUTINE GET_HBONDS_SKINNER (diff_set,pbc_opt,acc_A,acc_sH_A,acc_H_A,don_A,don_sH_A,don_H_A,coors,box,ortho, &
     acc_B,acc_sH_B,acc_H_B,don_B,don_sH_B,don_H_B,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, &
     nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H,numat_glob)

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  INTEGER,INTENT(IN)::diff_set,pbc_opt,ortho,numat_glob
  INTEGER,INTENT(IN)::nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H
  INTEGER,INTENT(IN)::nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H
  INTEGER,DIMENSION(nA_acc),INTENT(IN)    ::acc_A
  INTEGER,DIMENSION(nB_acc),INTENT(IN)    ::acc_B
  INTEGER,DIMENSION(nA_don),INTENT(IN)    ::don_A
  INTEGER,DIMENSION(nB_don),INTENT(IN)    ::don_B
  INTEGER,DIMENSION(nA_acc_sH),INTENT(IN) ::acc_sH_A
  INTEGER,DIMENSION(nB_acc_sH),INTENT(IN) ::acc_sH_B
  INTEGER,DIMENSION(nA_acc_H),INTENT(IN)  ::acc_H_A
  INTEGER,DIMENSION(nB_acc_H),INTENT(IN)  ::acc_H_B
  INTEGER,DIMENSION(nA_don_sH),INTENT(IN) ::don_sH_A
  INTEGER,DIMENSION(nB_don_sH),INTENT(IN) ::don_sH_B
  INTEGER,DIMENSION(nA_don_H),INTENT(IN)  ::don_H_A
  INTEGER,DIMENSION(nB_don_H),INTENT(IN)  ::don_H_B
  DOUBLE PRECISION,DIMENSION(numat_glob,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box


  TYPE(iarray_pointer),DIMENSION(:),POINTER::hbs_a_ind,hbs_b_ind 
  TYPE(darray_pointer),DIMENSION(:),POINTER::hbs_a_val,hbs_b_val 
  INTEGER,DIMENSION(:),ALLOCATABLE::num_hbs_a,num_hbs_b

  INTEGER::ii,jj,gg
  INTEGER::don,don_h,acc
  INTEGER::lim_hbs

  DOUBLE PRECISION,DIMENSION(3)::pos_acc,pos_don,pos_h
  DOUBLE PRECISION,DIMENSION(3)::vect_don_acc,vect_h_acc,vect_don_h
  DOUBLE PRECISION,DIMENSION(3)::aux_vect_1,aux_vect_2,aux_vect_3
  DOUBLE PRECISION::dist_h_acc,dist_don_acc,dist_don_h,aux_cos,sk_val
  DOUBLE PRECISION::dist2_h_acc,dist2_don_acc,dist2_don_h
  DOUBLE PRECISION:: cutdist2hacc,pi

  INTEGER::acc_H1,acc_H2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_perp

  INTEGER,DIMENSION(:),ALLOCATABLE::aux_box_ind,aux2_box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_box_val,aux2_box_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  pi=3.14159265358979
  cutdist2hacc=(-0.3430d0*log(sk_param/7.10d0)+0.010d0)**2 !2.30770d0 cut off distance

  lim_hbs=3

  ALLOCATE(aux_box_ind(lim_hbs),aux_box_val(lim_hbs),filtro(lim_hbs))
  filtro=.FALSE.

  ALLOCATE(hbs_a_ind(nA_don_H),hbs_a_val(nA_don_H))
  ALLOCATE(hbs_b_ind(nB_don_H),hbs_b_val(nB_don_H))
  ALLOCATE(num_hbs_a(nA_don_H),num_hbs_b(nB_don_H))

  num_hbs_a=0
  num_hbs_b=0

  !!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)

  ALLOCATE(vect_perp(nB_acc,3))
  DO ii=1,nB_acc
     acc=acc_B(ii)+1
     jj=acc_sH_B(ii)+1
     acc_H1=acc_H_B(jj)+1
     acc_H2=acc_H_B(jj+1)+1
     pos_acc=coors(acc,:)
     aux_vect_1=coors(acc_H1,:)-pos_acc(:)
     aux_vect_2=coors(acc_H2,:)-pos_acc(:)
     IF (pbc_opt==1) CALL PBC(aux_vect_1,box,ortho)
     IF (pbc_opt==1) CALL PBC(aux_vect_2,box,ortho)
     CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)
     vect_perp(ii,:)=aux_vect_3
  END DO

  DO ii=1,nA_don
     don=don_A(ii)+1
     DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
        don_h=don_H_A(hh)+1
        pos_h=coors(don_h,:)
        gg=0
        DO jj=1,nB_acc
           acc=acc_B(jj)+1
           IF (don/=acc) THEN
              vect_h_acc=pos_h(:)-coors(acc,:)
              IF (pbc_opt==1) CALL PBC(vect_h_acc,box,ortho)
              dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)
              IF (dist2_h_acc<cutdist2hacc) THEN
                 dist_h_acc=sqrt(dist2_h_acc)
                 aux_cos=dot_product(vect_perp(jj,:),vect_h_acc(:))/dist_h_acc
                 if (aux_cos>=1.0d0) aux_cos=1.0d0
                 if (aux_cos<=-1.0d0) aux_cos=-1.0d0
                 aux_cos=acos(abs(aux_cos))
                 aux_cos=aux_cos*(180.0d0/pi)
                 IF (aux_cos>90.0d0) THEN
                    print*,'aquiii error 3.14',aux_cos
                    STOP
                 END IF
                 IF (aux_cos<0.0d0) THEN
                    print*,'aquiii error 3.14',aux_cos
                    STOP
                 END IF
                 sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)
                 IF (sk_val>sk_param) THEN
                    gg=gg+1
                    IF (gg>lim_hbs) THEN
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                       aux2_box_ind=aux_box_ind
                       aux2_box_val=aux_box_val
                       DEALLOCATE(aux_box_ind,aux_box_val)
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                       aux_box_ind(1:lim_hbs)=aux2_box_ind
                       aux_box_val(1:lim_hbs)=aux2_box_val
                       DEALLOCATE(aux2_box_ind,aux2_box_val)
                       lim_hbs=gg
                    END IF
                    aux_box_ind(gg)=acc-1
                    aux_box_val(gg)=sk_val
                 END IF
              END IF
           END IF
        END DO
        IF (gg>0) THEN
           ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))
           filtro(1:gg)=.TRUE.
           DO jj=1,gg
              ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
              filtro(ll)=.FALSE.
              hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)
              hbs_a_val(hh)%d1(jj)=aux_box_val(ll)
           END DO
        END IF

        num_hbs_a(hh)=gg

     END DO
  END DO

  IF (diff_set==1) THEN

     DEALLOCATE(vect_perp)
     ALLOCATE(vect_perp(nA_acc,3))
     DO ii=1,nA_acc
        acc=acc_A(ii)+1
        jj=acc_sH_A(ii)+1
        acc_H1=acc_H_A(jj)+1
        acc_H2=acc_H_A(jj+1)+1
        pos_acc=coors(acc,:)
        aux_vect_1=coors(acc_H1,:)-pos_acc(:)
        aux_vect_2=coors(acc_H2,:)-pos_acc(:)
        IF (pbc_opt==1) CALL PBC(aux_vect_1,box,ortho)
        IF (pbc_opt==1) CALL PBC(aux_vect_2,box,ortho)
        CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)
        vect_perp(ii,:)=aux_vect_3
     END DO

     DO ii=1,nB_don
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           don_h=don_H_B(hh)+1
           pos_h=coors(don_h,:)
           gg=0
           DO jj=1,nA_acc
              acc=acc_A(jj)+1
              vect_h_acc=pos_h(:)-coors(acc,:)
              IF (pbc_opt==1) CALL PBC(vect_h_acc,box,ortho)
              dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)
              IF (dist2_h_acc<cutdist2hacc) THEN
                 dist_h_acc=sqrt(dist2_h_acc)
                 aux_cos=dot_product(vect_perp(jj,:),vect_h_acc(:))/dist_h_acc
                 if (aux_cos>=1.0d0) aux_cos=1.0d0
                 if (aux_cos<=-1.0d0) aux_cos=-1.0d0
                 aux_cos=acos(abs(aux_cos))
                 aux_cos=aux_cos*(180.0d0/pi)
                 IF (aux_cos>90.0d0) THEN
                    print*,'aquiii error 3.14',aux_cos
                    STOP
                 END IF
                 IF (aux_cos<0.0d0) THEN
                    print*,'aquiii error 3.14',aux_cos
                    STOP
                 END IF
                 sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)
                 IF (sk_val>sk_param) THEN
                    gg=gg+1
                    IF (gg>lim_hbs) THEN
                       ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                       aux2_box_ind=aux_box_ind
                       aux2_box_val=aux_box_val
                       DEALLOCATE(aux_box_ind,aux_box_val)
                       ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                       aux_box_ind(1:lim_hbs)=aux2_box_ind
                       aux_box_val(1:lim_hbs)=aux2_box_val
                       DEALLOCATE(aux2_box_ind,aux2_box_val)
                       lim_hbs=gg
                    END IF
                    aux_box_ind(gg)=acc-1
                    aux_box_val(gg)=sk_val
                 END IF
              END IF
           END DO
           IF (gg>0) THEN
              ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))
              filtro(1:gg)=.TRUE.
              DO jj=1,gg
                 ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
                 filtro(ll)=.FALSE.
                 hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)
                 hbs_b_val(hh)%d1(jj)=aux_box_val(ll)
              END DO
           END IF
           
           num_hbs_b(hh)=gg
           
        END DO
     END DO
  END IF

  DEALLOCATE(vect_perp)
  IF (ALLOCATED(hbs_out)) DEALLOCATE(hbs_out)
  IF (ALLOCATED(hbs_vals_out)) DEALLOCATE(hbs_vals_out)

  IF (diff_set==1) THEN
     ii=SUM(num_hbs_a)+SUM(num_hbs_b)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=hbs_a_val(hh)%d1(jj)
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
     DO ii=1,nB_don
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           IF (num_hbs_b(hh)>0) THEN
              DO jj=1,num_hbs_b(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_B(ii)
                 hbs_out(gg,2)=don_H_B(hh)
                 hbs_out(gg,3)=hbs_b_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=hbs_b_val(hh)%d1(jj)
              END DO
              DEALLOCATE(hbs_b_ind(hh)%i1,hbs_b_val(hh)%d1)
           END IF
        END DO
     END DO
  ELSE
     ii=SUM(num_hbs_a)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=hbs_a_val(hh)%d1(jj)
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
  END IF

  DEALLOCATE(aux_box_ind,aux_box_val,filtro)
  DEALLOCATE(hbs_a_ind,hbs_a_val)
  DEALLOCATE(hbs_b_ind,hbs_b_val)
  DEALLOCATE(num_hbs_a,num_hbs_b)

END SUBROUTINE GET_HBONDS_SKINNER


SUBROUTINE GET_HBONDS_SKINNER_NS_LIST (diff_set,pbc_opt,acc_A,acc_sH_A,acc_H_A,don_A,don_sH_A,don_H_A,coors,box,ortho, &
     acc_B,acc_sH_B,acc_H_B,don_B,don_sH_B,don_H_B,nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H, &
     nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H,numat_glob)

  TYPE iarray_pointer
     INTEGER,DIMENSION(:),POINTER::i1
  END TYPE iarray_pointer
  TYPE darray_pointer
     DOUBLE PRECISION,DIMENSION(:),POINTER::d1
  END TYPE darray_pointer

  INTEGER,INTENT(IN)::diff_set,pbc_opt,ortho,numat_glob
  INTEGER,INTENT(IN)::nA_acc,nA_acc_sH,nA_acc_H,nA_don,nA_don_sH,nA_don_H
  INTEGER,INTENT(IN)::nB_acc,nB_acc_sH,nB_acc_H,nB_don,nB_don_sH,nB_don_H
  INTEGER,DIMENSION(nA_acc),INTENT(IN)    ::acc_A
  INTEGER,DIMENSION(nB_acc),INTENT(IN)    ::acc_B
  INTEGER,DIMENSION(nA_don),INTENT(IN)    ::don_A
  INTEGER,DIMENSION(nB_don),INTENT(IN)    ::don_B
  INTEGER,DIMENSION(nA_acc_sH),INTENT(IN) ::acc_sH_A
  INTEGER,DIMENSION(nB_acc_sH),INTENT(IN) ::acc_sH_B
  INTEGER,DIMENSION(nA_acc_H),INTENT(IN)  ::acc_H_A
  INTEGER,DIMENSION(nB_acc_H),INTENT(IN)  ::acc_H_B
  INTEGER,DIMENSION(nA_don_sH),INTENT(IN) ::don_sH_A
  INTEGER,DIMENSION(nB_don_sH),INTENT(IN) ::don_sH_B
  INTEGER,DIMENSION(nA_don_H),INTENT(IN)  ::don_H_A
  INTEGER,DIMENSION(nB_don_H),INTENT(IN)  ::don_H_B
  DOUBLE PRECISION,DIMENSION(numat_glob,3),intent(in)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box


  TYPE(iarray_pointer),DIMENSION(:),POINTER::hbs_a_ind,hbs_b_ind 
  TYPE(darray_pointer),DIMENSION(:),POINTER::hbs_a_val,hbs_b_val 
  INTEGER,DIMENSION(:),ALLOCATABLE::num_hbs_a,num_hbs_b

  INTEGER::ii,jj,gg
  INTEGER::don,don_h,acc
  INTEGER::lim_hbs

  DOUBLE PRECISION,DIMENSION(3)::pos_acc,pos_don,pos_h
  DOUBLE PRECISION,DIMENSION(3)::vect_don_acc,vect_h_acc,vect_don_h
  DOUBLE PRECISION,DIMENSION(3)::aux_vect_1,aux_vect_2,aux_vect_3
  DOUBLE PRECISION::dist_h_acc,dist_don_acc,dist_don_h,aux_cos,sk_val
  DOUBLE PRECISION::dist2_h_acc,dist2_don_acc,dist2_don_h
  DOUBLE PRECISION:: cutdist2hacc,pi

  INTEGER::acc_H1,acc_H2
  DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE::vect_perp
  INTEGER,DIMENSION(:),ALLOCATABLE::inds_perp

  INTEGER,DIMENSION(:),ALLOCATABLE::aux_box_ind,aux2_box_ind
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::aux_box_val,aux2_box_val

  LOGICAL,DIMENSION(:),ALLOCATABLE::filtro

  pi=3.14159265358979
  cutdist2hacc=(-0.3430d0*log(sk_param/7.10d0)+0.010d0)**2 !2.30770d0 cut off distance

  lim_hbs=3

  ALLOCATE(aux_box_ind(lim_hbs),aux_box_val(lim_hbs),filtro(lim_hbs))
  filtro=.FALSE.

  ALLOCATE(hbs_a_ind(nA_don_H),hbs_a_val(nA_don_H))
  ALLOCATE(hbs_b_ind(nB_don_H),hbs_b_val(nB_don_H))
  ALLOCATE(num_hbs_a(nA_don_H),num_hbs_b(nB_don_H))

  num_hbs_a=0
  num_hbs_b=0

  !!!! Source: A. Luzar, D. Chandler. Phys. Rev. Lett. 76, 928-931 (1996)


  CALL EXTRACT_NS_LIST_SETS(diff_set,don_H_A,acc_B,nA_don_H,nB_acc,numat_glob)

  ALLOCATE(vect_perp(nB_acc,3),inds_perp(numat_glob))
  DO ii=1,nB_acc
     acc=acc_B(ii)+1
     inds_perp(acc)=ii
     jj=acc_sH_B(ii)+1
     acc_H1=acc_H_B(jj)+1
     acc_H2=acc_H_B(jj+1)+1
     pos_acc=coors(acc,:)
     aux_vect_1=coors(acc_H1,:)-pos_acc(:)
     aux_vect_2=coors(acc_H2,:)-pos_acc(:)
     IF (pbc_opt==1) CALL PBC(aux_vect_1,box,ortho)
     IF (pbc_opt==1) CALL PBC(aux_vect_2,box,ortho)
     CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)
     vect_perp(ii,:)=aux_vect_3
  END DO


  DO ii=1,nA_don
     don=don_A(ii)+1
     DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
        don_h=don_H_A(hh)+1
        pos_h=coors(don_h,:)
        gg=0
        DO jj=1,ver_ic_dim(don_h)
           IF (filt_sets_ns_ind(don_h,jj)) THEN
              acc=ver_ic_ind(don_h,jj)
              IF (don/=acc) THEN
                 vect_h_acc=pos_h(:)-coors(acc,:)
                 IF (pbc_opt==1) CALL PBC(vect_h_acc,box,ortho)
                 dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)
                 IF (dist2_h_acc<cutdist2hacc) THEN
                    dist_h_acc=sqrt(dist2_h_acc)
                    aux_cos=dot_product(vect_perp(inds_perp(acc),:),vect_h_acc(:))/dist_h_acc
                    if (aux_cos>=1.0d0) aux_cos=1.0d0
                    if (aux_cos<=-1.0d0) aux_cos=-1.0d0
                    aux_cos=acos(abs(aux_cos))
                    aux_cos=aux_cos*(180.0d0/pi)
                    IF (aux_cos>90.0d0) THEN
                       print*,'aquiii error 3.14',aux_cos
                       STOP
                    END IF
                    IF (aux_cos<0.0d0) THEN
                       print*,'aquiii error 3.14',aux_cos
                       STOP
                    END IF
                    sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)
                    IF (sk_val>sk_param) THEN
                       gg=gg+1
                       IF (gg>lim_hbs) THEN
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                          aux2_box_ind=aux_box_ind
                          aux2_box_val=aux_box_val
                          DEALLOCATE(aux_box_ind,aux_box_val)
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                          aux_box_ind(1:lim_hbs)=aux2_box_ind
                          aux_box_val(1:lim_hbs)=aux2_box_val
                          DEALLOCATE(aux2_box_ind,aux2_box_val)
                          lim_hbs=gg
                       END IF
                       aux_box_ind(gg)=acc-1
                       aux_box_val(gg)=sk_val
                    END IF
                 END IF
              END IF
           END IF
        END DO
        IF (gg>0) THEN
           ALLOCATE(hbs_a_ind(hh)%i1(gg),hbs_a_val(hh)%d1(gg))
           filtro(1:gg)=.TRUE.
           DO jj=1,gg
              ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
              filtro(ll)=.FALSE.
              hbs_a_ind(hh)%i1(jj)=aux_box_ind(ll)
              hbs_a_val(hh)%d1(jj)=aux_box_val(ll)
           END DO
        END IF

        num_hbs_a(hh)=gg

     END DO
  END DO

  IF (diff_set==1) THEN

     CALL EXTRACT_NS_LIST_SETS(diff_set,don_H_B,acc_A,nB_don_H,nA_acc,numat_glob)

     DEALLOCATE(vect_perp,inds_perp)
     ALLOCATE(vect_perp(nA_acc,3),inds_perp(numat_glob))
     DO ii=1,nA_acc
        acc=acc_A(ii)+1
        inds_perp(acc)=ii
        jj=acc_sH_A(ii)+1
        acc_H1=acc_H_A(jj)+1
        acc_H2=acc_H_A(jj+1)+1
        pos_acc=coors(acc,:)
        aux_vect_1=coors(acc_H1,:)-pos_acc(:)
        aux_vect_2=coors(acc_H2,:)-pos_acc(:)
        IF (pbc_opt==1) CALL PBC(aux_vect_1,box,ortho)
        IF (pbc_opt==1) CALL PBC(aux_vect_2,box,ortho)
        CALL PERPENDICULAR_NORMED_VECT (aux_vect_1,aux_vect_2,aux_vect_3)
        vect_perp(ii,:)=aux_vect_3
     END DO

     DO ii=1,nB_don
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           don_h=don_H_B(hh)+1
           pos_h=coors(don_h,:)
           gg=0
           DO jj=1,ver_ic_dim(don)
              IF (filt_sets_ns_ind(don,jj)) THEN
                 acc=ver_ic_ind(don,jj)
                 vect_h_acc=pos_h(:)-coors(acc,:)
                 IF (pbc_opt==1) CALL PBC(vect_h_acc,box,ortho)
                 dist2_h_acc=dot_product(vect_h_acc,vect_h_acc)
                 IF (dist2_h_acc<cutdist2hacc) THEN
                    dist_h_acc=sqrt(dist2_h_acc)
                    aux_cos=dot_product(vect_perp(inds_perp(acc),:),vect_h_acc(:))/dist_h_acc
                    if (aux_cos>=1.0d0) aux_cos=1.0d0
                    if (aux_cos<=-1.0d0) aux_cos=-1.0d0
                    aux_cos=acos(abs(aux_cos))
                    aux_cos=aux_cos*(180.0d0/pi)
                    IF (aux_cos>90.0d0) THEN
                       print*,'aquiii error 3.14',aux_cos
                       STOP
                    END IF
                    IF (aux_cos<0.0d0) THEN
                       print*,'aquiii error 3.14',aux_cos
                       STOP
                    END IF
                    sk_val=exp(-dist_h_acc/0.3430d0)*(7.10d0-0.050d0*aux_cos+0.000210d0*aux_cos**2)
                    IF (sk_val>sk_param) THEN
                       gg=gg+1
                       IF (gg>lim_hbs) THEN
                          ALLOCATE(aux2_box_ind(lim_hbs),aux2_box_val(lim_hbs))
                          aux2_box_ind=aux_box_ind
                          aux2_box_val=aux_box_val
                          DEALLOCATE(aux_box_ind,aux_box_val)
                          ALLOCATE(aux_box_ind(gg),aux_box_val(gg))
                          aux_box_ind(1:lim_hbs)=aux2_box_ind
                          aux_box_val(1:lim_hbs)=aux2_box_val
                          DEALLOCATE(aux2_box_ind,aux2_box_val)
                          lim_hbs=gg
                       END IF
                       aux_box_ind(gg)=acc-1
                       aux_box_val(gg)=sk_val
                    END IF
                 END IF
              END IF
           END DO
           IF (gg>0) THEN
              ALLOCATE(hbs_b_ind(hh)%i1(gg),hbs_b_val(hh)%d1(gg))
              filtro(1:gg)=.TRUE.
              DO jj=1,gg
                 ll=MAXLOC(aux_box_val(:),DIM=1,MASK=filtro)
                 filtro(ll)=.FALSE.
                 hbs_b_ind(hh)%i1(jj)=aux_box_ind(ll)
                 hbs_b_val(hh)%d1(jj)=aux_box_val(ll)
              END DO
           END IF
           
           num_hbs_b(hh)=gg
           
        END DO
     END DO
  END IF

  IF (ALLOCATED(hbs_out)) DEALLOCATE(hbs_out)
  IF (ALLOCATED(hbs_vals_out)) DEALLOCATE(hbs_vals_out)


  IF (diff_set==1) THEN
     ii=SUM(num_hbs_a)+SUM(num_hbs_b)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=hbs_a_val(hh)%d1(jj)
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
     DO ii=1,nB_don
        DO hh=don_sH_B(ii)+1,don_sH_B(ii+1)
           IF (num_hbs_b(hh)>0) THEN
              DO jj=1,num_hbs_b(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_B(ii)
                 hbs_out(gg,2)=don_H_B(hh)
                 hbs_out(gg,3)=hbs_b_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=hbs_b_val(hh)%d1(jj)
              END DO
              DEALLOCATE(hbs_b_ind(hh)%i1,hbs_b_val(hh)%d1)
           END IF
        END DO
     END DO
  ELSE
     ii=SUM(num_hbs_a)
     ALLOCATE(hbs_out(ii,3),hbs_vals_out(ii))
     gg=0
     DO ii=1,nA_don
        DO hh=don_sH_A(ii)+1,don_sH_A(ii+1)
           IF (num_hbs_a(hh)>0) THEN
              DO jj=1,num_hbs_a(hh)
                 gg=gg+1
                 hbs_out(gg,1)=don_A(ii)
                 hbs_out(gg,2)=don_H_A(hh)
                 hbs_out(gg,3)=hbs_a_ind(hh)%i1(jj)
                 hbs_vals_out(gg)=hbs_a_val(hh)%d1(jj)
              END DO
              DEALLOCATE(hbs_a_ind(hh)%i1,hbs_a_val(hh)%d1)
           END IF
        END DO
     END DO
  END IF

  DEALLOCATE(vect_perp,inds_perp)
  DEALLOCATE(aux_box_ind,aux_box_val,filtro)
  DEALLOCATE(hbs_a_ind,hbs_a_val)
  DEALLOCATE(hbs_b_ind,hbs_b_val)
  DEALLOCATE(num_hbs_a,num_hbs_b)

END SUBROUTINE GET_HBONDS_SKINNER_NS_LIST

END MODULE HBONDS
