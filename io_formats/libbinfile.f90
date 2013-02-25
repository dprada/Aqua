SUBROUTINE open_read(len_ch,file_name,funit,o_natom,o_cell,pos_o)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::len_ch
  CHARACTER(80),INTENT(IN)::file_name
  INTEGER,INTENT(OUT)::funit,o_natom
  INTEGER(KIND=8),INTENT(OUT)::pos_o
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::o_cell

  LOGICAL:: UNITOP
  REAL::cell(9)

  UNITOP=.False.
  funit=500

  do while (UNITOP)
     inquire (unit=funit,opened=UNITOP)
     if (UNITOP) funit=funit+1
  end do
  
  if (len_ch>80) then
     PRINT*, '# Error: Name of file too long.'
  end if

  OPEN (unit=funit,file=TRIM(file_name),status='old',action='read', &
       form='unformatted',access='stream')

  READ(funit) o_natom,cell

  o_cell(1,:)=10.0d0*dble(cell(1:3))
  o_cell(2,:)=10.0d0*dble(cell(4:6))
  o_cell(3,:)=10.0d0*dble(cell(7:9))
  o_cell(1,2)=90.0d0
  o_cell(1,3)=90.0d0
  o_cell(2,3)=90.0d0

  INQUIRE(funit,pos=pos_o)

END SUBROUTINE open_read

SUBROUTINE read (funit,natom,pos_i,pos_o,step,time,prec,cell,coors,io_err,io_end)

  IMPLICIT NONE
  INTEGER,INTENT(IN)::funit,natom
  INTEGER(KIND=8),INTENT(IN)::pos_i
  INTEGER,INTENT(OUT)::io_err,io_end
  INTEGER(KIND=8),INTENT(OUT)::pos_o
  DOUBLE PRECISION,DIMENSION(natom,3),INTENT(OUT)::coors
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT)::cell
  INTEGER,INTENT(OUT)::step
  DOUBLE PRECISION,INTENT(OUT)::prec,time

  REAL,DIMENSION(:),ALLOCATABLE::box_buffer
  REAL,DIMENSION(:),ALLOCATABLE::x_buffer
  REAL::prec_buffer,time_buffer
  INTEGER::num_atoms,ii,jj,gg

  pos_o=pos_i
  io_err=0
  io_end=1

  ALLOCATE(box_buffer(9),x_buffer(natom*3))


  READ(funit,pos=pos_o,end=600) num_atoms,step,time_buffer,box_buffer,x_buffer,prec_buffer
  time=dble(time_buffer)
  prec=dble(prec_buffer)
  
  cell(1,:)=10.0d0*dble(box_buffer(1:3))
  cell(2,:)=10.0d0*dble(box_buffer(4:6))
  cell(3,:)=10.0d0*dble(box_buffer(7:9))
  cell(1,2)=90.0d0
  cell(1,3)=90.0d0
  cell(2,3)=90.0d0

  gg=0
  DO ii=1,num_atoms
     DO jj=1,3
        gg=gg+1
        coors(ii,jj)=10.0d0*dble(x_buffer(gg))
     END DO
  END DO

  io_end=0

  INQUIRE(funit,pos=pos_o)

600 DEALLOCATE(box_buffer,x_buffer)

END SUBROUTINE read

SUBROUTINE close(funit,io_err)

  IMPLICIT NONE
INTEGER,INTENT(IN)::funit
  INTEGER,INTENT(OUT)::io_err

  CLOSE(funit)
  io_err=0 !good

END SUBROUTINE close
