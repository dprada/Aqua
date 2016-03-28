!
SUBROUTINE NORMALIZE_VECT (a)

  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::a
  DOUBLE PRECISION::norm

  norm=sqrt(dot_product(a,a))
  a=a/norm

END SUBROUTINE NORMALIZE_VECT

SUBROUTINE PRODUCT_VECT(a,b,normal)

  DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::a,b
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT)::normal
  
  normal(1)=a(2)*b(3)-a(3)*b(2)
  normal(2)=-a(1)*b(3)+a(3)*b(1)
  normal(3)=a(1)*b(2)-a(2)*b(1)

END SUBROUTINE PRODUCT_VECT


SUBROUTINE PBC(vector,box,inv,ortho)
 
  IMPLICIT NONE
 
  DOUBLE PRECISION,DIMENSION(3),INTENT(INOUT)::vector
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box,inv
  DOUBLE PRECISION,DIMENSION(3)::vaux,vaux2
  INTEGER,INTENT(IN)::ortho
  INTEGER::i,j,k
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
 
     IF (.TRUE.) THEN
        !vaux(1)=inv(1,1)*vector(1)
        !vaux(2)=inv(2,1)*vector(1)+inv(2,2)*vector(2)
        !vaux(3)=inv(3,1)*vector(1)+inv(3,2)*vector(2)+inv(3,3)*vector(3)
        vaux(1)=inv(1,1)*vector(1)+inv(2,1)*vector(2)+inv(3,1)*vector(3)
        vaux(2)=                   inv(2,2)*vector(2)+inv(3,2)*vector(3)
        vaux(3)=                                      inv(3,3)*vector(3)
        vaux(1)=vaux(1)-NINT(vaux(1))*1.0
        vaux(2)=vaux(2)-NINT(vaux(2))*1.0
        vaux(3)=vaux(3)-NINT(vaux(3))*1.0
        vector(1)=box(1,1)*vaux(1)+box(2,1)*vaux(2)+box(3,1)*vaux(3)
        vector(2)=                 box(2,2)*vaux(2)+box(3,2)*vaux(3)
        vector(3)=                                  box(3,3)*vaux(3)
     ELSE
        L=1000000.0d0
        vaux2=0.0d0
        DO i=-1,1,1
           DO j=-1,1,1
              DO k=-1,1,1
                 vaux(:)=vector(:)+i*box(1,:)+j*box(2,:)+k*box(3,:)
                 x=(vaux(1)*vaux(1)+vaux(2)*vaux(2)+vaux(3)*vaux(3))
                 IF (L>x) THEN
                    vaux2=vaux
                    L=x
                 END IF
              END DO
           END DO
        END DO
        vector(:)=vaux2(:)
     END IF

  END IF
  
END SUBROUTINE PBC

SUBROUTINE DISTANCE_POINTS (points1,points2,n1,n2,pbc_opt,box,inv,ortho,matrix)

INTEGER,INTENT(IN)::pbc_opt,ortho
INTEGER,INTENT(IN)::n1,n2
INTEGER,DIMENSION(n1,3),INTENT(IN)::points1
INTEGER,DIMENSION(n2,3),INTENT(IN)::points2
DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box,inv
DOUBLE PRECISION,dimension(n1,n2),intent(out)::matrix

INTEGER::ii,jj
DOUBLE PRECISION,DIMENSION(3)::vect_aux,vect_aux2

matrix=0.0d0

DO ii=1,n1

   vect_aux=points1(ii,:)

   DO jj=1,n2

      vect_aux2(:)=points2(jj,:)-vect_aux(:)
      CALL PBC(vect_aux2,box,inv,ortho)
      matrix(ii,jj)=sqrt(dot_product(vect_aux2,vect_aux2))

   END DO

END DO

END SUBROUTINE DISTANCE_POINTS

SUBROUTINE DISTANCE_TWO_POINTS (point1,point2,pbc_opt,box,inv,ortho,dd)
 
INTEGER,INTENT(IN)::pbc_opt,ortho
DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::point1
DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::point2
DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box,inv
DOUBLE PRECISION,DIMENSION(3)::vect_aux
DOUBLE PRECISION,INTENT(OUT)::dd
 
vect_aux(:)=point2(:)-point1(:)
CALL PBC(vect_aux,box,inv,ortho)
dd=sqrt(dot_product(vect_aux,vect_aux))
 
END SUBROUTINE DISTANCE_TWO_POINTS


SUBROUTINE ANGLE_THREE_POINTS_PROJXY (point1,point2,point3,pbc_opt,box,inv,ortho,ang)

INTEGER,INTENT(IN)::pbc_opt,ortho
DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::point1,point2,point3
DOUBLE PRECISION,DIMENSION(3,3),INTENT(IN)::box,inv
DOUBLE PRECISION,DIMENSION(3)::vect12,vect32,cross_prod
DOUBLE PRECISION::dot_prod
DOUBLE PRECISION,INTENT(OUT)::ang


vect12(:) = point1(:) - point2(:)
vect32(:) = point3(:) - point2(:)

CALL PBC(vect12,box,inv,ortho)
CALL PBC(vect32,box,inv,ortho)

vect12(3) = 0.0d0
vect32(3) = 0.0d0

CALL NORMALIZE_VECT (vect12)
CALL NORMALIZE_VECT (vect32)
CALL PRODUCT_VECT (vect12,vect32,cross_prod)
dot_prod=dot_product(vect12,vect32)

IF (abs(dot_prod)>1.0d0) THEN
   IF (abs(dot_prod)<1.0010d0) THEN
      dot_prod=SIGN(1.0d0,dot_prod)
   END IF
END IF

ang=ACOS(dot_prod)
ang=SIGN(ang,cross_prod(3))


END SUBROUTINE ANGLE_THREE_POINTS_PROJXY

SUBROUTINE ELEGIDOS (distances,helices,num_heads,num_ats_hx,eleg_val,eleg_ind)

INTEGER,INTENT(IN)::num_heads,num_ats_hx
DOUBLE PRECISION,DIMENSION(num_heads,num_ats_hx),INTENT(IN)::distances
INTEGER,DIMENSION(num_ats_hx,3),INTENT(IN)::helices
DOUBLE PRECISION,DIMENSION(num_heads,128),INTENT(OUT)::eleg_val
INTEGER,DIMENSION(num_heads,128),INTENT(OUT)::eleg_ind

INTEGER::ii,jj,destino

eleg_val(:,:)=1000.0

DO ii=1,num_heads
   DO jj=1,num_ats_hx
      destino=helices(jj,2)*8+helices(jj,3)+1
      IF (eleg_val(ii,destino)>distances(ii,jj)) THEN
         eleg_val(ii,destino)=distances(ii,jj)
         eleg_ind(ii,destino)=jj-1
      END IF
   END DO
END DO

END SUBROUTINE ELEGIDOS
