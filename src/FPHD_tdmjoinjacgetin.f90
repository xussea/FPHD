module tdmjoinjacgetin
        implicit none
contains
subroutine Get_TDM(Nbasis,MO,S,AM_X,AM_Y,mo_o,mo_v,TD,P)
implicit none
integer :: i,j,k,a,b,t
integer, intent(in) :: Nbasis,mo_o,mo_v
real*8 :: PS
real*8, dimension(:,:), intent(out) :: P
real*8, dimension(:,:), intent(in) :: S,MO,AM_X,AM_Y
character (len=5) :: TD

P(:,:)=0.00d00
if ( TD == "CIS") then

do i=1,mo_o
        do j=mo_o+1,mo_o+mo_v
                P(i,j)=AM_X(i,j)
                P(j,i)=AM_X(i,j)
        end do
end do

else if (TD == "STDA") then

do i=1,mo_o
        do j=mo_o+1,mo_o+mo_v
                P(i,j)=AM_X(i,j)
                P(j,i)=AM_X(i,j)
        end do
end do


else

do i=1,mo_o
        do j=mo_o+1,mo_o+mo_v
                P(i,j)=AM_X(i,j)+AM_Y(i,j)
                P(j,i)=AM_X(i,j)+AM_Y(i,j)
        end do
end do

end if

!P=matmul(matmul(MO,P),transpose(MO))
!do i=1,Nbasis
!       do j=1,Nbasis
!                PS=(P(i,j)+P(j,i))/2.0
!                P(i,j)=PS
!                P(j,i)=PS
!        end do
! end do


return
end subroutine Get_TDM

subroutine Get_PET_PHT(Nbasis,mo_o,mo_v,MO,P,ROET,ROHT,tre,trh)
implicit none
real(8), dimension(:,:), intent(in) :: P, MO
real(8), dimension(:,:), allocatable :: PPT
real(8), dimension(:,:), intent(out) :: ROET, ROHT
integer :: i, j
real(8) :: tre, trh
integer, intent(in) :: Nbasis, mo_o, mo_v
ROET=0.0d00
ROHT=0.0d00
trh=0.0
tre=0.0

allocate(PPT(Nbasis,Nbasis))
PPT=0.0d00
PPT=matmul(P,transpose(P))

do i=1, mo_o
        do j=1, mo_o
                ROHT(i,j)=PPT(i,j)
        enddo
enddo

do i=mo_o+1,mo_o+mo_v
        do j=mo_o+1,mo_o+mo_v
                ROET(i,j)=PPT(i,j)
        enddo
enddo
deallocate(PPT)

do i=1, mo_o
        trh=trh+ROHT(i,i)
enddo

do i=mo_o+1, mo_o+mo_v
        tre=tre+ROET(i,i)
enddo

ROHT=matmul(MO,matmul(ROHT,transpose(MO)))
ROET=matmul(MO,matmul(ROET,transpose(MO)))

return
end subroutine

subroutine JoinJac(A,V,D,thr)
implicit none
real*8, dimension(:,:), intent(in) :: A
real*8, dimension(:,:), intent(out) :: V,D
real*8, intent(in) :: thr
integer :: i,j,k,p,q,m,nm
real*8 :: encore,ton,toff,theta,c,s,tton,ttoff
real*8 :: itol
real*8, dimension(:,:), allocatable :: g,gt,mp,mq
real*8, dimension(:,:), allocatable :: tmp
real*8, dimension(:), allocatable :: rowp,rowq,Vq,Vp
real*8, parameter :: pi=3.14159265359

m=size(A,1)
nm=size(A,2)

V(:,:)=0.0
do i=1,m
        V(i,i)=1.0
end do

allocate(g(2,nm/m))
allocate(gt(2,2))
allocate(mp(m,nm/m))
allocate(mq(m,nm/m))
allocate(tmp(m,nm))
allocate(rowp(m))
allocate(rowq(m))
allocate(Vp(m))
allocate(Vq(m))
encore=1.0

tmp=A
do while (encore > 0.0 )
        itol=0.0
        do p=1,m-1
        do q=p+1,m
                g(1,:)=tmp(p,p:nm:m)-tmp(q,q:nm:m)
                g(2,:)=tmp(p,q:nm:m)+tmp(q,p:nm:m)
                gt=matmul(g,transpose(g))
                ton=gt(1,1)-gt(2,2)
                toff=gt(1,2)+gt(2,1)

                theta=0.5*atan2(toff,ton+sqrt(ton*ton+toff*toff))
!               print*, theta
                c=cos(theta)
                s=sin(theta)

                itol=itol+abs(s)
                if ( abs(s) > thr ) then
                        Mp(:,:)=tmp(:,p:nm:m)
                        Mq(:,:)=tmp(:,q:nm:m)
                        tmp(:,p:nm:m)=c*Mp+s*Mq
                        tmp(:,q:nm:m)=c*Mq-s*Mp
                        rowp=tmp(p,:)
                        rowq=tmp(q,:)
                        tmp(p,:)=c*rowp+s*rowq;
                        tmp(q,:)=c*rowq-s*rowp;
                        Vp=V(:,p)
                        Vq=V(:,q)
                        V(:,p)=c*Vp+s*Vq
                        V(:,q)=c*Vq-s*Vp
                end if
        end do
        end do
        itol=itol/(5*m**2)
!print*, itol,thr
        if ( itol < thr ) then
                encore=-1.0
        end if
end do

D=tmp
!print*, D
return
end subroutine JoinJac

subroutine getin(k,N,i,j)
implicit none
integer, intent(in) :: k,N
integer, intent(out) :: i,j
integer :: z

buc: do i=1,N+1
z=(i*(i-1))/2+i
if (z >= k ) then
        exit buc
end if
end do buc

j=k-(i*(i-1))/2

return
end subroutine getin
end module tdmjoinjacgetin

