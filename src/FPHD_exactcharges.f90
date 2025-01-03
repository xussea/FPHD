module exactcharges
      implicit none
contains
subroutine Ex_ATC(NA,Nfrag,NAfrag,Nbasis,MO,AM_X1,AM_Y1,AM_X2,AM_Y2,mo_o,mo_v,AI,AB,AN,Z,S,Vh,Ve,TD)
implicit none
integer :: i,j,k,a,b,t
integer, intent(in) :: NA,Nbasis,mo_o,mo_v,Nfrag,AI
integer, dimension(:), intent(in) :: AB,AN,Z,NAfrag
Real*8, Dimension (:,:), intent(in):: AM_X1,AM_Y1,AM_X2,AM_Y2,S
real*8, allocatable, dimension(:,:) :: Ph,Pe,P1,P2
Real*8, Dimension (:,:),intent(inout) :: MO
character (len=5) :: TD
real*8, dimension(:), intent(out) :: Vh,Ve

allocate(Ph(mo_o,mo_o))        ! CIS holes density matrix
allocate(Pe(mo_v,mo_v))        ! CIS electron density matrix
allocate(P1(Nbasis,Nbasis))
allocate(P2(Nbasis,Nbasis))

Ph(:,:)=0.
Pe(:,:)=0.
P1(:,:)=0.
P2(:,:)=0.
Vh(:)=0.0
Ve(:)=0.0

if ( TD == "CIS" ) then
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X1(i,mo_o+a)*AM_X2(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X1(i,mo_o+a)*AM_X2(i,mo_o+b)
                end do
        end do
end do

else if ( TD == "STDA" ) then
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X1(i,mo_o+a)*AM_X2(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X1(i,mo_o+a)*AM_X2(i,mo_o+b)
                end do
        end do
end do


else

do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                       !Ph(i,j)=Ph(i,j)+(AM_X1(i,mo_o+a)+AM_Y1(i,mo_o+a))*(AM_X2(j,mo_o+a)+AM_Y2(j,mo_o+a))
                       Ph(i,j)=Ph(i,j)+(AM_X1(i,mo_o+a)*AM_X2(j,mo_o+a))+(AM_Y2(j,mo_o+a)*AM_Y1(i,mo_o+a))
                end do
        end do
end do


do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        !Pe(a,b)=Pe(a,b)+(AM_X1(i,mo_o+a)+AM_Y1(i,mo_o+a))*(AM_X2(i,mo_o+b)+AM_Y2(i,mo_o+b))
                        Pe(a,b)=Pe(a,b)+(AM_X1(i,mo_o+a)*AM_X2(i,mo_o+b))+(AM_Y1(i,mo_o+a)*AM_Y2(i,mo_o+b))
                end do
        end do
end do
end if

do i=1,mo_o
        do j=i+1,mo_o
                Ph(i,j)=(Ph(i,j)+Ph(j,i))*0.5
                !Ph(i,j)=(Ph(i,j)+Ph(j,i))
                Ph(j,i)=Ph(i,j)
        end do
end do

do i=1,mo_v
        do j=1,mo_v
                Pe(i,j)=(Pe(i,j)+Pe(j,i))*0.5
                !Pe(i,j)=(Pe(i,j)+Pe(j,i))
                Pe(j,i)=Pe(i,j)
        end do
end do


P1=matmul(MO(:,1:mo_o),matmul(Ph,transpose(MO(:,1:mo_o))))
P2=matmul(MO(:,mo_o+1:mo_o+mo_v),matmul(Pe,transpose(MO(:,mo_o+1:mo_o+mo_v))))
!P1=matmul(P1,S)
!P2=matmul(P2,S)
!P1=matmul(S,P1)
!P2=matmul(S,P2)

k=0
t=0
do i=1,Nfrag
        do j=1,NAfrag(i)
                t=t+1
NFB:            do a=1,AI
                if ( Z(t) == AN(a) ) then
                        exit NFB
                end if
                end do NFB

                do b=1,AB(a)
                        k=k+1
                        Vh(i)=Vh(i)+P1(k,k)
                        Ve(i)=Ve(i)+P2(k,k)
                end do
        end do
end do

deallocate(Ph)        ! CIS holes density matrix
deallocate(Pe)        ! CIS electron density matrix
deallocate(P1)        ! CIS holes density matrix
deallocate(P2)        ! CIS electron density matrix

return
end subroutine Ex_ATC


subroutine Ex_charges(NA,Nfrag,NAfrag,Nbasis,MO,AM_X,AM_Y,mo_o,mo_v,AI,AB,AN,Z,S,Vh,Ve,TD)
implicit none
integer :: i,j,k,a,b,t
integer, intent(in) :: NA,Nbasis,mo_o,mo_v,Nfrag,AI
integer, dimension(:), intent(in) :: AB,AN,Z,NAfrag
Real*8, Dimension (:,:), intent(in):: AM_X,AM_Y,S
real*8, allocatable, dimension(:,:) :: Ph,Pe,P1,P2
Real*8, Dimension (:,:),intent(inout) :: MO
character (len=3) :: TD
real*8, dimension(:), intent(out) :: Vh,Ve

allocate(Ph(mo_o,mo_o))        ! CIS holes density matrix
allocate(Pe(mo_v,mo_v))        ! CIS electron density matrix
allocate(P1(Nbasis,Nbasis))
allocate(P2(Nbasis,Nbasis))

Ph(:,:)=0.
Pe(:,:)=0.
P1(:,:)=0.
P2(:,:)=0.
Vh(:)=0.0
Ve(:)=0.0

if (TD=="CIS") then

!Holes density matrix
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X(i,mo_o+a)*AM_X(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X(i,mo_o+a)*AM_X(i,mo_o+b)
                end do
        end do
end do

else if ( TD == "STDA") then
do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+AM_X(i,mo_o+a)*AM_X(j,mo_o+a)
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+AM_X(i,mo_o+a)*AM_X(i,mo_o+b)
                end do
        end do
end do
        
else

do i=1,mo_o
        do j=1,mo_o
                do a=1,mo_v
                        Ph(i,j)=Ph(i,j)+(AM_X(i,mo_o+a)+AM_Y(i,mo_o+a))*(AM_X(j,mo_o+a)-AM_Y(j,mo_o+a))
                end do
        end do
end do

do a=1,mo_v
        do b=1,mo_v
                do i=1,mo_o
                        Pe(a,b)=Pe(a,b)+(AM_X(i,mo_o+a)+AM_Y(i,mo_o+a))*(AM_X(i,mo_o+b)-AM_Y(i,mo_o+b))
                end do
        end do
end do
end if

P1=matmul(MO(:,1:mo_o),matmul(Ph,transpose(MO(:,1:mo_o))))
P2=matmul(MO(:,mo_o+1:mo_o+mo_v),matmul(Pe,transpose(MO(:,mo_o+1:mo_o+mo_v))))
!P1=matmul(P1,S)
!P2=matmul(P2,S)
!P1=matmul(S,P1)
!P2=matmul(S,P2)

k=0
t=0
do i=1,Nfrag
        do j=1,NAfrag(i)
                t=t+1
NFB2:           do a=1,AI
                if ( Z(t) == AN(a) ) then
                        exit NFB2
                end if
                end do NFB2

                do b=1,AB(a)
                        k=k+1
                        Vh(i)=Vh(i)+P1(k,k)
                        Ve(i)=Ve(i)+P2(k,k)
                end do
        end do
end do

deallocate(Ph)        ! CIS holes density matrix
deallocate(Pe)        ! CIS electron density matrix
deallocate(P1)        ! CIS holes density matrix
deallocate(P2)        ! CIS electron density matrix
return
end subroutine Ex_charges
end module exactcharges
