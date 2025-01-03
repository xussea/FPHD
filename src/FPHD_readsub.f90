module readdiab
      implicit none
contains
subroutine DR(fname,NS,Nbasis,NA,Nelec,Z,MO,mo_o,mo_v,AM_X,AM_Y,AE,TD,S,dip,DPE,DPM,DVPE)
!subroutine DR(fname,NS,Nbasis,NA,Nelec,Z,MO,mo_o,mo_v,AM_X,AM_Y,AE,TD,S)
implicit none
character(len=200)::fname, line
character(len=5)::TD,Spin
character(len=1)::dip
integer,intent(out):: NA, Nelec, Nbasis, mo_o, mo_v
integer::ic,i,j,NX,NY,NS,a,b,st
integer,dimension(:),allocatable,intent(out)::Z
real(8),dimension(:),allocatable::ZZ
real(8),dimension(:,:),allocatable,intent(out)::MO
real(8),dimension(:,:),allocatable,intent(out)::S
real(8),dimension(:),allocatable,intent(out)::AE
real(8),dimension(:,:,:),allocatable,intent(out)::AM_X, AM_Y
real(8), dimension(:,:), allocatable,intent(out) :: DPE, DPM, DVPE

open(unit=1,file=trim(fname)//'.diab.bin',form='unformatted',access='stream', action='read',status='old')

read(1)NA
!write(*,*)NA
read(1)Nelec
!write(*,*)Nelec
read(1)Nbasis
!write(*,*)Nbasis

mo_o=int(Nelec/2)+mod(Nelec,2)
mo_v=Nbasis-mo_o

allocate(Z(NA))

read(1) Z(1:NA)
!write(*,*) Z
!allocate(ZZ(NA))
!117   read(999,'(A)')line
!ic=index(line,'Atomic Number:')
!if(ic.eq.0)go to 117
!read(999,*) ZZ(1:NA)
!Z=int(ZZ)

allocate(MO(Nbasis,Nbasis))

MO(:,:)=0
do i=1, Nbasis
        read(1) (MO(i,j),j=1,Nbasis)
!        write(*,*) (MO(i,j),j=1,Nbasis)
end do

!do i=1, Nbasis
!        write(*,*) MO(i,1:Nbasis)
!end do


!MO(:,:)=0.
!118   read(999,'(A)') line
!ic=index(line,'Molecular Orbitals:')
!if(ic.eq.0)go to 118
!do i=1, Nbasis
!        read(999,*) (MO(i,j),j=1,Nbasis)
!end do

allocate(S(Nbasis,Nbasis))

S(:,:)=0
do i=1, Nbasis
        read(1) (S(i,j),j=1, Nbasis)
        !write(*,*) (S(i,j),j=1,Nbasis)
end do
!S(:,:)=0
!119     read(999,'(A)') line
!ic=index(line,'Overlap Matrix:')
!if(ic.eq.0)go to 119
!do i=1,Nbasis
!        write(*,*) (S(i,j),j=1,Nbasis)
!end do

allocate(AE(NS))
allocate(AM_X(Nbasis,Nbasis,NS))
allocate(AM_Y(Nbasis,Nbasis,NS))
AM_X=0.00d00
AM_Y=0.00d00
do i=1, NS
 read(1)AE(i)
 !write(*,*) AE(i)
end do
!do i=1,NS
!write(*,700) AE(1:i)
!enddo
!stop

do i=1, NS
 if(TD == "CIS") then
         read(1) NX
         do j=1,NX
                read(1) a, b, AM_X(a,b,i)
                !write(*,*) a, b, AM_X(a,b,i)
         end do
 else if(TD == "STDA") then
         read(1) NX
         do j=1,NX
                read(1) a, b, AM_X(a,b,i)
         end do
 else if(TD == "TD") then
         read(1)NX
         read(1)NY
         do j=1,NX
                read(1) a, b, AM_X(a,b,i)
         end do
         do j=1,NY
                read(1) a, b, AM_Y(a,b,i)
         end do
 end if
end do

if (dip == "Y") then
 allocate(DPE(NS,3),DPM(NS,3),DVPE(NS,3))

 do i=1,NS
   read(1) DPE(i,1), DPE(i,2), DPE(i,3)
   read(1) DPM(i,1), DPM(i,2), DPM(i,3)
   read(1) DVPE(i,1), DVPE(i,2), DVPE(i,3)
 end do
end if

!write(*,*) AE

!700 format(5(E13.5,2X))
!stop
close(1)
end subroutine
subroutine unDR(fname,NS,Nbasis,NA,Neleca,Nelecb,Z,MOa,MOb,mo_oa,mo_va,mo_ob,mo_vb,AM_Xa,AM_Ya,AM_Xb,AM_Yb,AE,TD,S,dip,DPE,DPM,DVPE)
!subroutine unDR(fname,NS,Nbasis,NA,Neleca,Nelecb,Z,MOa,MOb,mo_oa,mo_va,mo_ob,mo_vb,AM_Xa,AM_Ya,AM_Xb,AM_Yb,AE,TD,S)
implicit none
character(len=200)::fname, line
character(len=5)::TD,Spin
character(len=1):: dip
integer,intent(out):: NA, Nbasis
integer, intent(out):: Neleca, Nelecb, mo_oa, mo_va, mo_ob, mo_vb
integer::ic,i,j,NXa,NYa,NS,a,b,st,NXb, NYb
integer,dimension(:),allocatable,intent(out)::Z
real(8),dimension(:),allocatable::ZZ
real(8),dimension(:,:),allocatable,intent(out)::MOa, MOb
real(8),dimension(:,:),allocatable,intent(out)::S
real(8),dimension(:),allocatable,intent(out)::AE
real(8),dimension(:,:,:),allocatable,intent(out)::AM_Xa, AM_Ya, AM_Xb, AM_Yb
real(8), dimension(:,:), allocatable,intent(out) :: DPE, DPM, DVPE

open(unit=1,file=trim(fname)//'.diab.bin',form='unformatted',access='stream', action='read',status='old')

read(1)NA
!write(*,*)NA
read(1)Neleca,Nelecb
!write(*,*)Nelec
read(1)Nbasis
!write(*,*)Nbasis
!mo_oa=int(Neleca/2)+mod(Neleca,2)
!mo_ob=int(Nelecb/2)+mod(Nelecb,2)
mo_oa=Neleca
mo_ob=Nelecb
mo_va=Nbasis-mo_oa
mo_vb=Nbasis-mo_ob
!!!!!!
!write(*,*)'MO_OA', mo_oa, 'MO_OB', mo_ob
!write(*,*)'MO_VA', mo_va, 'MO_VB', mo_vb
!stop
!!!!!!
allocate(Z(NA))

read(1) Z(1:NA)
!write(*,*) Z
!allocate(ZZ(NA))
!117   read(999,'(A)')line
!ic=index(line,'Atomic Number:')
!if(ic.eq.0)go to 117
!read(999,*) ZZ(1:NA)
!Z=int(ZZ)

allocate(MOa(Nbasis,Nbasis),MOb(Nbasis,Nbasis))

MOa(:,:)=0
MOb(:,:)=0
do i=1, Nbasis
        read(1) (MOa(i,j),j=1,Nbasis)
!        write(*,*) (MO(i,j),j=1,Nbasis)
enddo
do i =1, Nbasis
        read(1) (MOb(i,j),j=1,Nbasis)
enddo
write(62,*)"MOa:"
do i=1, Nbasis
        write(62,*)(MOa(i,j),j=1,Nbasis)
enddo
write(62,*)"MOb:"
do i=1, Nbasis
        write(62,*)(MOb(i,j),j=1,Nbasis)
enddo

!do i=1, Nbasis
!        write(*,*) MO(i,1:Nbasis)
!end do


!MO(:,:)=0.
!118   read(999,'(A)') line
!ic=index(line,'Molecular Orbitals:')
!if(ic.eq.0)go to 118
!do i=1, Nbasis
!        read(999,*) (MO(i,j),j=1,Nbasis)
!end do

allocate(S(Nbasis,Nbasis))

S(:,:)=0
do i=1, Nbasis
        read(1) (S(i,j),j=1, Nbasis)
        !write(*,*) (S(i,j),j=1,Nbasis)
end do
!S(:,:)=0
!119     read(999,'(A)') line
!ic=index(line,'Overlap Matrix:')
!if(ic.eq.0)go to 119
!do i=1,Nbasis
!        write(*,*) (S(i,j),j=1,Nbasis)
!end do

allocate(AE(NS))
allocate(AM_Xa(Nbasis,Nbasis,NS))
allocate(AM_Ya(Nbasis,Nbasis,NS))
allocate(AM_Xb(Nbasis,Nbasis,NS))
allocate(AM_Yb(Nbasis,Nbasis,NS))
AM_Xa=0.00d00
AM_Ya=0.00d00
AM_Xb=0.00d00
AM_Yb=0.00d00
do i=1, NS
 read(1)AE(i)
end do
write(34,700)AE(1:i)
!stop

do i=1, NS
 if(TD == "CIS") then
         read(1) NXa, NXb
         write(34,*)NXa,NXb
         do j=1,NXa
                read(1) a, b, AM_Xa(a,b,i)
                write(34,*) a, b, AM_Xa(a,b,i)
         end do
         do j=1,NXb
                read(1) a, b, AM_Xb(a,b,i)
                write(34,*) a, b, AM_Xb(a,b,i)
         end do
 else if(TD == "TD") then
         read(1)NXa,NXb
         read(1)NYa,NYb
         !write(*,*)'XA STATE',i
         do j=1,NXa
                read(1) a, b, AM_Xa(a,b,i)
                !write(*,*)a,b,AM_Xa(a,b,i)
         end do
         !write(*,*)'XB STATE',i
         do j=1,NXb
                read(1) a, b, AM_Xb(a,b,i)
                !write(*,*)a,b,AM_Xb(a,b,i)
         end do
         !write(*,*)'YA STATE',i
         do j=1,NYa
                read(1) a, b, AM_Ya(a,b,i)
                !write(*,*)a,b,AM_Ya(a,b,i)
         end do
         !write(*,*)'YB STATE',i
         do j=1,NYb
                read(1) a, b, AM_Yb(a,b,i)
                !write(*,*)a,b,AM_Yb(a,b,i)
         end do
 end if
end do

if (dip == "Y") then
 allocate(DPE(NS,3),DPM(NS,3),DVPE(NS,3))

 do i=1,NS
   read(1) DPE(i,1), DPE(i,2), DPE(i,3)
   read(1) DPM(i,1), DPM(i,2), DPM(i,3)
   read(1) DVPE(i,1), DVPE(i,2), DVPE(i,3)
 end do
end if
!do i =1, NS
!  write(*,*) AM_Xa(:,:,i)
!enddo

!stop
!write(*,*) AE

700 format(5(E13.5,2X))
!stop
close(1)
end subroutine
subroutine SOCDR(nm,NS,Nbasis,NA,Nelca,Nelcb,Z,MOa,MOb,mo_oa,mo_va,mo_ob,mo_vb,AM_Xa,AM_Ya,AM_Xb,AM_Yb,AE,TD,S,dip,DPE,DPM,DVPE,SOC)
implicit none
character(len=200)::nm, line
character(len=5)::TD,Spin
character(len=1):: dip
integer,intent(out):: NA, Nbasis
integer, intent(out):: Nelca, Nelcb, mo_oa, mo_va, mo_ob, mo_vb
integer::ic,i,j,NXa,NYa,NS,a,b,st,NXb, NYb
integer,dimension(:),allocatable,intent(out)::Z
real(8),dimension(:),allocatable::ZZ
real(8),dimension(:,:),allocatable,intent(out)::MOa, MOb
real(8),dimension(:,:),allocatable,intent(out)::S, SOC
real(8),dimension(:),allocatable,intent(out)::AE
real(8),dimension(:,:,:),allocatable,intent(out)::AM_Xa, AM_Ya, AM_Xb, AM_Yb
real(8), dimension(:,:), allocatable,intent(out) :: DPE, DPM, DVPE

open(unit=1,file=trim(nm)//'.diab.bin',form='unformatted',access='stream', action='read',status='old')

read(1)NA
!write(*,*)NA
read(1)Nelca,Nelcb
!write(*,*)Nelec
read(1)Nbasis
mo_oa=Nelca
mo_ob=Nelcb
mo_va=Nbasis-mo_oa
mo_vb=Nbasis-mo_ob
!!!!!!
!write(*,*)'MO_OA', mo_oa, 'MO_OB', mo_ob
!write(*,*)'MO_VA', mo_va, 'MO_VB', mo_vb
!stop
!!!!!!
allocate(Z(NA))

read(1) Z(1:NA)
!write(*,*) Z
!allocate(ZZ(NA))
!117   read(999,'(A)')line
!ic=index(line,'Atomic Number:')
!if(ic.eq.0)go to 117
!read(999,*) ZZ(1:NA)
!Z=int(ZZ)

allocate(MOa(Nbasis,Nbasis),MOb(Nbasis,Nbasis))

MOa(:,:)=0
MOb(:,:)=0
do i=1, Nbasis
        read(1) (MOa(i,j),j=1,Nbasis)
enddo
do i =1, Nbasis
        read(1) (MOb(i,j),j=1,Nbasis)
enddo
!do i=1, Nbasis
!        write(62,*)(MOa(i,j),j=1,Nbasis)
!enddo
!do i=1, Nbasis
!        write(62,*)(MOb(i,j),j=1,Nbasis)
!enddo

allocate(S(Nbasis,Nbasis))

S(:,:)=0
do i=1, Nbasis
        read(1) (S(i,j),j=1, Nbasis)
        !write(*,*) (S(i,j),j=1,Nbasis)
end do
!S(:,:)=0
!119     read(999,'(A)') line
!ic=index(line,'Overlap Matrix:')
!if(ic.eq.0)go to 119
!do i=1,Nbasis
!        write(*,*) (S(i,j),j=1,Nbasis)
!end do

allocate(AE(NS))
allocate(AM_Xa(Nbasis,Nbasis,NS))
allocate(AM_Ya(Nbasis,Nbasis,NS))
allocate(AM_Xb(Nbasis,Nbasis,NS))
allocate(AM_Yb(Nbasis,Nbasis,NS))
AM_Xa=0.00d00
AM_Ya=0.00d00
AM_Xb=0.00d00
AM_Yb=0.00d00
do i=1, NS
 read(1)AE(i)
end do
write(34,700)AE(1:i)

allocate(SOC(NS,NS))
SOC = 0.00d00
do i=1, NS
 read(1) (SOC(i,j), j=1,NS)
enddo

do i=1, NS
 if(TD == "CIS") then
         read(1) NXa, NXb
         write(34,*)NXa,NXb
         do j=1,NXa
                read(1) a, b, AM_Xa(a,b,i)
                write(34,*) a, b, AM_Xa(a,b,i)
         end do
         do j=1,NXb
                read(1) a, b, AM_Xb(a,b,i)
                write(34,*) a, b, AM_Xb(a,b,i)
         end do
 else if(TD == "TD") then
         read(1)NXa,NXb
         read(1)NYa,NYb
         write(34,*)'XA STATE',i
         do j=1,NXa
                read(1) a, b, AM_Xa(a,b,i)
                write(34,*)a,b,AM_Xa(a,b,i)
         end do
         write(34,*)'XB STATE',i
         do j=1,NXb
                read(1) a, b, AM_Xb(a,b,i)
                write(34,*)a,b,AM_Xb(a,b,i)
         end do
         write(34,*)'YA STATE',i
         do j=1,NYa
                read(1) a, b, AM_Ya(a,b,i)
                write(34,*)a,b,AM_Ya(a,b,i)
         end do
         write(34,*)'YB STATE',i
         do j=1,NYb
                read(1) a, b, AM_Yb(a,b,i)
                write(34,*)a,b,AM_Yb(a,b,i)
         end do
 end if
enddo

if (dip == "Y") then
 allocate(DPE(NS,3),DPM(NS,3),DVPE(NS,3))

 do i=1,NS
   read(1) DPE(i,1), DPE(i,2), DPE(i,3)
   read(1) DPM(i,1), DPM(i,2), DPM(i,3)
   read(1) DVPE(i,1), DVPE(i,2), DVPE(i,3)
 end do
end if

700 format(5(E13.5,2X))
!stop
close(1)
end subroutine
end module readdiab

