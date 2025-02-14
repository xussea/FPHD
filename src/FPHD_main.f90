program main
use readdiab
use stateanalysis
use exactcharges
use tdmjoinjacgetin
implicit none
integer :: IOS
integer :: i,j,k,a,b,AI
integer :: NA,Nbasis,Nelec,Neleca,Nelecb
integer :: mo_o,mo_v,mo_oa,mo_va,mo_ob,mo_vb
integer :: NStates,NS,Nfrag,kstates
integer, dimension(:), allocatable :: NAFrag,AB,AN,Z,state_vec,alpbetsp
integer, dimension(:,:), allocatable :: state_mat
real*8, dimension(:), allocatable :: AE,hh,ee,auxh,auxe, hha,eea,hhb,eeb
real*8, dimension(:,:), allocatable :: S,S_tmp,S2,MO,MOS,JJ,JV,JD,H
real*8, dimension(:,:), allocatable :: MOa,MOb,MOSa,MOSb
real*8, dimension(:,:), allocatable :: trot,rot,rot_tmp,HD_tmp,H_tmp
real*8, dimension(:,:,:), allocatable :: AM_X,AM_Y,Dh,De,TDMA,TDMD,TDMAa,TDMAb,TDMDa,TDMDb
real*8,dimension(:,:,:), allocatable :: AM_Xa,AM_Ya,AM_Xb,AM_Yb,Dha,Dhb,Dea,Deb
real*8,dimension(:,:,:), allocatable :: Pet,Pht,Peta,Petb,Phta,Phtb
real*8 :: Dthr,tol, tre, trh, trea, treb, trha, trhb
character (len=200) :: Inp
character (len=5) :: TD
character (len=1) :: Spin,dip,DenMat,soc

!real*8, dimension(:,:), allocatable :: dipo

real*8, dimension(:), allocatable :: work,eigen
integer :: lwork, info

character(len=100):: pythonex, cpFPHDpy, rmFPHDpy
integer :: exits

integer :: nargs
character(len=100) :: arg

real(8),dimension(:,:),allocatable:: DPE, DPM, DVPE
real(8), dimension(:,:), allocatable:: SOC_M

!Input example
! J. Phys. Chem. Lett. 2021, 12, 3, 1032–1039
! $DIAB
! Inp = prefix of the calculation ; state vectors must have the name prefix.state"i" and adiabatic transition dipole moments in prefix.dipo
! Nstates = number of states to diabatize must be equal or lower than the total excited states computed in the QM calculation. It is recomended ask 10 or 20 states more than the ones of interest
! TD = TD or CIS
! Nfrag = Number of fragments
! AI = number of different type of atoms
! tol = threshold for jacobi diagonalization
! Dthr = threshold between 0 and 1 for clasify the states
! dip = Y or N computes the diabatic transition dipole moment, requires .dipo file
! DenMat = Y or N writes binary file of transition density matrix
! $END
!
! $FRAG
! NAFrag = array of number of atoms for each fragment, in the QM calculation the atoms should be ordered by fragments
! AN = array of the atomic number of each atom
! AB = array with the number of basis functions of each type of atom, the order must be the same as the atomic numbers in AN
! $END
!
!
! Compile as:  gfortran -O2 -fopenmp -fdefault-real-8 main.f90 -llapack -lblas -L PATH-TO-LAPACK-BLAS -o EXEC_NAME
! Excute as:  ./progran < input > output
namelist /DIAB/ Inp,Spin,Nstates,TD,Nfrag,AI,tol,Dthr,DenMat,dip,soc
namelist /FRAG/ NAFrag,AN,AB

tol=1.0e-08
TD='CIS'
Nfrag=2
Nstates=2
DenMat='N'
dip='N'
Spin='R'
SOC='N'
nargs = command_argument_count()
if (nargs > 0) then
call get_command_argument(1,arg)
else
       write(*,*)"WARNING: You don't put a input file"
end if

open(unit=10, file=arg, status='old', action='read')
read(10,NML=DIAB)
allocate(NAFrag(Nfrag))
allocate(AB(AI))
allocate(AN(AI))
read(10,NML=FRAG)

if ( Dthr == 0.0 ) then
 Dthr=0.8
end if

!cpFPHDpy = "cp /home/pau/.soft/FPHD/FPHD_Unrestricted/FPHD_file_reader_xtb.py ." 
!call execute_command_line(cpFPHDpy, wait=.TRUE., exitstat=exits)
!        if(exits /= 0) then
!                print *, "Error al copiar el script de Python"
!        end if

!pythonex="python3 FPHD_file_reader.py <" // trim(Inp) // "_diab.inp"

pythonex="python3 $FPHDPATH/src/FPHD_file_reader_cis.py <" // trim(arg)
call execute_command_line(pythonex, wait=.TRUE., exitstat=exits)
    if (exits /= 0) then
        print*, "Warning: There's something wrong in the input file of this program"
        stop
    end if

!rmFPHDpy = "rm FPHD_file_reader_xtb.py"
!call execute_command_line(rmFPHDpy, wait=.TRUE., exitstat=exits)
!        if(exits /= 0) then
!                print *, "Error al eliminar el script de Python"
!        end if

!call fchkreader(Inp,NA,Z,Nbasis,Nelec,MO,mo_o,mo_v)
!call overlap_reader(Inp,Nbasis,S)

!if ( NAFrag(Nfrag) < 0 ) then
!       NAFrag(:)=NA/Nfrag
!end if
!AM_X=0.0d00
!AM_Y=0.0d00
!

!call readatmelebasis(Inp,Nbasis,NA,Nelec)

!allocate(Z(NA))
!allocate(MO(Nbasis,Nbasis), S(Nbasis,Nbasis))
!allocate(AE(Nstates))
!allocate(AM_X(Nbasis,Nbasis,Nstates))
!allocate(AM_Y(Nbasis,Nbasis,Nstates))

!call readrest(Inp,Nstates,Nbasis,NA,Nelec,Z,MO,mo_o,mo_v,AM_X,AM_Y,AE,TD,S)
!stop
!write(*,*)"LIST"
if (soc == "Y") then
call SOCDR(Inp,Nstates,Nbasis,NA,Neleca,Nelecb,Z,MOa,MOb,mo_oa,mo_va,mo_ob,mo_vb, &
 & AM_Xa,AM_Ya,AM_Xb,AM_Yb,AE,TD,S,dip,DPE,DPM,DVPE,SOC_M)
else        
 if (Spin == "U") then
  call unDR(Inp,Nstates,Nbasis,NA,Neleca,Nelecb,Z,MOa,MOb,mo_oa,mo_va,mo_ob,mo_vb,AM_Xa,AM_Ya,AM_Xb,AM_Yb,AE,TD,S,dip,DPE,DPM,DVPE)
 else
  call DR(Inp,Nstates,Nbasis,NA,Nelec,Z,MO,mo_o,mo_v,AM_X,AM_Y,AE,TD,S,dip,DPE,DPM,DVPE)
 end if
end if
!allocate(SOC(Nstates,Nstates))
!write(*,*)'SOC:'
!do i=1, Nstates
!        write(*,'(20(E13.5,2x))') SOC_M(i,1:Nstates)
!end do
!write(*,*)'MOa'
!do i=1, Nbasis
!        write(*,'(752(E13.5,2x))') (MOa(i,j),j=1,Nbasis)
!enddo
!write(*,*)'MOb'
!do i=1, Nbasis
!        write(*,'(752(E13.5,2x))') (MOb(i,j),j=1,Nbasis)
!enddo
 !AM_X=AM_X*dsqrt(2.0d0)
 !AM_y=AM_y*:wdsqrt(2.0d0)

!AM_X=0.0d00
!AM_Y=0.0d00
!allocate(AE(Nstates))
!allocate(AM_X(Nbasis,Nbasis,Nstates))
!allocate(AM_Y(Nbasis,Nbasis,Nstates))
!allocate(TDMA(Nbasis,Nbasis,Nstates))
!allocate(TDMD(Nbasis,Nbasis,Nstates))
!allocate(auxe((Nbasis*(Nbasis-1))/2+Nbasis))
!allocate(auxh((Nbasis*(Nbasis-1))/2+Nbasis))

!do NS=1,Nstates
!	call state_reader(Inp,NS,Nbasis,mo_o,mo_v,AM_X(:,:,NS),AM_Y(:,:,NS),AE(NS),TD,'Y')
!end do

!print*, Nbasis
!stop

 allocate(S2(Nbasis,Nbasis))
 allocate(eigen(Nbasis))
 S2(:,:)=0.00
 
    allocate(work(1))
    lwork=-1
    call DSYEV('V', 'L', Nbasis, S, size(S,1), eigen, work, lwork, info)
    if (info /= 0) then
       stop 'dsyev 1 failed'
    endif
    lwork = int(work(1))
    deallocate(work)
    
    allocate(work(lwork))
    call DSYEV('V', 'U', Nbasis, S, size(S,1), eigen, work, lwork, info)
    if (info /= 0) then
       stop 'dsyev 2 failed'
    endif
    deallocate(work)
 
 do i=1,Nbasis
 S2(i,i)=sqrt(eigen(i))
 end do
 
 S2=matmul(matmul(S,S2),transpose(S))
! if (Spin == "U") then
if (Spin == "U" .or. soc == "Y") then
         allocate(MOSa(Nbasis,Nbasis))
         allocate(MOSb(Nbasis,Nbasis))
         MOSa=matmul(S2,MOa)
         MOSb=matmul(S2,MOb)
 else
         allocate(MOS(Nbasis,Nbasis))
         MOS=matmul(S2,MO)
 end if

 
!if (Spin == "U") then
if (Spin == "U" .or. soc == "Y") then
 allocate(Dh(Nstates,Nstates,Nfrag))
 allocate(De(Nstates,Nstates,Nfrag))
 allocate(Dha(Nstates,Nstates,Nfrag))
 allocate(Dhb(Nstates,Nstates,Nfrag))
 allocate(Dea(Nstates,Nstates,Nfrag))
 allocate(Deb(Nstates,Nstates,Nfrag))
 Dh = 0.0
 De = 0.0
 Dha = 0.0
 Dhb = 0.0
 Dea = 0.0
 Deb = 0.0
 !write(34,*)"MO_Oa",mo_oa, "MO_Va", mo_va
 !write(34,*)"MO_Ob",mo_ob, "MO_Vb", mo_vb
 !Fill diagonal of matrices

 !$OMP PARALLEL DEFAULT(SHARED) private(i,j,k)
 !$OMP DO
 do k=1,Nstates*(Nstates-1)/2+Nstates

 call getin(k,Nstates,i,j)

 if ( i == j ) then
 call Ex_charges(NA,Nfrag,NAfrag,Nbasis,MOSa,AM_Xa(:,:,i),AM_Ya(:,:,i), &
 & mo_oa,mo_va,AI,AB,AN,Z,S,Dha(i,i,:),Dea(i,i,:),TD)
 call Ex_charges(NA,Nfrag,NAfrag,Nbasis,MOSb,AM_Xb(:,:,i),AM_Yb(:,:,i), &
 & mo_ob,mo_vb,AI,AB,AN,Z,S,Dhb(i,i,:),Deb(i,i,:),TD)
 else
 call Ex_ATC(NA,Nfrag,NAfrag,Nbasis,MOSa,AM_Xa(:,:,i),AM_Ya(:,:,i),AM_Xa(:,:,j), &
 & AM_Ya(:,:,j),mo_oa,mo_va,AI,AB,AN,Z,S,Dha(i,j,:),Dea(i,j,:),TD)
 Dha(j,i,:)=Dha(i,j,:)
 Dea(j,i,:)=Dea(i,j,:)
 call Ex_ATC(NA,Nfrag,NAfrag,Nbasis,MOSb,AM_Xb(:,:,i),AM_Yb(:,:,i),AM_Xb(:,:,j), &
 & AM_Yb(:,:,j),mo_ob,mo_vb,AI,AB,AN,Z,S,Dhb(i,j,:),Deb(i,j,:),TD)
 Dhb(j,i,:)=Dhb(i,j,:)
 Deb(j,i,:)=Deb(i,j,:)
 end if

 end do
 !$OMP END DO
 !$OMP BARRIER
 !$OMP END PARALLEL

!Dh = Dha + Dhb
!De = Dea + Deb
Dh = (Dha + Dhb)
De = (Dea + Deb)
!Dh = Dha
!De = Dea

else

 allocate(Dh(Nstates,Nstates,Nfrag)) 
 allocate(De(Nstates,Nstates,Nfrag))
 Dh = 0.0
 De = 0.0
!Fill diagonal of matrices

 !$OMP PARALLEL DEFAULT(SHARED) private(i,j,k)
 !$OMP DO
 do k=1,Nstates*(Nstates-1)/2+Nstates

 call getin(k,Nstates,i,j)

 if ( i == j ) then
 call Ex_charges(NA,Nfrag,NAfrag,Nbasis,MOS,AM_X(:,:,i),AM_Y(:,:,i), &
 & mo_o,mo_v,AI,AB,AN,Z,S,Dh(i,i,:),De(i,i,:),TD)
 else
 call Ex_ATC(NA,Nfrag,NAfrag,Nbasis,MOS,AM_X(:,:,i),AM_Y(:,:,i),AM_X(:,:,j), &
 & AM_Y(:,:,j),mo_o,mo_v,AI,AB,AN,Z,S,Dh(i,j,:),De(i,j,:),TD)	
 Dh(j,i,:)=Dh(i,j,:)
 De(j,i,:)=De(i,j,:)
 end if

 end do
 !$OMP END DO
 !$OMP BARRIER
 !$OMP END PARALLEL
end if
!!!!!!
if (Spin == "U" .or. soc == "Y") then
 do i=1,Nfrag
! 	print*, "Dha",i
! do j=1,Nstates
! 	write(*,'(5(2X,F9.5))') Dha(j,1:j,i)
!end do
! 	print*,
!  	print*, "Dhb",i
! do j=1,Nstates
! 	write(*,'(5(2X,F9.5))') Dhb(j,1:j,i)
! end do
! 	print*,
! 	print*, "Dea",i
! do j=1,Nstates
! 	write(*,'(5(2X,F9.5))') Dea(j,1:j,i)
! end do
! 	print*,
! 	print*, "Deb",i
! do j=1,Nstates
! 	write(*,'(5(2X,F9.5))') Deb(j,1:j,i)
! end do
! 	print*,
!        write(11,*) "Dh",i
! do j=1,Nstates
!        write(11,'(5(2X,F9.5))') Dh(j,1:j,i)
! end do
!        write(11,*)""
!        write(11,*) "De",i
! do j=1,Nstates
!        write(11,'(5(2X,F9.5))') De(j,1:j,i)
! end do
!        write(11,*)""
! end do
!else
! do i=1,Nfrag
!          write(12,*) "Dh",i
! do j=1,Nstates
!        write(12,'(5(2X,F9.5))') Dh(j,1:j,i)
! end do
!        write(12,*)""
!        write(12,*) "De",i
! do j=1,Nstates
!        write(12,'(5(2X,F9.5))') De(j,1:j,i)
! end do
!        write(12,*)""
!      print*, "Dh",i
! do j=1,Nstates
!        write(*,'(5(2X,F9.5))') Dh(j,1:j,i)
!end do
!        print*,
!        print*, "De",i
! do j=1,Nstates
!        write(*,'(5(2X,F9.5))') De(j,1:j,i)
! end do
!        print*,
enddo
end if
!stop
!!!!!!
!if (Spin == "U") then
! allocate(JJ(Nstates,4*Nfrag*Nstates))
! allocate(JD(Nstates,4*Nfrag*Nstates))
! allocate(JV(Nstates,Nstates))
! allocate(H(Nstates,Nstates))
! do i=1,Nfrag
!        JJ(:,(i-1)*Nstates+1:i*Nstates)=Dha(:,:,i)
!        JJ(:,Nfrag*Nstates+(i-1)*Nstates+1:Nfrag*Nstates+i*Nstates)=Dhb(:,:,i)
!        JJ(:,2*Nfrag*Nstates+(i-1)*Nstates+1:2*Nfrag*Nstates+i*Nstates)=Dea(:,:,i)
!        JJ(:,3*Nfrag*Nstates+(i-1)*Nstates+1:3*Nfrag*Nstates+i*Nstates)=Deb(:,:,i)
! enddo
!else
 allocate(JJ(Nstates,2*Nfrag*Nstates))
 allocate(JD(Nstates,2*Nfrag*Nstates))
 allocate(JV(Nstates,Nstates))
 allocate(H(Nstates,Nstates))
 do i=1,Nfrag
        JJ(:,(i-1)*Nstates+1:i*Nstates)=Dh(:,:,i)
        JJ(:,Nfrag*Nstates+(i-1)*Nstates+1:Nfrag*Nstates+i*Nstates)=De(:,:,i)
 end do
!end if
!print*, "start JoinJac"
call JoinJac(JJ,JV,JD,tol) 
!print*, "End JoinJac"
!!!!!!
!do i=1,Nstates
!	print*, "Dh1"
!	write(*,'(5(2X,F9.5))') JD(i,1:Nstates)
!	print*, 
!end do
!do i=1,Nstates
!	print*, "Dh2"
!	write(*,'(5(2X,F9.5))') JD(i,Nstates+1:2*Nstates)
!	print*, 
!end do
!do i=1,Nstates
!	print*, "De1"
!	write(*,'(5(2X,F9.5))') JD(i,2*Nstates+1:3*Nstates)
!	print*, 
!end do
!do i=1,Nstates
!	print*, "De2"
!	write(*,'(5(2X,F9.5))') JD(i,3*Nstates+1:4*Nstates)
!	print*, 
!end do
!stop
!!!!!!
H(:,:)=0.00
do i=1,Nstates
        H(i,i)=AE(i)
end do
!print*, "BIEN"
!print*,"SOC sin rotar"
!do i=1, Nstates
!        write(*,'(10(E13.5,2x))') SOC_M(i,1:Nstates)
!end do
H=matmul(matmul(transpose(JV),H),(JV))
if (soc == "Y") then
        SOC_M = matmul(matmul(transpose(JV),SOC_M),(JV))
end if

!print*,"SOC rotada"
!do i=1, Nstates
!        write(*,'(10(E13.5,2x))') SOC_M(i,1:Nstates)
!end do
!print*, "H primera"
!do i=1, Nstates
!        write(*,'(10(E13.5,2x))') H(i,1:Nstates)
!end do
!print*, "H+SOC"
!do i=1, Nstates
!        write(*,'(10(E13.5,2x))') H(i,1:Nstates)
!end do
!print*, "Requetebien"
!stop

!!!!!
! allocate(hha(NStates))
! allocate(eea(NStates))
! allocate(hhb(NStates))
! allocate(eeb(NStates))
! print*,
! do i=1,Nfrag
!        do j=1,Nstates
!                hha(j)=Dha(j,j,i)
!                eea(j)=Dea(j,j,i)
!                hhb(j)=Dhb(j,j,i)
!                eeb(j)=Deb(j,j,i)
!        end do
!         write(90,*) "Hole population of Fragment with Alpha Spin", i
!        write(90,400) hha(1:Nstates)
!        write(90,*) " "
!        write(90,*) "Hole population of Fragment with Beta Spin", i
!        write(90,400) hhb(1:Nstates)
!        write(90,*) " "
!         write(90,*) "Particle population of Fragment with Alpha Spin", i
!        write(90,400) eea(1:NStates)
!        write(90,*) " "
!        write(90,*) "Particle population of Fragment with Beta Spin", i
!        write(90,400) eeb(1:NStates)
!        write(90,*) " "
! end do
! deallocate(hha,eea,hhb,eeb)
!stop
!!!!!

!print*,
!print*, "Diabatic Hamiltonian before subspace diagonalization"
!do i=1,Nstates
!	write(*,200) H(i,1:i)
!end do
!print*,
if (Spin == "U" .or. soc == "Y") then
 do j=1,NFrag
        Dha(:,:,j)=matmul(matmul(transpose(JV),Dha(:,:,j)),JV)
        Dhb(:,:,j)=matmul(matmul(transpose(JV),Dhb(:,:,j)),JV)
        Dea(:,:,j)=matmul(matmul(transpose(JV),Dea(:,:,j)),JV)
        Deb(:,:,j)=matmul(matmul(transpose(JV),Deb(:,:,j)),JV)
 end do
 !Dh=Dha+Dhb
 !De=Dea+Deb
Dh = (Dha + Dhb)
De = (Dea + Deb)
!Dh = Dha
!De = Dea
else
 do j=1,NFrag
        Dh(:,:,j)=matmul(matmul(transpose(JV),Dh(:,:,j)),JV)
        De(:,:,j)=matmul(matmul(transpose(JV),De(:,:,j)),JV)
 enddo
end if
!!!!!!
!If (Spin == "U" .or. soc == "Y") then
! do i=1,Nfrag
!         write(13,*) "Dh",i
! do j=1,Nstates
!        write(13,'(5(2X,F9.5))') Dh(j,1:j,i)
! end do
!        write(13,*)""
!        write(13,*) "De",i
! do j=1,Nstates
!        write(13,'(5(2X,F9.5))') De(j,1:j,i)
! end do
!        write(13,*)""
! end do
!Else
! do i=1,Nfrag
!          write(14,*) "Dh",i
! do j=1,Nstates
!        write(14,'(5(2X,F9.5))') Dh(j,1:j,i)
! end do
!        write(14,*)""
!        write(14,*) "De",i
! do j=1,Nstates
!        write(14,'(5(2X,F9.5))') De(j,1:j,i)
! end do
!        write(14,*)""
! enddo
!End if

!        write(11,*)"Dh",i 
!do j=1,Nstates
!       write(11,'(5(2X,F9.5))') Dh(j,1:j,i)
!       write(11,*)""
!end do
!        write(11,*)"De",i
       !print*, "De",i
!do j=1,Nstates
!       write(11,'(5(2X,F9.5))') De(j,1:j,i)
!       !print*,
!       write(11,*)""
!end do
!end do
!stop
!!!!!!


kstates=NFrag+NFrag*(NFrag-1)+1
allocate(state_vec(kstates))
allocate(state_mat(kstates,Nstates))


call analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dh,De,Dthr)
!if (Spin == "U") then
! allocate(alpbetsp(Nstates))
 !print*, "ALPHA SPIN"
 !call analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dha,Dea,Dthr)
 !print*, "BETHA SPIN"
 !call analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dhb,Deb,Dthr)
!Print*, "OKKK1" ,Dha(2,2,1),Dhb(2,2,1),dea(2,2,1),deb(2,2,1)
! call analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dha,Dea,Dhb,Deb,Dthr)
!Print*, "OKKK3" ,Dha(2,2,1),Dhb(2,2,1),dea(2,2,1),deb(2,2,1)
!else
 !call analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dh,De,Dthr)
 !end if
!!!!!!!!!!!!!!!!!!!!!!!
!stop
!!!!!!!!!!!!!!!!!!!!!!!
allocate(trot(Nstates,Nstates))
allocate(rot(Nstates,Nstates))

trot(:,:)=0.0
rot(:,:)=0.0
do i=1,Nstates
	trot(i,i)=1.0
	rot(i,i)=1.0
end do

do i=1,kstates

	rot(:,:)=0.0
	do j=1,Nstates
		rot(j,j)=1.0
	end do

	if (state_vec(i) /= 0 ) then

	allocate(H_tmp(state_vec(i),state_vec(i)))
	allocate(HD_tmp(state_vec(i),state_vec(i)))
	allocate(rot_tmp(state_vec(i),state_vec(i)))

	do j=1,state_vec(i)
	do k=1,state_vec(i)
		H_tmp(j,k)=H(state_mat(i,j),state_mat(i,k))
	end do
	end do

	call JoinJac(H_tmp,rot_tmp,HD_tmp,1.0E-08)

	H_tmp=matmul(transpose(rot_tmp),matmul(H_tmp,rot_tmp))

	do j=1,state_vec(i)
	do k=1,state_vec(i)
		rot(state_mat(i,j),state_mat(i,k))=rot_tmp(j,k)
	end do
	end do
	trot=matmul(trot,rot)

	print*, "subspace",i,"out of",kstates,"diagonalized"

	deallocate(H_tmp)
	deallocate(HD_tmp)
	deallocate(rot_tmp)
	end if
end do

H=matmul(transpose(trot),matmul(H,trot))
!write(57,*)"Hamiltoniano rotado sin SOC"
!do i =1, Nstates
!        write(57,'(E13.5,2X)') (H(i,j),j=1,Nstates)
!end do
if (soc == "Y") then
        SOC_M = matmul(transpose(trot),matmul(SOC_M,trot))
        !write(57,*)"SOC rotada:"
        !do i =1, Nstates
        !        write(57,'(E13.5,2X)') (SOC_M(i,j),j=1,Nstates)
        !end do
        !H=H+SOC_M
        !write(57,*)"Hamiltoniano rotado con SOC rotada"
        !do i =1, Nstates
        !        write(57,'(E13.5,2X)') (H(i,j),j=1,Nstates)
        !end do
end if

if (Spin == "U" .or. soc == "Y") then
 do j=1,NFrag
        Dha(:,:,j)=matmul(matmul(transpose(trot),Dha(:,:,j)),trot)
        Dhb(:,:,j)=matmul(matmul(transpose(trot),Dhb(:,:,j)),trot)
        Dea(:,:,j)=matmul(matmul(transpose(trot),Dea(:,:,j)),trot)
        Deb(:,:,j)=matmul(matmul(transpose(trot),Deb(:,:,j)),trot)
 enddo
 !Dh=Dha+Dhb
 !De=Dea+Deb
else
 do j=1,NFrag
        Dh(:,:,j)=matmul(matmul(transpose(trot),Dh(:,:,j)),trot)
        De(:,:,j)=matmul(matmul(transpose(trot),De(:,:,j)),trot)
 end do
end if
print*, " "
print*, "Diabatic energies"
!if (Spin == "U")then
! do i=1,kstates
!  print*,"States of subpace", i
!  do j=1,state_vec(i)
!    if (alpbetsp(state_mat(i,j))==1) then
!    print*, "Alpha State",state_mat(i,j),"Enenergy:",H(state_mat(i,j),state_mat(i,j))
!    else if (alpbetsp(state_mat(i,j))==2) then
!    print*, "Beta  State",state_mat(i,j),"Enenergy:",H(state_mat(i,j),state_mat(i,j))
!    else
!    print*, "State",state_mat(i,j),"Enenergy:",H(state_mat(i,j),state_mat(i,j))
!    end if
!  end do
! end do
!else
 do i=1,kstates
 	print*,"States of subpace", i
 	do j=1,state_vec(i)
 		print*, "State",state_mat(i,j),"Enenergy:",H(state_mat(i,j),state_mat(i,j))
 	end do
 end do
!end if
print*, ""

! call JoinJac(H,rot,JD,1.0e-8)
! H=matmul(transpose(rot),matmul(H,rot))

print*, ""
print*, "Diabatic Hamiltonian after subspace diagonalization"
do i=1,Nstates
	write(*,200) H(i,1:i)
end do

JV=transpose(matmul(JV,trot))
print*, ""
print*, "Eigenvectors"
do i=1,Nstates
	print*, "state",i
	write(*,200) JV(1:Nstates,i)
end do
print*, ""
!!!!
!do i=1,Nstates
! print*, sum(abs(JV(i,:))**2)
!end do
!print*, 
!!!!!
if (Spin == "U" .or. soc == "Y") then
 allocate(hha(NStates))
 allocate(eea(NStates))
 allocate(hhb(NStates))
 allocate(eeb(NStates))
 print*, ""
 do i=1,Nfrag
        do j=1,Nstates
                hha(j)=Dha(j,j,i)
                eea(j)=Dea(j,j,i)
                hhb(j)=Dhb(j,j,i)
                eeb(j)=Deb(j,j,i)
        end do
        print*, "Hole population of Fragment with Alpha Spin", i
        write(*,400) hha(1:Nstates)
        print*, ""
        print*, "Hole population of Fragment with Beta Spin", i
        write(*,400) hhb(1:Nstates)
        print*, ""
        print*, "Particle population of Fragment with Alpha Spin", i
        write(*,400) eea(1:NStates)
        print*, ""
        print*, "Particle population of Fragment with Beta Spin", i
        write(*,400) eeb(1:NStates)
        print*, ""
 end do
else
 allocate(hh(NStates))
 allocate(ee(NStates))
 print*, ""
 do i=1,Nfrag
 	do j=1,Nstates
 		hh(j)=Dh(j,j,i)
 		ee(j)=De(j,j,i)
 	end do
        print*, "Hole population of Fragment", i
 	write(*,400) hh(1:Nstates)
 	print*, ""
        print*, "Particle population of Fragment", i
 	write(*,400) ee(1:NStates)
 	print*, ""
 end do
end if
 print*, ""

if (dip == 'Y') then
! print*, "Adiabatic velTEDM"
! do  i=1, Nstates
!        write(*,'(I5,3X,3(F9.5,3X))') i, DVPE(i,1:3)
! end do
! print*,
!
! print*, "Adiabatic TMDM"
! do  i=1, Nstates
!        write(*,'(I5,3X,3(F9.5,3X))') i, DPM(i,1:3)
! end do
! print*,

 DPE=matmul(JV,DPE)
 DPM=matmul(JV,DPM)
 DVPE=matmul(JV,DVPE)
 if (soc == 'Y') then
  print*,"The Transition Dipole Moments (TDM) are from singlets excited states:"
  print*, ""
  print*, "Diabatic eTDM"
  print*, ""
  write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
  do  i=1, Nstates/2
         write(*,'(I5,3X,3(F9.5,3X))') i, DPE(i,1:3)
  end do
  print*, ""

  print*, "Diabatic eTDM (vel)"
  print*, ""
  write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
  do  i=1, Nstates/2
        write(*,'(I5,3X,3(F9.5,3X))') i, DVPE(i,1:3)
  end do
  print*, ""

  print*, "Diabatic mTDM"
  print*, ""
  write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
  do  i=1, Nstates/2
         write(*,'(I5,3X,3(F9.5,3X))') i, DPM(i,1:3)
  end do
  print*, ""
 else
  print*, "Diabatic eTDM"
  print*, ""
  write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
  do  i=1, Nstates
         write(*,'(I5,3X,3(F9.5,3X))') i, DPE(i,1:3)
  end do
  print*, ""
 
  print*, "Diabatic eTDM (vel)"
  print*, ""
  write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
  do  i=1, Nstates
        write(*,'(I5,3X,3(F9.5,3X))') i, DVPE(i,1:3)
  end do
  print*, ""
 
  print*, "Diabatic mTDM"
  print*, ""
  write(*,'(1X,A6,6X,A1,11X,A1,11X,A1)')'States','X','Y','Z'
  do  i=1, Nstates
         write(*,'(I5,3X,3(F9.5,3X))') i, DPM(i,1:3)
  end do
  print*, ""
 end if
end if

if (trim(DenMat) == 'Y') then
 if (Spin == 'U' .or. soc == "Y") then
  allocate(TDMAa(Nbasis,Nbasis,Nstates),TDMAb(Nbasis,Nbasis,Nstates))
  allocate(TDMDa(Nbasis,Nbasis,Nstates),TDMDb(Nbasis,Nbasis,Nstates))
  allocate(auxe((Nbasis*(Nbasis-1))/2+Nbasis))
  allocate(auxh((Nbasis*(Nbasis-1))/2+Nbasis))
  TDMAa(:,:,:)=0.0
  TDMAb(:,:,:)=0.0
  TDMDa(:,:,:)=0.0
  TDMDb(:,:,:)=0.0
  !$OMP PARALLEL DEFAULT(SHARED) private(i)
  !$OMP DO
  do i=1,Nstates
        call Get_TDM(Nbasis,MOa,S,AM_Xa(:,:,i),AM_Ya(:,:,i),mo_oa,mo_va,TD,TDMAa(:,:,i))
        call Get_TDM(Nbasis,MOb,S,AM_Xb(:,:,i),AM_Yb(:,:,i),mo_ob,mo_vb,TD,TDMAb(:,:,i))
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(SHARED) private(i,j)
  !$OMP DO
  do i=1,Nstates
   do j=1,Nstates
    TDMDa(:,:,i)=TDMDa(:,:,i)+JV(i,j)*TDMAa(:,:,j)
    TDMDb(:,:,i)=TDMDb(:,:,i)+JV(i,j)*TDMAb(:,:,j)
   end do
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  allocate(Peta(Nbasis,Nbasis,Nstates),Petb(Nbasis,Nbasis,Nstates))
  allocate(Phta(Nbasis,Nbasis,Nstates),Phtb(Nbasis,Nbasis,Nstates))

  Peta(:,:,:)=0.0
  Phta(:,:,:)=0.0
  Petb(:,:,:)=0.0
  Phtb(:,:,:)=0.0

  do i=1, Nstates
        trha=0.0
        trea=0.0
        trhb=0.0
        treb=0.0
        call Get_PET_PHT(Nbasis,mo_oa,mo_va,MOa,TDMDa(:,:,i),Peta(:,:,i),Phta(:,:,i), trea, trha)
        !print*,"Trace state", i, "hole dens alpha:",trha,"particle dens alpha:",trea
        call Get_PET_PHT(Nbasis,mo_ob,mo_vb,MOb,TDMDb(:,:,i),Petb(:,:,i),Phtb(:,:,i), treb, trhb)
        !print*,"Trace state",i,"hole dens beta:",trhb,"particle dens beta:",treb
        !trh=trha+trhb
        !tre=trea+treb
        !print*,"Trace state",i,"hole dens:",trh,"particle dens:",tre
  end do

  open(13,file=trim(Inp)//"_hole.dens",form='unformatted')
  open(14,file=trim(Inp)//"_elec.dens",form='unformatted')
  
  Nelec=Neleca+Nelecb
  write(13) Nbasis, Nstates, Nelec
  write(14) Nbasis, Nstates, Nelec
  
  do a = 1, Nstates
        k=0
        auxe(:)=0.00
        auxh(:)=0.00
        !auxeb(:)=0.00
        !auxhb(:)=0.00
        do i =1, Nbasis
                do j=1,i
                        k=k+1
                        auxh(k)=Phta(i,j,a)+Phtb(i,j,a)
                        auxe(k)=Peta(i,j,a)+Petb(i,j,a)
                enddo
        enddo
        
  write(13) auxh
  write(14) auxe
  enddo
  close(13)
  close(14)
!Aqui acaba el Unrestricted
else

  allocate(TDMA(Nbasis,Nbasis,Nstates))
  allocate(TDMD(Nbasis,Nbasis,Nstates))
  allocate(auxe((Nbasis*(Nbasis-1))/2+Nbasis))
  allocate(auxh((Nbasis*(Nbasis-1))/2+Nbasis))

  TDMA(:,:,:)=0.0
  TDMD(:,:,:)=0.0
  !$OMP PARALLEL DEFAULT(SHARED) private(i)
  !$OMP DO
  do i=1,Nstates
        call Get_TDM(Nbasis,MO,S,AM_X(:,:,i),AM_Y(:,:,i),mo_o,mo_v,TD,TDMA(:,:,i))
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL
  !TDMA=TDMA/sqrt(2.0)

  !$OMP PARALLEL DEFAULT(SHARED) private(i,j)
  !$OMP DO
  do i=1,Nstates
   do j=1,Nstates
    TDMD(:,:,i)=TDMD(:,:,i)+JV(i,j)*TDMA(:,:,j)
   end do
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  allocate(Pet(Nbasis,Nbasis,Nstates))
  allocate(Pht(Nbasis,Nbasis,Nstates))

  Pet(:,:,:)=0.0
  Pht(:,:,:)=0.0

!!$omp parallel default(shared) private(i,j)
!!$omp do
  do i=1, Nstates
        trh=0.0
        tre=0.0
        call Get_PET_PHT(Nbasis,mo_o,mo_v,MO,TDMD(:,:,i),Pet(:,:,i),Pht(:,:,i), tre, trh)
        !print*,"Trace state",i,"hole dens:",trh,"particle dens:",tre
  enddo
!!$omp end do
!!$omp barrier
!!$omp end parallel
!do a=1,Nstates
!trh=0.0
!tre=0.0
!        do i= 1, mo_o
!                trh=trh+Pht(i,i,a)
!        enddo
!        do i=mo_o+1,mo_o+mo_v
!                tre=tre+Pet(i,i,a)
!        enddo
!        print*,"Trace state",a,"hole",trh,"particle",tre
!enddo
  open(13,file=trim(Inp)//"_hole.dens",form='unformatted')
  open(14,file=trim(Inp)//"_elec.dens",form='unformatted')

  write(13) Nbasis, Nstates, Nelec
  write(14) Nbasis, Nstates, Nelec

  do a = 1, Nstates
        k=0
        auxe(:)=0.00
        auxh(:)=0.00
        do i =1, Nbasis
                do j=1,i
                        k=k+1
                        auxh(k)=Pht(i,j,a)
                        auxe(k)=Pet(i,j,a)
                enddo
        enddo
  write(13) auxh
  write(14) auxe
  enddo
  close(13)
  close(14)
 end if
end if

200 format(5(E13.5,2X))
300 format(10(E13.5,2X))
400 format(5(F13.5,2X))

end program main
