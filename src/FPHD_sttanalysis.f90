module stateanalysis
      implicit none
contains
subroutine analyze_statesun(Nfrag,Nstates,kstates,state_vec,state_mat,Dha,Dhb,Dea,Deb,Dthr,alpbetsp)
implicit none
integer :: i,j,k,t
integer, intent(in) :: Nfrag,Nstates,kstates
integer, dimension(:), intent(out) :: state_vec,alpbetsp
integer, dimension(:,:), intent(out) :: state_mat
real*8, dimension(:,:,:), intent(in) :: Dha,Dea,Dhb,Deb
real*8 :: hea,heb,thr,Dthr

integer, dimension(:), allocatable :: cou

allocate(cou(kstates))

alpbetsp(:)=0
cou(:)=0
thr=2*Dthr
state_vec(:)=0.0
state_mat(:,:)=0

print*, "KEY OF SUBSPACES"
do j=1,NFrag
        print*, "Subspace",j,"is Frenkel localized in frag",j
end do
do j=1,Nfrag
do k=1,Nfrag
if ( j /= k ) then
if ( k > j ) then
        t=NFrag+(j-1)*(NFrag-1)+k-1
else
       t=NFrag+(j-1)*(NFrag-1)+k
end if
print*, "Subspace",t,"is CT from frag",j,"to Frag",k
end if
end do
end do
print*, "Subspace",kstates,"is undefined subspace"
print*, " "

print*, "NATURE ANALISIS"
do i=1,Nstates
        hea=0.0
        heb=0.0
        do j=1,Nfrag
                hea=abs(Dea(i,i,j))+abs(Dha(i,i,j))
                heb=abs(Deb(i,i,j))+abs(Dhb(i,i,j))
                if ( hea > thr ) then
                        state_vec(j)=state_vec(j)+1
                        cou(j)=cou(j)+1
                        state_mat(j,cou(j))=i
                        alpbetsp(i)=1

                       print*, "State",i,"located of frag",j, "for Alpha Spin"

                        go to 150
                else if ( heb > thr ) then
                        state_vec(j)=state_vec(j)+1
                        cou(j)=cou(j)+1
                        state_mat(j,cou(j))=i
                        alpbetsp(i)=2

                       print*, "State",i,"located of frag",j, "for Beta Spin"

                        go to 150
                end if
        end do

        do j=1,Nfrag
        do k=1,Nfrag
        if ( j /= k ) then
                hea=abs(Dea(i,i,k))+abs(Dha(i,i,j))
                heb=abs(Deb(i,i,k))+abs(Dhb(i,i,j))
                if ( hea > thr ) then
                        if ( k > j ) then
                                t=NFrag+(j-1)*(NFrag-1)+k-1
                        else
                                t=NFrag+(j-1)*(NFrag-1)+k
                        end if
                        print*, "State",i,"CT from Frag",j,"to Frag",k, "for Alpha Spin"
                        state_vec(t)=state_vec(t)+1
                        cou(t)=cou(t)+1
                        state_mat(t,cou(t))=i
                        alpbetsp(i)=1
                        go to 150
                else if ( heb > thr ) then
                        if ( k > j ) then
                                t=NFrag+(j-1)*(NFrag-1)+k-1
                        else
                                t=NFrag+(j-1)*(NFrag-1)+k
                        end if
                        print*, "State",i,"CT from Frag",j,"to Frag",k, "for Beta Spin"
                        state_vec(t)=state_vec(t)+1
                        cou(t)=cou(t)+1
                        state_mat(t,cou(t))=i
                        alpbetsp(i)=2
                        go to 150
                end if
        end if
        end do
        end do

        print*, "State",i,"has undefined nature for both Alpha and Beta Spin"
        state_vec(kstates)=state_vec(kstates)+1
        cou(kstates)=cou(kstates)+1
        state_mat(kstates,cou(kstates))=i
        alpbetsp(i)=0

150 continue
end do
print*, " "

return
end subroutine analyze_statesun

subroutine analyze_states(Nfrag,Nstates,kstates,state_vec,state_mat,Dh,De,Dthr)
implicit none
integer :: i,j,k,t
integer, intent(in) :: Nfrag,Nstates,kstates
integer, dimension(:), intent(out) :: state_vec
integer, dimension(:,:), intent(out) :: state_mat
real*8, dimension(:,:,:), intent(in) :: Dh,De
real*8 :: he,thr,Dthr

integer, dimension(:), allocatable :: cou

allocate(cou(kstates))

cou(:)=0
thr=2*Dthr
state_vec(:)=0.0
state_mat(:,:)=0

print*, "KEY OF SUBSPACES"
do j=1,NFrag
        print*, "Subspace",j,"is Frenkel localized in frag",j
end do
do j=1,Nfrag
do k=1,Nfrag
if ( j /= k ) then
if ( k > j ) then
        t=NFrag+(j-1)*(NFrag-1)+k-1
else
        t=NFrag+(j-1)*(NFrag-1)+k
end if
print*, "Subspace",t,"is CT from frag",j,"to Frag",k
end if
end do
end do
print*, "Subspace",kstates,"is undefined subspace"
print*, " "


print*, "NATURE ANALISIS"
do i=1,Nstates
        he=0.0
        do j=1,Nfrag
                he=abs(De(i,i,j))+abs(Dh(i,i,j))
                if ( he > thr ) then
                        state_vec(j)=state_vec(j)+1
                        cou(j)=cou(j)+1
                        state_mat(j,cou(j))=i

                        print*, "state",i,"located of frag",j

                        go to 150
                end if
        end do

        do j=1,Nfrag
        do k=1,Nfrag
        if ( j /= k ) then
                he=abs(De(i,i,k))+abs(Dh(i,i,j))
                if ( he > thr ) then
                        if ( k > j ) then
                                t=NFrag+(j-1)*(NFrag-1)+k-1
                        else
                                t=NFrag+(j-1)*(NFrag-1)+k
                        end if
                        print*, "State",i,"CT from Frag",j,"to Frag",k
                        state_vec(t)=state_vec(t)+1
                        cou(t)=cou(t)+1
                        state_mat(t,cou(t))=i
                        go to 150
                end if
        end if
        end do
        end do

        print*, "State",i,"has undefined nature"
        state_vec(kstates)=state_vec(kstates)+1
        cou(kstates)=cou(kstates)+1
        state_mat(kstates,cou(kstates))=i

150 continue
end do
print*, " "

return
end subroutine analyze_states


end module stateanalysis
