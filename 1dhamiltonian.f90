!Auther Mukesh prakash khanore (email: mukeshkhanore@gmail.com)
! program to evaluate 1d Hamiltonian for Hubbard model (electron) with periodic boundary conditions and nearest neighbour interaction
!ns=nos of sites,ne=nos of elctron,	ts=total states, t= hoping amplitude, u=interaction term, sup= total upstates, sdn= total downstates
!it needs lapack library 'apt-get install libblas-dev liblapack-dev'
!	
program one_d
	integer::i,j,ne,ns,ts,temp,c1_up,a2_up,a1_up,k,l,m,c1_dn,a2_dn,a1_dn,inf,l1,u1
	real::u,t
	real*8,allocatable::hmatrix(:,:),eig(:),work(:),uterm(:),temp_h(:,:)
	integer,allocatable::states(:),sup(:),sdn(:)
	open(unit=150,file="input.txt",status="old",action="read")
	write(*,*) "give nos of sites,total nos of electron, total nos of states"
	read(150,*)ns,ne,ts
	write(*,*) "given nos of sites=",ns,"total nos of electron= ",ne, "total nos of states=",ts
	u=1.d0
	t=1.d0
	!read(*,*)U,t
	write(*,*)"U=",u,"t=",t

	t=-t
	
	allocate(states(ts),sup(ts),sdn(ts),hmatrix(ts,ts),eig(ts),work(ts*(3+ts/2)),uterm(ts),temp_h(ts,ts))
	open(unit=100,file='states_dni.txt',status='old',action='read')
	open(unit=110,file='states_upi.txt',status='old',action='read')
	!	open(unit=120,file='1dh_matrix.txt')
	do i=1,ts
		read(100,*)sdn(i)
		read(110,*)sup(i)	
		!states(i)=(2**ns)*sdn(i)+sup(i)
	enddo
	!write(*,*)"sup(ts)=",sup(ts)
	write(*,*)"---------------procecing----------------------------------------"
	hmatrix=0.d0
	uterm=0.d0
	!----------------------------------------------------------
	! calculating u term
	do i=1,ts
		temp=0
		do j=0,ns-1
			if(btest(sup(i),j).eqv..True.)then	
				if(btest(sup(i),j).eqv.btest(sdn(i),j))then
				temp=temp+1
				endif
			endif
		enddo
		hmatrix(i,i)=temp+hmatrix(i,i)
		uterm(i)=temp+uterm(i)
		!if(i==18)then
		!	write(*,*)uterm(i),sup(i),sdn(i)
		!endif
	enddo
	!----------------------------------------------------------------------   
	!calculating t term
	do k=1,ts
		do i=0,ns-1
	!---------calculating for up states---------------------------------------------------------------------      
			if(btest(sup(k),i).eqv..false.)then	!c_dagar_up "i" c_up "j"
				c1_up=ibset(sup(k),i)
				if(i+1>=ns)then
					a1_up=ibclr(c1_up,0)
				else
					a1_up=ibclr(c1_up,i+1)
				endif

				if(i-1<0)then
					a2_up=ibclr(c1_up,ns-1)
				else
					a2_up=ibclr(c1_up,i-1)
				endif

				do l=1,ts
					if(l.ne.k)then
						if(a1_up==sup(l).and.sdn(l)==sdn(k))then
							hmatrix(k,l)=t!+hmatrix(k,l)
						endif
						if(a2_up==sup(l).and.sdn(l)==sdn(k))then
							hmatrix(k,l)=t!+hmatrix(k,l)
						endif
					endif 
				enddo!l
			endif
	!---------calculating for down states-------------------------------------------------------------------
			if(btest(sdn(k),i).eqv..false.)then	!c_dagar_dn "i" c_dn "j"
				c1_dn=ibset(sdn(k),i)

				if(i+1>=ns)then
					a1_dn=ibclr(c1_dn,0)
				else
					a1_dn=ibclr(c1_dn,i+1)
				endif
				if(i-1<0)then
					a2_dn=ibclr(c1_dn,ns-1)
				else
					a2_dn=ibclr(c1_dn,i-1)
				endif

				do l=1,ts
					if(l.ne.k)then
						if(a1_dn==sdn(l).and.sup(l)==sup(k))then
							hmatrix(k,l)=t!+hmatrix(k,l)
						endif
						if(a2_dn==sdn(l).and.sup(l)==sup(k))then
							hmatrix(k,l)=t!+hmatrix(k,l)
						endif
					endif
				enddo
			endif
		enddo!i
	enddo!k

	!writing the hamiltonian matrix
	temp_h=hmatrix
	l1=ts*(3+ts/2)
	hmatrix=0.d0
	open(unit=1000,file="eig.txt")
	u1=u!50.d0	
	write(1000,*)u1
	write(1000,*)
	hmatrix=temp_h

	do i=1,ts
		hmatrix(i,i)=u1*uterm(i)
		!if(i==18)then
		!	write(*,*)u1*uterm(i),hmatrix(i,i),sup(i),sdn(i)
		!endif
	enddo

!	do i=1,ts
!		write(120,*)(int(hmatrix(i,j)),j=1,ts)
!	enddo

	write(*,*)"values of U=",u1

	call dsyev('V','U',ts,hmatrix,ts,eig,work,l1,inf)!lapack subroutine
	write(*,*)"INFO=", inf
	do j=1,ts	
		write(1000,*)eig(j)
	enddo
!enddo
!111 format(80I2,/)
end program one_d
