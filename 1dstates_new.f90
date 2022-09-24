! program to genrate 1d states for fermionic Hubbard Model
! mukesh prakash khanore mukeshkhanore (email: mpk22privacy@duck.com)
! this code genrate all possible states for fermions in binary system  k= nos of sites, ne = nos of electrons
!e.g. k=6,ne=3,d=56(111000),sr is a state configration of single row both up and down

Module global
	integer*4::i,j,temp,d,k,ne,rs,count,idn,iup
	integer*4,allocatable::Sr(:),sr_up(:),sr_dn(:)
endmodule global
    
program D_to_B
	use global

	write(*,*) "give nos of sites in single row,nos of electron"
	read(*,*)k,ne

	open(unit=80,file="input.txt",action='write')
	open(unit=90,file='states_i.txt',action='write')
	open(unit=100,file='states_dni.txt',action='write')
	open(unit=110,file='states_upi.txt',action='write')

	d=0
	count=0

	do i=0,2*k-1
		d=ibset(d,i)
	enddo

	write(*,*)d

	do j=1,d
		rs=0
		do i=0,2*k-1
			if(btest(j,i).eqv..true.)then
				rs=rs+1
			endif
		enddo
		if(rs==ne)then
			write(90,*)j
			count=1+count
		endif
	enddo

	close(unit=90)

	allocate(sr(count),sr_up(count),sr_dn(count))

	open(unit=120,file='states_i.txt',status='old',action='read')

	do i=1,count
		read(120,*)sr(i)
	enddo

	write(*,*)count

	do i=1,count
		idn=int(sr(i)/2**k)
		iup=mod(sr(i),2**k)
		sr_up(i)=iup
		sr_dn(i)=idn
		write(100,*)sr_dn(i)
		write(110,*)sr_up(i)
	enddo

	close(unit=100)
	close(unit=110)
	write(80,*)k,ne,count

end program D_to_B
