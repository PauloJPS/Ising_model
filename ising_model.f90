program ising
		implicit none

		real(8) :: B, J, T
		integer :: N, r
		integer :: n_ther, n_int

		real(8) :: M, E, X, C, Bin, chi    !thermodynamics variables
		real(8) :: sigma_M, sigma_E, sigma_X, sigma_C !thermodynamics std

		integer, allocatable :: lattice(:,:)
		integer :: i, s, len !counters

		character(len=20) :: string_N

		
		read(*,*) string_N
		read(string_N,*) N

		string_N = trim('metropolis'//trim(string_N)//'.txt')
		open(1, file=string_N)

		B = 0.
		J = 1.
		n_ther = 10000
		n_int  = 500000
		E = 0d0
		M = 0d0
		T = 1
		r = 2

		allocate(lattice(N,N))
		call create_lattice(N, lattice)
		do i=1, 30
			call thermalization(N, B, J, T, lattice, n_ther)
			call metropolis(N, B, J, T, lattice, n_int, E, M, C, X, Bin, sigma_E, sigma_M, sigma_C, sigma_X)
			write(1,*) E, M, C, X, Bin, sigma_E, sigma_M, sigma_C, sigma_X, T
			write(*,*)  T
			T = T + 0.1d0
		enddo
	
	contains 

	subroutine metropolis(N, B, J, T, lattice, n_int, E, M, C, X, Bin, sigma_E, sigma_M, sigma_C, sigma_X)
		implicit none

		integer, intent(in) :: N, n_int
		real(8), intent(in) :: B, J, T
		integer, dimension(N,N), intent(in out) :: lattice

		real(8), intent(out) :: M, E, X, C, Bin    !thermodynamics variables
		real(8), intent(out):: sigma_M, sigma_E, sigma_X, sigma_C !thermodynamics std
		
		integer :: i, k
		integer :: left, right, up, down
		integer :: si, sj
		real(8) :: lattice_sum
		real(8) :: E2, M2, C2, X2

		real(8) :: delta_E, E0, N2
		real(8) :: rand

		E0 = 0d0
		E = 0d0
		M = 0d0
		M2 = 0d0
		X = 0d0
		X2 = 0d0
		C = 0d0
		C2 = 0d0
		Bin = 0d0
		
		sigma_E = 0d0
		sigma_M = 0d0
		sigma_C = 0d0
		sigma_X = 0d0

		N2 = N*N

		do i=1,n_int

			call random_number(rand)
			si =  1 + floor(N*rand)
			call random_number(rand)
			sj =  1 + floor(N*rand)
			if (si == 1) then
      			left = lattice(N,sj)
      			right = lattice(2,sj)
   			else if (si == N) then
      			left = lattice(N-1,sj)
      			right = lattice(1,sj)
   			else
      			left = lattice(si-1,sj)
      			right = lattice(si+1,sj)
   			end if
			if (sj == 1) then
				up = lattice(si,2)
				down = lattice(si,N)
			else if (sj == N) then
				up = lattice(si,1)
				down = lattice(si,N-1)
			else
				up = lattice(si,sj+1)
				down = lattice(si,sj-1)
			end if

			delta_E = 2d0*lattice(si, sj)*(up + down + left + right)
			call random_number(rand)	
			if( rand < exp(-delta_E/T) ) then
				lattice(si, sj) = -lattice(si,sj)
			endif
			call energy(N, B, J, lattice, E0)

			lattice_sum = abs(sum(lattice))/N2
			E0 = E0/N2

			E = E + E0
			E2 = E2 + E0**2d0
			M = M + lattice_sum
			M2 = M2 + lattice_sum**2d0
			C = C + (E2/(i) - (E/(i))**2d0)/T/T
			C2 = C2 + (( E2/(i) - (E/(i))**2d0)/T/T)**2d0
			X = X + ((M2/(i) - (M/(i))**2d0)/T)
			X2 = X2 + (((M2/(i) - (M/(i))**2d0)/T))**2d0
			Bin = Bin + lattice_sum**4d0

		enddo


		E = E/n_int
		E2 = E2/n_int
		M = M/n_int
		M2 = M2/n_int
		X = X/n_int
		X2 = X2/n_int
		C = C/n_int
		C2 = C2/n_int
		Bin = 1d0 - Bin/M2/M2/3d0

		sigma_E = sqrt((E2 - E*E)/N_int)
		sigma_M = sqrt((M2 - M*M)/N_int)
		sigma_C = sqrt((C2 - C*C)/N_int)
		sigma_X = sqrt((X2 - X*X)/N_int)

	end subroutine metropolis
	
	subroutine energy(N, B, J, lattice, E)
		implicit none

		integer, intent(in) :: N
		real(8), intent(in) :: B, J
		integer, dimension(N,N), intent(in out) :: lattice

		real(8), intent(out) :: E

		integer :: si, sj, up, down, left, right

		E = 0d0

		do si=1, N
			do sj=1, N
				if (si == 1) then
      				left = lattice(N,sj)
      				right = lattice(2,sj)
				else if (si == N) then
					left = lattice(N-1,sj)
					right = lattice(1,sj)
				else
					left = lattice(si-1,sj)
					right = lattice(si+1,sj)
				end if
				if (sj == 1) then
					up = lattice(si,2)
					down = lattice(si,N)
				else if (sj == N) then
					up = lattice(si,1)
					down = lattice(si,N-1)
				else
					up = lattice(si,sj+1)
					down = lattice(si,sj-1)
				end if
				E = E - lattice(si,sj)*(left + right + up + down)
			enddo
		enddo
		E = E/2

	end	subroutine energy

	subroutine thermalization(N, B, J, T, lattice, n_ther)
		implicit none

		integer, intent(in) :: N, n_ther
		real(8), intent(in) :: B, J, T
		integer, dimension(N,N), intent(in out) :: lattice

		integer :: i, k
		integer :: left, right, up, down
		integer :: si, sj

		real(8) :: delta_E
		real(8) :: rand

		do i=1,n_ther
			call random_number(rand)
			si =  1 + floor(N*rand)
			call random_number(rand)
			sj =  1 + floor(N*rand)
			if (si == 1) then
      			left = lattice(N,sj)
      			right = lattice(2,sj)
   			else if (si == N) then
      			left = lattice(N-1,sj)
      			right = lattice(1,sj)
   			else
      			left = lattice(si-1,sj)
      			right = lattice(si+1,sj)
   			end if
			if (sj == 1) then
				up = lattice(si,2)
				down = lattice(si,N)
			else if (sj == N) then
				up = lattice(si,1)
				down = lattice(si,N-1)
			else
				up = lattice(si,sj+1)
				down = lattice(si,sj-1)
			end if

			delta_E = 2d0*lattice(si, sj)*(up + down + left + right)
			call random_number(rand)	
			if( rand < exp(-delta_E/T) ) then
				lattice(si, sj) = -lattice(si,sj)
			endif
		enddo
		

	end subroutine thermalization

	subroutine create_lattice(N, lattice)
		implicit none

		integer, intent(in) :: N
		integer, dimension(N,N), intent(in out) :: lattice

		real(8) :: rand
		integer :: i, j

		do i=1,N
			do j=1, n
				call random_number(rand)
				if( rand < 0.5)then
					lattice(i,j) = 1
				else
					lattice(i,j) = -1
				endif
			enddo
		enddo

	end subroutine create_lattice

	subroutine correlation_lenght(N, T, B, J, lattice,  r, N_int,  chi)
		implicit none

		real(8), intent(in) :: B, J, T
		integer, intent(in) :: N, r
		integer, intent(in) :: N_int
		integer, dimension(N,N), intent(inout) :: lattice

		real(8), intent(out) :: chi

		real(8) :: delta_E, sis, sisj
		integer :: i, k, w ! counters 
		integer :: si, sj

		real :: rand ! random number
		integer :: up, down, left, right

		chi = 0.
		sis = 0.
		sisj = 0.

		do i=1,N_int
			call random_number(rand)
			si =  1 + floor(N*rand)
			call random_number(rand)
			sj =  1 + floor(N*rand)
			if (si == 1) then
      			left = lattice(N-r,sj)
      			right = lattice(1+r,sj)
   			else if (si == N) then
      			left = lattice(N-r,sj)
      			right = lattice(r,sj)
   			else
      			left = lattice(si-r,sj)
      			right = lattice(si+r,sj)
   			end if
			if (sj == 1) then
				up = lattice(si,1+r)
				down = lattice(si,N-r)
			else if (sj == N) then
				up = lattice(si,r)
				down = lattice(si,N-r)
			else
				up = lattice(si,sj+r)
				down = lattice(si,sj-r)
			end if

			delta_E = 2d0*lattice(si, sj)*(up + down + left + right)
			call random_number(rand)	
			if( rand < exp(-delta_E/T) ) then
				lattice(si, sj) = -lattice(si, sj)
			endif
			
			sisj = sisj + lattice(si,sj)*(up + down + left + right)/4d0
			sis = sis + lattice(si, sj)
			
		enddo
	
		sisj = sisj/N_int
		sis = sis/N_int
		chi = abs(sisj - sis**2.)

	end subroutine correlation_lenght



end program ising
