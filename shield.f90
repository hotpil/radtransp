program shield
	use globals
	use constants
	use inter
	implicit none

	integer :: k, ind_mu
	real(8) :: a, b, old_diff, new_diff, rel_diff
	real :: start, finish
	character(len=40) :: name

	!  Initialize program:
	open(unit=unit_log_file, file='log.txt')
	call get_command_argument(1, name)
	if (len_trim(name) == 0) then
		write(*,*) 'Usage: transp inp_file'
		call exit(0)
	endif

	call read_input(name)
	call print_input()
	call init()
    new_diff = 1.0d0

	! Save some parameters:
	call save_before()

	!   Solve the equation:
	do k=1,niterations
		call cpu_time(start)
		do ind_mu = 1, nmu-1
			call update_intens(ind_mu)			
		enddo
		!intens_plus(:,:,1:ind_mu-1) = 2*intens_plus(:,:,1:ind_mu-1)
		!intens_minus(:,:,1:ind_mu-1) = 2*intens_minus(:,:,1:ind_mu-1)
		call update_intens(nmu)
		call get_Itilde() 
		call cpu_time(finish)
		
      old_diff = new_diff
      new_diff = conv_check()
		if (old_diff /= 0) then
			rel_diff = abs(old_diff - new_diff)/old_diff
		else
			rel_diff = 0.0
			write(*,*) 'Iterations converged'
            exit
		endif
		
		write(*,'(a, i5, a, es16.8e3)') 'After interation ', k, ' the rel diff is ', rel_diff
		print '("Iteration time = ",f6.3," seconds")', finish-start
        if (rel_diff < iter_accuracy) then
            write(*,*) 'Iterations converged'
            exit
        endif
	enddo    
	
	call update_total_intens()

	!  Write output:
	new_diff=conv_check_2() !DEBUG

	call save_results()	
	close(unit_log_file)	
	!	All done
contains

subroutine update_intens(ind_mu)
	use globals	
	use constants
	implicit none
	integer, intent(in) :: ind_mu
	integer i, j
	real(8) dx
	
	do i=1,nx-1
			do j=1,neps			
				dx = v_x(i+1) - v_x(i)
				call rhs(i+1, j, a, b, I_tilde(i+1), v_mu(ind_mu))
				intens_plus(i+1, j, ind_mu) = (intens_plus(i, j, ind_mu) + dx*b)/(1 - a*dx)
				dx = v_x(nx-i+1) - v_x(nx-i)
				call rhs(nx-i, j, a, b, I_tilde(nx-i), v_mu(ind_mu))
				intens_minus(nx-i, j, ind_mu) = (intens_minus(nx-i+1, j, ind_mu) + dx*b)/(1 - a*dx)		
			enddo	
		enddo
end subroutine

subroutine init()
	use globals
	use constants
	use line
	implicit none
	integer :: i, j
	integer, parameter :: file_descr = 13
	
	nx = 0
	!Read from input SOLPS file and count number of lines:	
	open(unit=file_descr, file=SOLPS_file_name)
	DO
    READ (file_descr,*, END=10)
		nx = nx + 1
	END DO
10	rewind file_descr
	!write(*,*) 'nx = ', nx
	
	allocate(v_dens_e(nx))
	allocate(v_dens_atom(nx))
	allocate(v_temp_atom(nx))
	allocate(v_temp_e(nx))
	allocate(v_x(nx))	
	
	do i=1, nx
	!x te_11 tatm_11 ne_11 natm_11
		read(file_descr, *) v_x(i), v_temp_e(i), v_temp_atom(i), v_dens_e(i), v_dens_atom(i)		
	enddo	
	close(unit=file_descr)	

	!write(*,'(a es16.8e3)') 'omega_0 = ', omega_0
	omega_0 = E/hpl_eV_s
	grel = real(g0)/g1 + 1.0d0
	B10 = A_emission*4*PI**3*LIGHT_SPEED**2/(hpl*omega_0**3)
	B01 = real(g1)/g0*B10

	allocate(intens_plus(nx, neps, nmu))
	intens_plus(1:nx, 1:neps, 1:nmu) = 0.0d0
	allocate(intens_minus(nx, neps, nmu))
	intens_minus(1:nx, 1:neps, 1:nmu) = 0.0d0
	allocate(total_intens_plus(nx))
	total_intens_plus(1:nx) = 0.0d0
	allocate(total_intens_minus(nx))
	total_intens_minus(1:nx) = 0.0d0	

	allocate(veps(neps))
	veps(1) = 1.0d0 - delta0*eps_integ_width
	deps = 2*delta0*eps_integ_width/(neps - 1)
	do i=1, neps - 1
		veps(i+1) = veps(i) + deps	
	enddo	
	
	dmu = (1.0d0 - mu0)/(nmu-1)
	allocate(v_mu(nmu))
	v_mu(1) = mu0
	do i=1, nmu-1
		v_mu(i+1) = v_mu(i) + dmu
	enddo		
	
	!A temporary vector used for I_tilde calculation:
	allocate(v_temp(neps, nmu)) 
	
	allocate(I_tilde(nx))
	I_tilde(1:nx) = 0.0d0
	
	!generate profile matrix
	write(*,'(a)', advance='no') "Generating profile matrix ..."
	allocate(v_profile(1:nx, 1:neps))
	do i=1, nx
		do j=1, neps
			v_profile(i, j) = profile(i, veps(j))
		enddo
	enddo
	write(*,*) "Done."
end subroutine init

subroutine update_total_intens()
	use globals
	use funcs
	use constants
	implicit none	
	integer i, j

	do i=1,nx
		do j=1, nmu
			v_temp(1:neps, j) = intens_plus(i,1:neps,j)*v_mu(j)
		enddo
		total_intens_plus(i) = trapz_2d(deps, dmu, v_temp)		
		do j=1, nmu
			v_temp(1:neps, j) = intens_minus(i,1:neps,j)*v_mu(j)
		enddo
		total_intens_minus(i) = trapz_2d(deps, dmu, v_temp)
	enddo	

	total_intens_plus = total_intens_plus*2*PI*omega_0
	total_intens_minus = total_intens_minus*2*PI*omega_0
end subroutine update_total_intens

subroutine get_Itilde()
	use globals
	use funcs
	implicit none	
	integer i, j

	!left boundary
	do j=1, neps
		v_temp(j, 1:nmu) = intens_minus(1, j, 1:nmu)*v_profile(1, j)	
	enddo
	I_tilde(1) = trapz_2d(deps, dmu, v_temp)*0.5d0
	
	!bulk
	do i=2, nx-1
		do j=1, neps
			v_temp(j, 1:nmu) = (intens_plus(i+1, j, 1:nmu) + intens_minus(i-1, j, 1:nmu))*v_profile(i, j)
		enddo
		I_tilde(i) = trapz_2d(deps, dmu, v_temp)*0.5d0
	enddo
	
	!rigth boundary
	do j=1, neps
		v_temp(nx, 1:nmu) = intens_plus(nx, j, 1:nmu)*v_profile(nx, j)	
	enddo
	I_tilde(nx) = trapz_2d(deps, dmu, v_temp)*0.5d0		
end subroutine


real(8) function V01(xind)
	use globals
	use funcs
	implicit none
	integer xind	

	V01 = 1.58114e-5*osc_force/(E*sqrt(v_temp_e(xind)))*my_exp_int(E/v_temp_e(xind))*gaunt		
	return
end function

real(8) function V10(xind)
	use globals
	implicit none	
	integer xind		
	
	V10 = real(g0)/g1*V01(xind)*exp(E/v_temp_e(xind))
	return
end function

real(8) function kappa(xind, ind_eps)
	use globals
	use constants
	implicit none
	integer :: xind, ind_eps
	
	kappa = real(g1)/g0*(PI*LIGHT_SPEED)**2/omega_0**3*A_emission &
				*v_dens_atom(xind)*v_profile(xind, ind_eps)
	return
end function

real(8) function get_n1(xind, iwave)
	use globals
	implicit none
	real(8), value :: iwave
	integer xind
	
	if (no_absorption == 1) then
		iwave = 0.0d0 
	endif
	
	get_n1 = (B01*iwave + v_dens_e(xind)*V01(xind))/ &
		(A_emission + iwave*(B01 + B10) + v_dens_e(xind)*(V10(xind) + V01(xind)))
	get_n1 = get_n1*v_dens_atom(xind)
end function

real(8) function get_n0(xind, iwave)
	use globals
	implicit none
	real(8), value :: iwave
	integer xind
		
	get_n0 = v_dens_atom(xind) - get_n1(xind, iwave)
end function

subroutine rhs(xind, ind_eps, a, b, iwave, mu)
	use globals
	use constants
	implicit none
	real(8), intent(in) :: iwave, mu 
	integer, intent(in):: xind, ind_eps	
	real(8), intent(out) :: a, b	
	real(8) tmp
	
	!equation is:
	! dy/dx = a*y + b
	
	tmp = hpl/(4*PI)*v_profile(xind, ind_eps)/mu
	a = -tmp*(get_n0(xind, iwave)*B01 - get_n1(xind, iwave)*B10)
	b = get_n1(xind, iwave)*A_emission*tmp	
	
	if (no_absorption == 1) then
		a = 0.0d0
	endif
end subroutine

function conv_check() result(diff_max)
	use globals
	implicit none
	integer j_middle, i
	real(8) diff_max
	real(8) lhs, a, b, diff_plus, diff_minus, dx
	integer, parameter :: file_descr = 13
	
	j_middle = neps/2+1
	diff_max = 0.0d0
	open(unit = file_descr, file = 'converg.dat')
	write(file_descr, *) "#==========================="

	do i=2, nx-1		
		dx = v_x(i+1) - v_x(i)
		call rhs(i, j_middle, a, b, I_tilde(i), v_mu(nmu))
		lhs = (intens_plus(i, j_middle, nmu) - intens_plus(i-1, j_middle, nmu))/dx
		diff_plus = lhs - a*intens_plus(i, j_middle, nmu) - b
		diff_plus = abs(diff_plus/intens_plus(i, j_middle, nmu))
		
		call rhs(nx-i, j_middle, a, b, I_tilde(nx-i), v_mu(nmu))
		lhs = -(intens_minus(nx-i+1, j_middle, nmu) - intens_minus(nx-i, j_middle, nmu))/dx
		diff_minus = lhs - a*intens_minus(nx-i, j_middle, nmu) - b	
		diff_minus = abs(diff_minus/intens_minus(nx-i, j_middle, nmu))
		
		if (diff_max < diff_plus) then
			diff_max = diff_plus
		endif
		if (diff_max < diff_minus) then
			diff_max = diff_minus
		endif
		
		write(file_descr, '(2(es16.8e3 a) es16.8e3)') v_x(i), ' ', diff_plus, ' ', diff_minus		
	enddo
	
	close(unit=file_descr)	
end function conv_check

function conv_check_2() result(diff_max)
	use globals
	implicit none
	real(8) diff_max
	real(8) :: diff, t1, t2, t3, t4
	integer :: i, k
	integer, parameter :: file_descr = 13
	
	open(unit = file_descr, file = 'converg.dat')
	write(file_descr, *) "#==========================="
	diff_max = 0.0d0

	do i=2, nx-1		
		t1 = get_n1(i, I_tilde(i))*A_emission
		t2 =  I_tilde(i)*(get_n0(i, I_tilde(i))*B01 - get_n1(i, I_tilde(i))*B10) 
		t3 = V01(i)*get_n0(i, I_tilde(i))*v_dens_e(i) 
		t4 = V10(i)*get_n1(i, I_tilde(i))*v_dens_e(i)			
		diff = t2 + t3 - t1 - t4		
			
		if (diff_max < diff) then
			diff_max = diff
			k = i
		endif			
	enddo
	
	write(6, '(a, es16.8e3, a, i15)') 'Maximum difference ', diff_max, ' at ', k
	write(6, '(4es16.8e3)') t1, t2, t3, t4 !DEBUG
	close(unit=file_descr)	
end function conv_check_2


subroutine transp_info(descr)
	use globals
	integer, intent(in) :: descr
	real(8) :: Q_e, Q_deex, Q_spont, Q_ind, Q_abs, tot_out
	
	call check_sources(Q_e, Q_deex, Q_spont, Q_ind, Q_abs)	
	tot_out = total_intens_plus(nx) + total_intens_minus(1)
	write(descr,'(a)') '==== Transport info ==='
	write(descr,'(a, es15.8)') 'Electron impact source ', Q_e
	write(descr,'(a, es15.8)') 'Electron deexcitation source ', Q_deex	
	write(descr,'(a, es15.8)') 'Spontanious emission ', Q_spont
	write(descr,'(a, es15.8)') 'Induced emission ', Q_ind	
	write(descr,'(a, es15.8)') 'ABsorbed radiation', Q_abs	
	write(descr,'(a, es15.8)') 'Check (zero) ', Q_e - Q_deex - Q_spont - Q_ind + Q_abs	
	write(descr,'(a, es15.8)') 'Intensity plus ', total_intens_plus(nx)
	write(descr,'(a, es15.8)') 'Intensity minus', total_intens_minus(1)
	write(descr,'(a, es15.8)') 'Rel Balance (zero) ', (Q_e - Q_deex - mu0*Q_spont - tot_out)/tot_out
	write(descr,'(a, es15.8)') 'Assymetry factor, plus/minus ', total_intens_plus(nx)/total_intens_minus(1)	
end subroutine transp_info

subroutine check_sources(Q_e, Q_deex, Q_spont, Q_ind, Q_abs)
	use globals
	use constants
	use funcs
	integer i
	real(8), intent(out) :: Q_e, Q_deex, Q_spont, Q_ind, Q_abs
	real(8), allocatable :: tmp(:)
	
	allocate(tmp(nx))
	do i=1, nx
		tmp(i) = V01(i)*get_n0(i, I_tilde(i))*v_dens_e(i)
	enddo
	Q_e = trapz_1d(v_x, tmp)*hpl*omega_0
		
	do i=1, nx
		tmp(i) = V10(i)*get_n1(i, I_tilde(i))*v_dens_e(i)
	enddo	
	Q_deex = trapz_1d(v_x, tmp)*hpl*omega_0
	
	do i=1, nx
		tmp(i) = get_n1(i, I_tilde(i))
	enddo	
	Q_spont = trapz_1d(v_x, tmp)*A_emission*hpl*omega_0	
	
	do i=1, nx
		tmp(i) = I_tilde(i)*get_n1(i, I_tilde(i))*B10
	enddo	
	Q_ind = trapz_1d(v_x, tmp)*hpl*omega_0
	
	do i=1, nx
		tmp(i) = I_tilde(i)*get_n0(i, I_tilde(i))*B01
	enddo	
	Q_abs = trapz_1d(v_x, tmp)*hpl*omega_0
end subroutine

real(8) function cond_elec(i)
	use constants
	use globals
	real(8) log_Coulomb
	integer, intent(in) :: i
	
	!See NRL plasma formulary
	
	!Coulomb logarihtm 
	!if (v_temp_e(i) < 50.0d0) then
	!Since we are always here
		log_Coulomb = 23.4d0 - 1.15d0*log10(v_dens_e(i)) + 3.45d0*log10(v_temp_e(i))
	!endif
	
	cond_elec =1.106d6*v_temp_e(i)**2.5*eV_erg/(log_Coulomb*E_MASS)
end function

real(8) function intens_cond(i)
	use constants
	use globals
	integer, intent(in) :: i
	
	intens_cond = - cond_elec(i)*eV_erg*(v_temp_e(i) - v_temp_e(i-1))*(v_x(i) - v_x(i-1))
end function 

subroutine save_results()
	use globals
	use constants
	use inter
	use line
	implicit none

	integer, parameter :: file_descr = 13
	integer i, k	
	!real(8), allocatable :: tmp(:)
	
	!allocate(tmp(neps))	
	
	!Save intensities along diferent directions
	!Left and right beams (k=nmu) are considered separately since are not always accessed in the loop.
	do k = 1, nmu-1, show_mu    				 
		call save_results_mu(k)				
	enddo	
	call save_results_mu(nmu)	

	!Total energy flux and other parameters
	open(unit=file_descr, file='total_intens.dat')
	write(file_descr, '(a)') '#1: x       2: intens_plus  3: intens_minus   4: n1 5: I_tilde' !6: u1_Boltzman 7: intens_Planck'
	do i=1,nx	
		write(file_descr, '(4(es16.8e3 a) es16.8e3)') v_x(i), ' ',  total_intens_plus(i), ' ', total_intens_minus(i), &
				' ', get_n1(i, I_tilde(i)), ' ', I_tilde(i)*B01
				!1.0d0/(1.0d0 + real(g0)/g1*exp(E/v_temp_atom(i))), ' ',  hpl*omega_0**3/(4*PI**3*LIGHT_SPEED**2)/(exp(E/v_temp_atom(i) - 1.0d0))
	enddo
	close(file_descr)
	
	open(unit=file_descr,file='coeffs.dat')
	write(file_descr, '(a)') '#1: x            &
	2: S10         &
	3: S01          &
	4: dens_e       &
	5: dens_a       &
	6: temp_e       &
	7: temp_a       &
	8: hw_Holtsmark_elec &
	9: hw_Holtsmark_ion &
	10: hw_Doppler &
    11: kappa(at eps=1) &
	12: Voigt hw'
	do i=1,nx
		!tmp= v_profile(i, :)
		write (file_descr,'(12(es16.8e3 a) es16.8e3)') &
		v_x(i), ' ', V10(i)*v_dens_e(i), ' ', V01(i)*v_dens_e(i), ' ', &
		v_dens_e(i), ' ', v_dens_atom(i), ' ', v_temp_e(i), ' ', &
		v_temp_atom(i), ' ', hw_Holtsmark_elec(i), ' ', &
		hw_Holtsmark_ion(i), ' ', hw_Doppler(i), ' ', -1.0d0, ' ', hw_Voigt(hw_Doppler(i), hw_Holtsmark_ion(i))
	enddo
	close(file_descr)
	
	!electron conductivity flux
	open(unit=file_descr, file='electron_cond_flux.dat')
	write(file_descr, '(a)') '#1: x       2: flux, erg/cm2s'
	do i=2,nx
		write(file_descr, '(2(es16.8e3 a) es16.8e3)') v_x(i), ' ', intens_cond(i), ' ', cond_elec(i)				
	enddo
	close(file_descr)
	
	call transp_info(6) !print on screen
	call transp_info(unit_log_file) !save to log file	
end subroutine save_results

end program
