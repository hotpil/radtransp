module inter

contains

subroutine read_input(name)
	use globals
	implicit none

	integer, parameter :: unit_ofile=13
	character(len=40) :: buf, number
	character(len=40) :: name
	
	open(unit=unit_ofile, file=name, STATUS='OLD', iostat=open_ind)

	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') E
	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') gaunt
	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') osc_force
	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') A_emission
	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') delta0
	read(unit_ofile, *) buf, number
	read(number, '(i15)') g0
	read(unit_ofile, *) buf, number
	read(number, '(i15)') g1
	read(unit_ofile, *) buf, number
	read(number, '(i15)') eps_integ_width
	read(unit_ofile, *) buf, number
	read(number, '(i15)') neps
	read(unit_ofile, *) buf, number
	read(number, '(i15)') nmu
	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') mu0
	read(unit_ofile, *) buf, number
	read(number, '(i15)') niterations
	read(unit_ofile, *) buf, number
	read(number, '(es15.8)') iter_accuracy
	read(unit_ofile, *) buf, number
	read(number, '(i15)') profile_type
	read(unit_ofile, *) buf, SOLPS_file_name
	read(unit_ofile, *) buf, number
	read(number, '(i15)') show_mu
	read(unit_ofile, *) buf, number
	read(number, '(i15)') hw_const
	read(unit_ofile, *) buf, number
	read(number, '(i15)') no_absorption
	close(unit_ofile)
	
end subroutine read_input

subroutine print_input()
	use globals
	use constants
	implicit none
	
	write (unit_log_file,*) 'Input parameters are:'

	write (unit_log_file,'(a, es15.8)') 'E = ', E
	write (unit_log_file,'(a, es15.8)') 'gaunt = ', gaunt
	write (unit_log_file,'(a, es15.8)') 'osc_force = ', osc_force
	write (unit_log_file,'(a, es15.8)') 'A = ', A_emission
	write (unit_log_file,'(a, es15.8)') 'delta0 = ', delta0
	write (unit_log_file,'(a, i5)') 'g0 = ', g0
	write (unit_log_file,'(a, i5)') 'g1 = ', g1
	write (unit_log_file,'(a, i5)') 'eps_integ_width = ', eps_integ_width
	write (unit_log_file,'(a, i5)') 'neps = ', neps
	write (unit_log_file,'(a, i5)') 'nmu = ', nmu
	write (unit_log_file,'(a, es15.8)') 'mu0 = ', mu0
	write (unit_log_file,'(a, i5)') 'niterations = ', niterations
	write (unit_log_file,'(a, es15.8)') 'iter_accuracy = ', iter_accuracy
	write (unit_log_file,'(a, i5)') 'profile_type = ', profile_type
	write (unit_log_file,'(a, a)') 'SOLPS_file_name = ', SOLPS_file_name
	write (unit_log_file,'(a, i5)') 'show_mu= ', show_mu
	if (profile_type == 8) then
		write (unit_log_file,'(a, i5)') 'Using constant half-width at = ', hw_const
	endif
	
	if (no_absorption == 1) then
		write (unit_log_file,'(a)') 'Absorption is off !'
		write (*,'(a)') 'Absorption is off !'
	endif	
		
end subroutine print_input

subroutine save_before()
	use globals
	implicit none

	integer, parameter :: file_descr = 13
	integer i, j

	open(unit=file_descr,file='profile.dat')
	write(file_descr, '(a)') '#1: eps       2: profile(i,eps)'
	do j=1,neps
		write(file_descr, '(es16.8e3 a)', advance='no') veps(j), ' '
		do i=1, nx-1
			write(file_descr, '(es16.8e3 a)', advance='no') v_profile(i, j), ' '
		enddo
		write(file_descr, '(es16.8e3)') v_profile(i, j)
	enddo
	close(file_descr)	
end subroutine save_before	

subroutine save_results_mu(k)
	use globals
	implicit none
	
	character(128) :: intens_fname
	character(10) mu_name
	integer i, j, k
	integer, parameter :: file_descr = 14

	if (k < nmu) then
		write(mu_name, '(i5)') k
	else
		mu_name = 'nmu'
	endif
	mu_name =  trim(adjustl(mu_name))	
	
	intens_fname = 'intens_'//trim(adjustl(mu_name))//adjustl('.dat')
	open(unit=file_descr,file=intens_fname)
	write(file_descr, '(a es16.8e3)') '# mu = ', v_mu(k)
	do i=1,nx-1
		write(file_descr, '(2(es16.8e3 a) es16.8e3)') v_x(i), ' ', intens_plus(i, neps/2+1, k), &
			' ', intens_minus(i, neps/2+1, k)			
	enddo
	close(file_descr)
       
    !The whole intens array            
    intens_fname = 'iplus_omega_'//trim(adjustl(mu_name))//adjustl('.dat')
    open(unit=file_descr,file=intens_fname)
    do j=1,neps
        write(file_descr, '(es16.8e3 a)', advance='no') veps(j), ' '
        do i=1, nx-1
            write(file_descr, '(es16.8e3 a)', advance='no') intens_plus(i, j, k), ' '
        enddo
        write(file_descr, '(es16.8e3)') intens_plus(nx, j, k)
    enddo
    close(file_descr)       
            
    intens_fname = 'iminus_omega_'//trim(adjustl(mu_name))//adjustl('.dat')
    open(unit=file_descr,file=intens_fname)
    write(file_descr, '(a)') '#1: eps       2: I(v_x(i),eps)'
    do j=1,neps
        write(file_descr, '(es16.8e3 a)', advance='no') veps(j), ' '
        do i=1, nx-1
            write(file_descr, '(es16.8e3 a)', advance='no') intens_minus(i, j, k)*v_mu(k), ' '
        enddo 
        write(file_descr, '(es16.8e3)') intens_minus(nx, j, k)
    enddo
    close(file_descr)  
end subroutine

end module