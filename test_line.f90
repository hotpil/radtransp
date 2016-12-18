program test_line
use line
use funcs
implicit none
integer i
real(8) beta, d_beta, t, dt, hwd, hwh, tmax
real(8) hw_Doppler, hw_Holtsmark, eps, deps, eps_max
real(8), allocatable :: array(:,:)
integer, parameter :: file_unit = 13
!test bed for different subroutines

beta = 0.0d0;
d_beta = 1.0

hw_Doppler = 2.0d0
hw_Holtsmark = 4.0d0
deps = 0.1d0
eps_max = 5*hw_Doppler 
eps = -eps_max

open(unit=file_unit,file='spectr.dat')
do while (eps < eps_max) 
	write(file_unit, '(4(es16.8e3 a) es16.8e3 )') eps, ' ', Doppler(eps, hw_Doppler), &
	' ', holtsmark_distr(eps, hw_Holtsmark), ' ', convol_dd(Doppler, hw_Doppler, holtsmark_distr, hw_Holtsmark, eps), &
	' ', Doppler(eps, sqrt(hw_Holtsmark**2+hw_Doppler**2))
	eps = eps + deps
enddo      	
close(file_unit)

contains

!function convol_dd(hw_Doppler, hw_Holtsmark, t) result(y)
function convol_dd(f1, h1, f2, h2, t) result(y)
	implicit none
	real(8) y
	real(8), intent(in) :: h1, h2, t
	real(8) t1, dt, term, acc, y_old, y_new
	integer n

	interface
		real(8) function f1(x, half)
			real(8), intent(in)::x, half
		end function
		
		real(8) function f2(x, half)
			real(8), intent(in)::x, half
		end function
	end interface
	
	n = 0
	t1 = 0.0d0
	term = 0.0d0
	dt = hw_Doppler/100
	acc = 10.0d0
	y_new = 0.0d0
	y_old = 0.0d0
	do while (acc > 1.0d-3)
		term = f1(t1, h1)*f2(t-t1, h2)
		term = term + f1(-t1, h1)*f2(t+t1, h2)
		y_new = y_new + term
		acc = abs(y_old - y_new)/y_new
		y_old = y_new
		t1 = t1 + dt
	enddo
	y = y_new*dt
end function

function Doppler(t, hw_Doppler) result(y)
	use constants
	implicit none
	real(8) y
	real(8), intent(in) :: t, hw_Doppler
	y = exp(-(t/hw_Doppler)**2)/(hw_Doppler*sqrt(pi))
end function

function holtsmark_distr(t, hwh) result(y)
	use funcs
	implicit none
	real(8) y
	real(8), intent(in) :: t, hwh
	y = holtsmark_aline(t/hwh)
end function


end program test_line
