module line
use constants
implicit none      

contains

real(8) function profile(xind, eps)
	use globals
	use constants
	implicit none
	real(8) eps
	integer xind	
	
	profile = 0.0d0
	select case (profile_type)
		case (1)
		!Lorentz:
			profile = 1.0/(PI)*hw_Doppler(xind)/(hw_Doppler(xind)**2 + (eps - 1.0)**2)
			return
		case (2)
		!Gauss:
			profile = Doppler(hw_Doppler(xind), eps-1.0d0)
			return			
		case (3)
		!Voigt
			profile = Voigt(hw_Doppler(xind), hw_Holtsmark_ion(xind), eps - 1.0d0) 
			return
		case (4)
			!profile = holtsmark_aline((eps-1.0)/hw_Holtsmark_ion(xind))/hw_Holtsmark_ion(xind)
			profile = holtsmark((eps-1.0)/hw_Holtsmark_ion(xind))/hw_Holtsmark_ion(xind)
		case (5)
			profile = convol_vh(hw_Doppler(xind), hw_Holtsmark_ion(xind), hw_Holtsmark_ion(xind), eps - 1.0d0)
		case(6)
			profile = convol_hd(hw_Doppler(xind), hw_Holtsmark_ion(xind), eps - 1.0d0)
		case (7)
			if (abs(eps-1.0d0) < delta0) then
				profile = 1.0d0/(2*delta0)
			else
				profile = 0.0d0
			endif
		case (8)
			profile = Voigt(hw_Doppler(hw_const), hw_Holtsmark_ion(hw_const), eps - 1.0d0)
		case (9)
			!profile = Voigt(hw_Doppler_W(xind), hw_Lor_W(xind), eps - 1.0d0)
			profile = Doppler(hw_Doppler_W(xind), eps-1.0d0)
	end select
end function

real(8) function hw_Holtsmark_elec(xind)
!Normalized to omega_0 hw for electron Stark
!Static approximation
	use globals
	implicit none	
	integer xind
	!See Kotov et al.
	hw_Holtsmark_elec = 1.6e-22*v_dens_e(xind)/sqrt(v_temp_e(xind))
end function

real(8) function hw_Holtsmark_ion(xind)
!Normalized to omega_0 hw for ion Stark
!Static approximation
	use globals
	implicit none	
	integer xind
	!See Kotov et al.
	hw_Holtsmark_ion = 1.8d-15*v_dens_e(xind)**(2.0/3)	
end function

real(8) function hw_Doppler(xind)
!Normalized to omega_0
	use globals
	use constants
	implicit none
	real(8) atom_mass
	integer xind
	atom_mass = 2*H_MASS	
	!Dopler broadening:
	hw_Doppler = sqrt(2*v_temp_atom(xind)*eV_erg/atom_mass)/LIGHT_SPEED
	return
end function

real(8) function hw_Doppler_W(xind)
	use globals
	use constants
	implicit none
	real(8) atom_mass
	integer xind
	atom_mass = W_MASS	
	!Dopler broadening:
	hw_Doppler_W= sqrt(2*v_temp_atom(xind)*eV_erg/atom_mass)/LIGHT_SPEED
	return
end function

real(8) function hw_Lor_W(xind)
	use globals
	
	implicit none
	integer xind
	!See Vainstein broadening for non-H lines. 
	!Density in cm-2, temperature in eV
	
	hw_Lor_W = 8.47d-8*v_dens_e(xind)*sqrt(v_temp_e(xind))
	return
end function

real(8) function hw_Zeeman(xind)
	use globals
	use constants
	implicit none
	integer xind
	real(8) B_field
	B_field = 6.0
	!No use w/o splitting ...
	!Dopler broadening:
	hw_Zeeman = 4.3e-6*B_field
	return
end function

real(8) function holtsmark(beta)
!Holtsmark function by pure integration. A bit slow.
	use constants
	implicit none
	real(8), value ::  beta
	real(8) term, x, dx, sum
	real(8), parameter:: eps = 1.0d-8
	
	if (beta .eq. 0.0d0) then
		holtsmark = 0.0d0
		return
	endif

	if (beta < 0.0d0) then
		beta = -beta
	endif
	
	term = 0.0
	dx = 0.1*1.0/beta
	sum = 0.0
	x = dx
	do while ((abs(term) > eps) .or. (x < 1.0/beta))
		term = x * sin(beta*x) * exp(-x**1.5)
		sum = sum + term
		x = x + dx
	enddo
	holtsmark = sum*dx*2*beta/pi
end function holtsmark

real(8) function holtsmark_aline(beta)
!Holtsmark function fit by Aline
	implicit none
	real(8) beta
	beta = abs(beta)
	if (beta <= 0.5) then
		holtsmark_aline = 0.42441*beta**2*(1-0.463*beta**2)
		return
	else if (beta >= 3.5) then
		holtsmark_aline = 1.496*beta**(-2.5)*(1+5.107*beta**(-1.5))
		return
	else 
		holtsmark_aline = - 0.00578 - 0.10313*beta + 0.89176*beta**2 - 0.68127*beta**3 + 0.18818*beta**4 - 0.01804*beta**5
		return
	endif
end function holtsmark_aline

real(8) function holtsmark_fit(beta)
!fir from [Deane M. Peterson, 2014, Rational Approximations for the Microfield Distributions ...]
!BS does not work
	implicit none
	real(8), value :: beta
	beta = beta + 1.0d-10
	holtsmark_fit = 1.0d0*beta**(-8) + 276.2085d0*beta**(-6) &
		+ 693.5800d0*beta**(-4) + 325.2479d0*beta**(-2) &
		+ 23.18481d0*beta**(-1) + 48.37532d0 + 28.58249d0*beta &
		+ 1.5520050d0*beta**2.5 + 2.380155d0*beta**4 - 0.1407808d0*beta**5.5 &
		+ 0.1708035d0*beta**7
end function

real(8) function Voigt(hw_Doppler, hw_Holtsmark, t)
	implicit none
	real(8) hw_Doppler, hw_Holtsmark, t
	real(8) U, V
	logical flag

	call WOFZ (t/hw_Doppler, hw_Holtsmark/hw_Doppler, U, V, FLAG)
	Voigt = U/(sqrt(PI)*hw_Doppler)
end function Voigt

real(8) function Doppler(hw, t)
	use constants
	implicit none
	real(8), intent(in) :: t, hw

	Doppler = exp(-(t/hw)**2)/(hw*sqrt(pi))
end function

real(8) function holtsmark_distr(t, hwh)
	use funcs
	implicit none
	real(8), intent(in) :: t, hwh

	holtsmark_distr = holtsmark_aline(t/hwh)
end function

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
	dt = h1/100
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

real(8) function convol_vh(hw_Doppler, hw_Holtsmark_1, hw_Holtsmark_2, t)
!convolution of Voigt with Holtsmark distribution
	implicit none
	real(8), intent(in) :: hw_Doppler, hw_Holtsmark_1, hw_Holtsmark_2, t
	real(8) t1, dt, term, acc, y_old, y_new
	integer n
	
	n = 0
	t1 = 0.0d0
	term = 0.0d0
	dt = min(hw_Doppler, hw_Holtsmark_1, hw_Holtsmark_2)/100
	acc = 10.0d0
	y_new = 0.0d0
	y_old = 0.0d0
	do while (acc > 1.0d-3)
		term = holtsmark_distr(t1, hw_Holtsmark_2)*Voigt(hw_Doppler, hw_Holtsmark_1, t-t1)		
		term = term + holtsmark_distr(-t1, hw_Holtsmark_2)*Voigt(hw_Doppler, hw_Holtsmark_1, t+t1)
		y_new = y_new + term
		if (term /= 0) then
				acc = abs(y_old - y_new)/y_new
		endif
		y_old = y_new
		t1 = t1 + dt
	enddo
	convol_vh = dt*y_new
end function

real(8) function convol_hd(hw_Doppler, hw_Holtsmark, t)
!convolution of Doppler and Holtsmark distributions
	implicit none
	real(8), intent(in) :: hw_Doppler, hw_Holtsmark, t
	real(8) t1, dt, term, acc, y_old, y_new
	integer n
	
	n = 0
	t1 = 0.0d0
	term = 0.0d0
	dt = min(hw_Doppler, hw_Holtsmark)/100
	acc = 10.0d0
	y_new = 0.0d0
	y_old = 0.0d0
	do while (acc > 1.0d-3)
		term = holtsmark_distr(t1, hw_Holtsmark)*Doppler(hw_Doppler, t-t1)		
		term = term + holtsmark_distr(-t1, hw_Holtsmark)*Doppler(hw_Doppler, t+t1)
		y_new = y_new + term
		if (term /= 0) then
				acc = abs(y_old - y_new)/y_new
		endif
		y_old = y_new
		t1 = t1 + dt
	enddo
	convol_hd = dt*y_new
end function

real(8) function hw_Voigt(hw_Gauss, hw_Lorentz)
	real(8), intent(in) :: hw_Gauss, hw_Lorentz

	hw_Voigt = 0.5346*hw_Voigt + sqrt(0.2166*hw_Lorentz**2 + hw_Gauss**2)
end function

end module line


