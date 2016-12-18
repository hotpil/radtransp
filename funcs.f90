module funcs
!functions which could be used anywhere

contains

real(8) function my_exp_int(u)  
![D.A. Barry /Journal of hydrology 227 (2000)]
	implicit none
	real(8) u
	real(8) :: G, b, hinf, q, h
	G = 5.61459484d-001
	b = 1.49906926d+000
	hinf = 7.50855456d-001
	q = 20.0/47*u**sqrt(31.0/26.0)
	h = 1.0/(1 + u*sqrt(u)) + hinf*q/(1+q)
	my_exp_int = exp(-u)*log(1+G/u-(1-G)/(h+b*u)**2)/ &
			 (G + (1-G)*exp(-u/(1-G)))
	return		 
end function

real(8) function integ(dt, v)
	!integrate v(:) with step dt
	implicit none
	real(8), allocatable, intent(in) :: v(:)
	real(8), intent(in) :: dt

	integer n
	n = size(v)	
	integ = dt*(sum(v(2:n-1)) + 0.5*v(1) + 0.5*v(n))	
end function

real(8) function integ_2(dt, x1, v)
	!integrate v(x1, x2) over x2 with step dt with constant x1
	implicit none
	real(8), allocatable, intent(in) :: v(:,:)
	real(8), intent(in) :: dt
	integer, intent(in) :: x1

	integ_2 = dt*(sum(v(x1,:)) - 0.5*v(x1, 1) - 0.5*v(x1, size(v, 2)))
end function

real(8) function trapz_2d(dx, dy, v)
!integrate v(:,:) over [a,b] x [c,d]
	implicit none
	real(8), allocatable, intent(in) :: v(:,:)
	real(8), intent(in) :: dx, dy
	integer :: nx, ny, i, j

	nx = size(v,1)
	ny = size(v,2)
	trapz_2d = 0.0d0
	
	!boundary
	trapz_2d = v(1,1) + v(1,ny) + v(nx,1) + v(nx, ny) &
		+ 2*sum(v(1, 2:ny-1) + v(nx, 2:ny-1)) &
		+ 2*sum(v(2:nx-1, 1) + v(2:nx-1, ny))

	!interior
	do i = 2, nx-1
		do j = 2, ny-1
			trapz_2d = trapz_2d + 4*v(i, j)
		enddo
	enddo

	trapz_2d = trapz_2d*0.25*dx*dy	
end function trapz_2d

real(8) function trapz_1d(x, y)
!integrate y over x
	implicit none
	real(8), allocatable, intent(in) :: x(:), y(:)
	integer nx, i
	
	nx = size(x)	
	trapz_1d = y(1)*(x(2)-x(1))*0.5d0
	do i=2, nx-1
		trapz_1d = trapz_1d + y(i)*(x(i+1) - x(i-1))*0.5d0
	enddo
	trapz_1d = trapz_1d + y(nx)*(x(nx) - x(nx-1))*0.5d0

end function trapz_1d	

real(8) function find_hw(x, y)
!define half-width of y(x) function
	implicit none
	real(8) , allocatable, intent(in) :: x(:), y(:)
	real(8) max_val
	integer nx, i

	nx = size(x)	
	max_val = 0.0d0
	do i = 1, nx	
		if (max_val < y(i)) then
			max_val = y(i)
		endif
	enddo
	
	max_val = max_val/2.0d0
	i = 1
	do while (y(i) < max_val)
		i = i + 1
	enddo
	find_hw = (1.0d0 - x(i))*2
end

	
	
end module funcs