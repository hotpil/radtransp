module globals

integer :: nx, nmu, open_ind, g0, g1, neps, niterations, eps_integ_width, profile_type, show_mu, no_absorption, hw_const
real(8) :: L, E, gaunt, B01, B10, &
	omega_0, osc_force, A_emission, delta0, grel, iter_accuracy, mu0
real(8) :: deps, dmu

!v_intens(x, eps, mu)
real(8), allocatable :: intens_plus(:,:,:), veps(:), intens_minus(:,:,:), total_intens_plus(:), total_intens_minus(:), &
	v_dens_e(:), v_dens_atom(:), v_temp_atom(:), v_temp_e(:), v_x(:), I_tilde(:), v_mu(:), v_temp(:, :), v_profile(:, :)
	
character(128) SOLPS_file_name
integer, parameter :: unit_log_file=12

end