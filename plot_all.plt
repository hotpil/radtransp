#plot input data:
#set terminal wxt 0
set terminal pdfcairo
set output 'plots.pdf'

set xlabel 'x, cm'
set ylabel 'dens, cm^{-3}'
plot 'coeffs.dat' u 1:4 w l title 'electrons', 'coeffs.dat' u 1:5 w l title 'atoms'

#set terminal wxt 1
set xlabel 'x, cm'
set ylabel 'temp, eV'
plot 'coeffs.dat' u 1:6 w l title 'electrons', 'coeffs.dat' u 1:7 w l title 'atoms'

#set terminal wxt 2
set xlabel 'x, cm'
set ylabel 'sources'
plot 'coeffs.dat' u 1:2 w l title 's10', 'coeffs.dat' u 1:3 w l title 's01'

#plot results:
#set terminal wxt 3
set xlabel 'eps'
set ylabel '\psi'
plot for [i=2:100:10] 'profile.dat' u 1:i w l title sprintf("%d", i) 

#set terminal wxt 4
set xlabel 'eps'
set ylabel 'intens, arb.units'
plot for [i=2:100:10] 'iplus_omega_nmu.dat' u 1:i w l title sprintf("%d", i)

#set terminal wxt 5
set xlabel 'x, cm'
set ylabel 'total flux, erg /(cm^{2} s)'
plot 'total_intens.dat' u 1:2 w l title 'plus', 'total_intens.dat' u 1:3 w l title 'minus'

set xlabel 'x, cm'
set ylabel 'excited atoms, arb. units'
plot 'total_intens.dat' u 1:4 w l

set xlabel 'x, cm'
set ylabel 'i tilde, arb. units'
plot 'total_intens.dat' u 1:5 w l

unset output
unset terminal
unset xlabel
unset ylabel










