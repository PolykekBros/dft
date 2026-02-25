# Gnuplot script file for plotting data in file "si.dos.dat"
# This file is called plot_dos.gnu
set terminal pngcairo enhanced font 'Verdana,10'
set output 'NiO_dos_dft_U.png'
set grid
set title "Density of states (DOS) of NiO crystal"
set xlabel "Energy (eV)"
set ylabel "DOS"
set format y "%.4f"
plot    "nio.dos_dft_U.dat" using 1:2 title 'DOS of spin-up' with line,\
	"nio.dos_dft_U.dat" using 1:(-$3) title 'DOS of spin-down' with line
