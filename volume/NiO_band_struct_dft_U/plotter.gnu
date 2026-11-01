    set terminal pngcairo enhanced font 'Verdana,12' size 800,600
    set output 'NiO_bands_relativistic.png'
    
    set title "Relativistic Band Structure of NiO (DFT+U+SOC)"
    set ylabel "Energy - E_F (eV)"
    set yrange [-8:8]  # Focus specifically on the Mott gap region
    set xzeroaxis ls -1
    set grid y
    
    set xtics ("L" 0.0000, "{/Symbol G}" 1.0606, "X" 1.5909, "W" 2.3408, "{/Symbol G}" 2.7738)
    
    # Add vertical dashed lines at the high-symmetry points for clarity
    set arrow from 1.0606, graph 0 to 1.0606, graph 1 nohead lt 0
    set arrow from 1.5909, graph 0 to 1.5909, graph 1 nohead lt 0
    set arrow from 2.3408, graph 0 to 2.3408, graph 1 nohead lt 0
    
    # Plotting the pre-formatted GNU output from bands.x
    plot 'NiO_bands.dat.gnu' using 1:($2 - 11.4934) with lines lw 1.5 lc rgb 'dark-blue' notitle
