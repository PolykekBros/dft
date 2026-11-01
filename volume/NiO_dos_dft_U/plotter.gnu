    set terminal pngcairo enhanced font 'Verdana,10' size 800,600
    set output 'NiO_pdos_sz_accurate.png'
    set grid
    set title "Exact S_z Projected DOS of Ni 3d in Relativistic NiO"
    set xlabel "Energy (eV)"
    set ylabel "Net S_z (States / eV)"
    set format y "%.2f"
    
    # Plotting Column 4 (Net Sz) against Column 1 (Energy)
    plot "Ni1_Sz_projected.dat" using 1:4 title 'Ni1 3d (Net S_z)' with line lw 2 lc rgb 'red', \
         "Ni2_Sz_projected.dat" using 1:4 title 'Ni2 3d (Net S_z)' with line lw 2 lc rgb 'blue'
