    set terminal pngcairo enhanced font 'Verdana,10' size 800,600
    set output 'NiO_pdos_mz.png'
    set grid
    set title "Projected Local Magnetization (m_z) of Ni 3d in NiO"
    set xlabel "Energy (eV)"
    set ylabel "PDOS (States / eV)"
    set format y "%.2f"
    
    plot "NiO.pdos_atm#1(Ni1)_wfc#4(d_j2.5)" using 1:5 title 'Ni1 3d j=5/2 (m_z)' with line lw 2 lc rgb 'red', \
         "NiO.pdos_atm#1(Ni1)_wfc#5(d_j1.5)" using 1:5 title 'Ni1 3d j=3/2 (m_z)' with line lw 2 lc rgb 'dark-red', \
         "NiO.pdos_atm#2(Ni2)_wfc#4(d_j2.5)" using 1:5 title 'Ni2 3d j=5/2 (m_z)' with line lw 2 lc rgb 'blue', \
         "NiO.pdos_atm#2(Ni2)_wfc#5(d_j1.5)" using 1:5 title 'Ni2 3d j=3/2 (m_z)' with line lw 2 lc rgb 'dark-blue'
