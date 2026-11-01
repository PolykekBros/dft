    set terminal pngcairo enhanced font 'Verdana,12' size 800,600
    set output 'NiO_fatbands.png'
    set title "Orbital-Projected Fat Bands of NiO (DFT+U+SOC)"
    set ylabel "Energy - E_F (eV)"
    set yrange [-8:8]
    set xzeroaxis ls -1
    set grid y
    
    # Standard Brillouin zone labels from your array
    set xtics ("L" 0.0000, "{/Symbol G}" 1.0606, "X" 1.5909, "W" 2.3408, "{/Symbol G}" 2.7738)
    set arrow from 1.0606, graph 0 to 1.0606, graph 1 nohead lt 0
    set arrow from 1.5909, graph 0 to 1.5909, graph 1 nohead lt 0
    set arrow from 2.3408, graph 0 to 2.3408, graph 1 nohead lt 0
    
    # Base bands in gray, orbital weights as dynamic circles
    plot 'NiO_bands.dat.gnu' using 1:($2 - 11.4934) with lines lw 1.5 lc rgb 'gray' notitle,          'fatbands.dat' using 1:($2 - 11.4934):($3 > 0.05 ? $3 * 2.5 : 1/0) with points pt 7 lc rgb '#AA0000' ps variable title 'Ni 3d',          'fatbands.dat' using 1:($2 - 11.4934):($4 > 0.05 ? $4 * 2.5 : 1/0) with points pt 7 lc rgb '#0000AA' ps variable title 'O 2p'
