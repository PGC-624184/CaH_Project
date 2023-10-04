set encoding iso_8859_1
set logscale y 10
set xlabel "Radial Distance ({\305})"
set ylabel "Wavelength Broadening ({\305})"
set title "Interatomic Calcium Broadening Interactions in the Photosphere of M Dwarf Stars"
plot "../Ca_H2_output.csv" using 1:($2*10) with points title "Ca H_2 System", "../Ca_H_output.csv" using 1:($2*10) with points title "Ca H System"

