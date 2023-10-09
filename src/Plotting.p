set xlabel "Radial Distance (angstrom)"
factor=45.56335252907954
far=422.673
set ylabel "4P to 4S Orbital Energy Difference (nm)"
set title "Interatomic Calcium Broadening Interactions in the Photosphere of M Dwarf Stars"
plot [3:20] "../Ca_H2_test.csv" using 1:(factor/$2) with linespoints title "Ca H_2 {/Symbol S}^{+} ", \
    "../Ca_H2_test.csv" using 1:(factor/$3) with linespoints title "Ca H_2 {/Symbol P}", \
    "../Ca_H2_test.csv" using 1:(factor/$4) with linespoints title "Ca H_2 {/Symbol P}^{*}", \
    "../Ca_H_test.csv" using 1:(factor/$2) with linespoints title "Ca H {/Symbol S}^{+}", \
    "../Ca_H_test.csv" using 1:(factor/$3) with linespoints title "Ca H {/Symbol P}"
