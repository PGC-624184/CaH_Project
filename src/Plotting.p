set xlabel "Radial Distance (angstrom)"
set ylabel "4P to 4S Orbital Energy Difference (nm)"
set title "Interatomic Calcium Broadening Interactions in the Photosphere of M Dwarf Stars"
plot "../Ca_H2_output.csv" using 1:2 with points title "Ca H_2 {/Symbol S}^{+}", \
    "../Ca_H2_output.csv" using 1:3 with points title "Ca H_2 {/Symbol P}^{*}", \
    "../Ca_H2_output.csv" using 1:4 with points title "Ca H_2 {/Symbol P}", \
    "../Ca_H_output.csv" using 1:2 with points title "Ca H {/Symbol S}^{+}", \
    "../Ca_H_output.csv" using 1:3 with points title "Ca H P_y", \
    "../Ca_H_output.csv" using 1:4 with points title "Ca H P_z"
