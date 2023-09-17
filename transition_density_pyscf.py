from pyscf import gto, scf, dft, tddft
mol = gto.Mole()
mol.build(
    atom = 'Ca 0 0 0; H 2.14 0 0 0',  # in Angstrom
    basis = 'ccpvdz',
    spin = 3,
    verbose = 4,
)

mf = dft.RKS(mol)
mf.xc = 'b3lyp'
mf.frozen = 20
mf.kernel()

mytd = tddft.TDDFT(mf)
mytd.nstates = 15
mytd.kernel()
mytd.analyze()