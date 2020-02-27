from ase import Atoms
from gpaw import GPAW, PW, FermiDirac
from ase.units import Bohr, Hartree
atoms = Atoms('Al', positions=[
    [      0.0000000000,       0.0000000000,       0.0000000000],
    ],
    pbc = [True, True, True],
    cell=[
    [     -2.0247899727,       0.0000000000,       2.0247899727],
    [      0.0000000000,       2.0247899727,       2.0247899727],
    [     -2.0247899727,       2.0247899727,       0.0000000000],
    ]
)
atoms.write('ATOMS.xsf')
ecutwfc = 15.000000*Hartree
calc = GPAW( mode=PW(ecutwfc),
             setups='hgh',
             xc='LDA_X+LDA_C_VWN',
             spinpol=False,
             occupations=FermiDirac(0.010000*Hartree),
             kpts={'size': (8, 8, 8), 'gamma':True},
             txt='-')
atoms.set_calculator(calc)
e1 = atoms.get_potential_energy()
print('\n!!!!!! Total energy      : %18.10f eV = %18.10f Hartree\n' % (e1,e1/Hartree))
