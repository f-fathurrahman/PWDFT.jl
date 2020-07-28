import ase.build
import ase.data.g2

for molname in ase.data.g2.data:
    atoms = ase.build.molecule(molname)
    atoms.set_pbc([True,True,True])
    atoms.center(vacuum=5.0)
    #cell = atoms.get_cell()
    #A = cell[0,0]
    #B = cell[1,1]
    #C = cell[2,2]
    #atoms.center(about=(0.5*A,0.5*B,0.5*C))
    atoms.write(molname + ".xyz")
    #atoms.write(molname + ".xsf")

