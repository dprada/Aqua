import parmed as pmd
import MDAnalysis as mda
import MyLab as mylab
import aqua as aqua


def MDAnalysis2Aqua (universe):

    if mylab.typename(universe) != 'MDAnalysis.core.AtomGroup.Universe':
        raise Exception('Argument type needed: MDAnalysis.core.AtomGroup.Universe')

    temp=aqua.msystem()
    
    for atom in universe.atoms:

        temp_atom=aqua.cl_unit()

        temp_atom.name=atom.name
        temp_atom.index=atom.index
        temp_atom.pdb_index=atom.index

        temp_atom.resid.name=atom.resname
        temp_atom.resid.index=atom.resid
        temp_atom.resid.pdb_index=atom.resid

        temp.atom.append(temp_atom)

    temp.num_atoms=len(temp.atom)
        
    for residue in universe.residues:

        temp_residue=aqua.cl_residue()

        temp_residue.index = residue.id
        temp_residue.pdb_index = residue.id
        temp_residue.name = residue.name
        temp_residue.list_atoms = []
        for atom in residue.atoms:
            temp_residue.list_atoms.append(atom.index)

        temp_residue.type = None

        temp.resid.append(temp_residue)

    if len(universe.bonds)>0:

        for bond in universe.bonds:
            at1=bond.atoms[0].index
            at2=bond.atoms[1].index
            temp.atom[at1].covalent_bonds.append(at2)
            temp.atom[at2].covalent_bonds.append(at1)

        cov_chains = temp.selection_covalent_chains(chain='ALL',select='ALL')
        temp.rebuild_chains(cov_chains)

    return temp

