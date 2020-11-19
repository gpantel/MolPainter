import MDAnalysis as mda
import numpy as np
from scipy.spatial import distance
import toml 
import argparse

def main():
    parser = argparse.ArgumentParser(description='Read in the TOML-format input file')
    parser.add_argument("-i", "--input", type=str, help='TOML-format input file')
    parser.add_argument("-centerc", action='store_true', help='Center system at (x,y,z) = (0,0,0) after solvation')
    parser.add_argument("-zeroz", action='store_true', help='Raise the system such that min(z) = solvent spacing / 2 after solvation')
    args = parser.parse_args()
    
    inputs = toml.load(args.input)
    
    if 'solute' not in inputs: raise Exception('solute missing from input file')
    if 'output' not in inputs: raise Exception('output missing from input file')
    if 'shape' not in inputs['solute_lattice']: raise Exception('shape of solute_lattice missing from input file')
    if 'width' not in inputs['solute_lattice']: raise Exception('width of solute_lattice missing from input file')
    if 'height' not in inputs['solute_lattice']: raise Exception('height of solute_lattice missing from input file')
    if 'spacing' not in inputs['solute_lattice']: raise Exception('spacing of solute_lattice missing from input file')
    if 'spacing' not in inputs['solvent_lattice']: raise Exception('spacing of solvent_lattice missing from input file')
    if 'buffer' not in inputs['solvent_lattice']: raise Exception('buffer of solvent_lattice missing from input file')
    if 'lower_z_position' not in inputs['solvent_lattice']: raise Exception('lower_z_position of solvent_lattice missing from input file')
    if 'upper_z_position' not in inputs['solvent_lattice']: raise Exception('upper_z_position of solvent_lattice missing from input file')
    if 'paths' not in inputs['solvent_molecules']: raise Exception('paths list of solvent_molecules missing from input file')
    if 'numbers' not in inputs['solvent_molecules']: raise Exception('numbers list of solvent_molecules missing from input file')
    if 'max_iterations' not in inputs['solvent_molecules']: raise Exception('max_iterations missing from input file')
    if 'rotate' not in inputs['solvent_molecules']: raise Exception('rotate missing from input file')
    if len(inputs['solvent_molecules']['paths']) != len(inputs['solvent_molecules']['numbers']):
        raise Exception('There should be an equal number of entries for paths and numbers of solent_molecules')
    
    # 1) Load the solute to be solvated into a universe
    u = mda.Universe(inputs['solute'])
    layers_positions = u.atoms.positions
    
    # 2) Define lattice_width and lattice_height coming from the canvas
    lattice_width = inputs['solute_lattice']['width'] # 32
    lattice_height = inputs['solute_lattice']['height'] # 16
    lattice_spacing = inputs['solute_lattice']['spacing'] # 8 angstroms
    
    # 3) Define solvent lattice spacing
    solvent_lattice_spacing = inputs['solvent_lattice']['spacing']# 5 angstroms
    solvent_lattice_buffer  = inputs['solvent_lattice']['buffer'] # 4 angstroms
    
    solvent_lattice_width = int((lattice_width*lattice_spacing)//solvent_lattice_spacing)
    if inputs['solute_lattice']['shape'] == 'hex': solvent_lattice_height = int(((lattice_height*lattice_spacing)+(lattice_spacing/2.))//solvent_lattice_spacing)
    elif inputs['solute_lattice']['shape'] == 'square': solvent_lattice_height = int((lattice_height*lattice_spacing)//solvent_lattice_spacing)
    solvent_z1 = inputs['solvent_lattice']['lower_z_position'] # 17 angstroms
    solvent_z2 = inputs['solvent_lattice']['upper_z_position'] # 72 angstroms
    solvent_lattice_depth = int(np.abs(solvent_z2 - solvent_z1)//solvent_lattice_spacing)
    
    # 4) Make and array to contain solvent molecule IDs, where "0" is Empty, generate the lattice
    # 3-dimensional array specifying molecule IDs of each lattice point
    solvent_lattice = np.zeros((solvent_lattice_width,solvent_lattice_height,solvent_lattice_depth)).astype('int')
    # 4-dimensional array specifying positions of each lattice point
    solvent_lattice_xyz_coordinates = np.zeros((solvent_lattice_width,solvent_lattice_height,solvent_lattice_depth,3))
    for i in range(solvent_lattice_width):
        for j in range(solvent_lattice_height):
            for k in range(solvent_lattice_depth):
                x = (i*solvent_lattice_spacing)-(solvent_lattice_spacing/2)
                y = (j*solvent_lattice_spacing)-(solvent_lattice_spacing/2)
                z = (k*solvent_lattice_spacing)+solvent_z1
                solvent_lattice_xyz_coordinates[i][j][k] = np.array([x,y,z])
    
    # 5) Mark cells in the array as -1000 if molecules in the layers are within "lattice_spacing" distance of points
    # Check to see if any lattice site is within solvent_lattice_buffer distance to the solute -- if it is, mark the cell with -1000
    solvent_lattice_sites = np.where(solvent_lattice == 0)
    for i in range(len(solvent_lattice_sites[0])):
        xyz = solvent_lattice_xyz_coordinates[solvent_lattice_sites[0][i]][solvent_lattice_sites[1][i]][solvent_lattice_sites[2][i]]
        if np.amin(distance.cdist(layers_positions, [xyz])) < solvent_lattice_buffer:
            solvent_lattice[solvent_lattice_sites[0][i]][solvent_lattice_sites[1][i]][solvent_lattice_sites[2][i]] = -1000
    
    def rand_rotation_matrix():
        # This is a pythonic implementation of the random 3D rotation method described
        # in the Graphics Gems book series by Jim Arvo in 1991
        # http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
        randnums = np.random.uniform(size=(3,))
        theta, phi, z = randnums
        theta = theta*2.0*np.pi  # Rotation about the pole (Z).
        phi   = phi*2.0*np.pi  # For direction of pole deflection.
        z     = z*2.0  # For magnitude of pole deflection.
    
        r = np.sqrt(z)
        V = (np.sin(phi) * r,
             np.cos(phi) * r,
             np.sqrt(2.0 - z))
        st = np.sin(theta)
        ct = np.cos(theta)
        R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
        M = (np.outer(V, V) - np.eye(3)).dot(R)
        return M
    
    def create_solvent_universe(path, solvent_molecule_id, n_residues):
        solvu = mda.Universe(path)
        solvuxyz = np.copy(solvu.atoms.positions)
        solvent_positions = solvu.atoms.positions[np.where(solvu.atoms.tempfactors == 1)[0]]  
        solvent_cog = np.mean(solvent_positions, axis=0)
        solvuxyz = solvuxyz - solvent_cog
    
        mol_n_atoms = solvu.atoms.n_atoms
        n_atoms = mol_n_atoms*n_residues
        resindices = np.repeat(range(n_residues), mol_n_atoms)
        segindices = [0]*n_residues
    
        # initialize an empty universe and populate it
        solvent_universe = mda.Universe.empty(n_atoms, n_residues=n_residues, atom_resindex=resindices,
                                              residue_segindex=segindices, trajectory=True)
        solvent_universe.add_TopologyAttr('name', list(solvu.atoms.names)*n_residues)
        solvent_universe.add_TopologyAttr('type', list(solvu.atoms.types)*n_residues)
        solvent_universe.add_TopologyAttr('resnames', list(solvu.atoms.residues.resnames)*n_residues)
        solvent_universe.add_TopologyAttr('resid', list(range(1, (n_residues)+1)))
        solvent_universe.add_TopologyAttr('segid', [solvu.atoms.segids[0]])
        solvent_coordinates = np.concatenate(np.array([solvu.atoms.positions]*n_residues))
        return solvent_universe
    
    def lookup_lattice_xyz(lattice_xyz, lattice_sites, site_index):
        x, y, z = lattice_xyz[lattice_sites[0][site_index]][lattice_sites[1][site_index]][lattice_sites[2][site_index]]
        return x, y, z
    
    def assign_lattice_positions(path, solvent_lattice, solvent_molecule_id, rotation_matrices):
        solvu = mda.Universe(path)
        solvxyz = solvu.atoms.positions
        solvcom = solvu.atoms.center_of_mass()
    
        solvent_lattice_sites = np.array(np.where(solvent_lattice == solvent_molecule_id))
        # assign coordinates
        solvent_coordinates = []
        for i in range(len(solvent_lattice_sites[0])):
            xyz = np.copy(solvxyz - solvcom) # centered coordinates
            if inputs['solvent_molecules']['rotate'] == True:
                xyz = np.dot(xyz, rotation_matrices[i])
    
            x, y, z = lookup_lattice_xyz(solvent_lattice_xyz_coordinates, solvent_lattice_sites, i)
            xyz[:,0] += x
            xyz[:,1] += y
            xyz[:,2] += z
        
            solvent_coordinates.append(xyz)
        if len(solvent_lattice_sites[0]) > 1: solvent_coordinates = np.concatenate(solvent_coordinates)
        return solvent_coordinates
    
    def check_nearby_lattice_sites(n_range_x, n_range_y, n_range_z, col, row, dep, solvent_lattice):
        for x in n_range_x:
            for y in n_range_y:
                for z in n_range_z:
                    if solvent_lattice[col+x][row+y][dep+z] != 0:
                        return False
        return True
    
    def assign_empty_lattice_site(empty_lattice_sites):
        randnum = np.random.randint(0, high=len(empty_lattice_sites[0]))
        randcol = empty_lattice_sites[0,randnum]
        randrow = empty_lattice_sites[1,randnum]
        randdep = empty_lattice_sites[2,randnum]
        return randcol, randrow, randdep
    
    def insert_to_random_empty_point(path, solvent_lattice, solvent_molecule_id):
        solvu = mda.Universe(path)
        solvxyz = solvu.atoms.positions
        solvcom = solvu.atoms.center_of_mass()
        xyz = np.copy(solvxyz - solvcom) # centered coordinates
        if inputs['solvent_molecules']['rotate'] == True:
            rotation_matrix = rand_rotation_matrix()
            xyz = np.dot(xyz, rand_rotation_matrix())
        else: rotation_matrix = 0
    
        n_pos_x =  int(((abs(max(xyz[:,0]))+solvent_lattice_buffer) // solvent_lattice_spacing))
        n_neg_x = int(-((abs(min(xyz[:,0]))+solvent_lattice_buffer) // solvent_lattice_spacing))
        n_pos_y =  int(((abs(max(xyz[:,1]))+solvent_lattice_buffer) // solvent_lattice_spacing))
        n_neg_y = int(-((abs(min(xyz[:,1]))+solvent_lattice_buffer) // solvent_lattice_spacing))
        n_pos_z =  int(((abs(max(xyz[:,2]))+solvent_lattice_buffer) // solvent_lattice_spacing))
        n_neg_z = int(-((abs(min(xyz[:,2]))+solvent_lattice_buffer) // solvent_lattice_spacing))
        n_range_x = list(range(n_neg_x, n_pos_x+1))
        n_range_y = list(range(n_neg_y, n_pos_y+1))
        n_range_z = list(range(n_neg_z, n_pos_z+1))
    
        empty_lattice_sites = np.array(np.where(solvent_lattice == 0))
        # prune empty lattice sites that can't fit the excess space
        invalid_entries = []
        invalid_entries.append(np.where(empty_lattice_sites[0] <= -n_neg_x)[0])
        invalid_entries.append(np.where(empty_lattice_sites[1] <= -n_neg_y)[0])
        invalid_entries.append(np.where(empty_lattice_sites[2] <= -n_neg_z)[0])
        invalid_entries.append(np.where(empty_lattice_sites[0] >= (solvent_lattice.shape[0] - n_pos_x))[0])
        invalid_entries.append(np.where(empty_lattice_sites[1] >= (solvent_lattice.shape[1] - n_pos_y))[0])
        invalid_entries.append(np.where(empty_lattice_sites[2] >= (solvent_lattice.shape[2] - n_pos_z))[0])
        invalid_entries = np.unique(np.concatenate(invalid_entries))
        empty_lattice_sites = np.delete(empty_lattice_sites, invalid_entries, axis=1)
    
        if len(empty_lattice_sites[0]) == 0:
            raise Exception('Lattice full!')
        else:
            # roll random number to determine cell to try to insert to
            randcol, randrow, randdep = assign_empty_lattice_site(empty_lattice_sites)
            # check surrounding cells to ensure they are also empty. If not, re-roll
            iteration = 0
            while iteration < inputs['solvent_molecules']['max_iterations']:
                iteration += 1
                choice_boolean = check_nearby_lattice_sites(n_range_x, n_range_y, n_range_z, randcol, randrow, randdep, solvent_lattice)
                if choice_boolean == True:
                    for x in n_range_x:
                        for y in n_range_y:
                            for z in n_range_z:
                                if (x == 0) and (y == 0) and (z == 0):
                                    solvent_lattice[randcol][randrow][randdep] = solvent_molecule_id
                                else:
                                    solvent_lattice[randcol+x][randrow+y][randdep+z] = -solvent_molecule_id
                    break
                elif choice_boolean == False:
                    randcol, randrow, randdep = assign_empty_lattice_site(empty_lattice_sites)
            if iteration == inputs['solvent_molecules']['max_iterations']:
                raise Exception('Failed to insert to lattice after max_interations!')
            return solvent_lattice, rotation_matrix
    
    solvent_groups = []
    for i in range(len(inputs['solvent_molecules']['paths'])):
        solvent_molecule_id = i+1
        path = inputs['solvent_molecules']['paths'][i]
        num_solvent = inputs['solvent_molecules']['numbers'][i]
        print('Inserting %i molecules of '%(num_solvent) + path)
        solvent_universe = create_solvent_universe(path, solvent_molecule_id, num_solvent)
        rotation_matrices = []
        for i in range(num_solvent):
            solvent_lattice, rotation_matrix = insert_to_random_empty_point(path, solvent_lattice, solvent_molecule_id)
            rotation_matrices.append(rand_rotation_matrix())
    
        solvent_universe.atoms.positions = assign_lattice_positions(path, solvent_lattice, solvent_molecule_id, rotation_matrices)
        solvent_groups.append(solvent_universe.atoms)
    
    all_solvent_universe = mda.Merge(*solvent_groups)
    full_universe = mda.Merge(u.atoms, all_solvent_universe.atoms)
    
    # translate output coordinates if flags are defined
    if args.centerc == True:
        xyz = np.copy(full_universe.atoms.positions)
        xyz[:,0] -= np.mean(xyz[:,0])
        xyz[:,1] -= np.mean(xyz[:,1])
        xyz[:,2] -= np.mean(xyz[:,2])
        full_universe.atoms.positions = xyz
    if args.zeroz == True:
        xyz = np.copy(full_universe.atoms.positions)
        xyz[:,2] = (xyz[:,2] - np.amin(xyz[:,2])) + (solvent_lattice_spacing/2.)
        full_universe.atoms.positions = xyz
    
    xedge = np.amax(full_universe.atoms.positions[:,0]) - np.amin(full_universe.atoms.positions[:,0]) + solvent_lattice_spacing/2.
    yedge = np.amax(full_universe.atoms.positions[:,1]) - np.amin(full_universe.atoms.positions[:,1]) + solvent_lattice_spacing/2.
    zedge = np.amax(full_universe.atoms.positions[:,2]) - np.amin(full_universe.atoms.positions[:,2]) + solvent_lattice_spacing/2.
    
    full_universe.dimensions = np.array([xedge, yedge, zedge, 90., 90., 90.])
    full_universe.atoms.write(inputs['output'], bonds=None)

if __name__ == "__main__":
    main()
