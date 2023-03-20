import pytest
from MolPainter.MolSolvator import solvate_and_export
import filecmp

class TestMolSolvator:
    
    def test_solvent_construction(self):
        args_tuple = ('./test/PDBfiles/DPPC.pdb',\
        ['./test/PDBfiles/NA+.pdb', './test/PDBfiles/CL-.pdb', './test/PDBfiles/WF.pdb', './test/PDBfiles/W.pdb'],\
        [1, 1, 9, 90],\
        './test/PDBfiles/solvator_test.pdb',\
        4, 4, 8, 4.8, 4.5, 6, 6, 0, 35, 7, False, False, False, 1000, 0)
        solvate_and_export(*args_tuple)
        refpdb = "./test/PDBfiles/solvator_ref.pdb"
        testpdb = "./test/PDBfiles/solvator_test.pdb"
        result = filecmp.cmp(refpdb, testpdb, shallow=False)
        assert result is True

