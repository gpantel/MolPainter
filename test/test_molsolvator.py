import pytest
from MolPainter.MolSolvator import solvate_and_export
import filecmp

class TestMolSolvator:
    
    def test_solvent_construction(self):
        args_tuple = ('./test/PDBfiles/solute+uppersolvent.pdb',\
        ['./test/PDBfiles/NA+.pdb', './test/PDBfiles/CL-.pdb', './test/PDBfiles/WF.pdb', './test/PDBfiles/W.pdb'],\
        [54, 54, 501, 4511],\
        './test/PDBfiles/solute+solvent_test.pdb',\
        32,16,8,5.3,5,48,24,-55,-16,7,False,False,False,1000,0)
        solvate_and_export(*args_tuple)
        refpdb = "./test/PDBfiles/solute+solvent_ref.pdb"
        testpdb = "./test/PDBfiles/solute+solvent_test.pdb"
        result = filecmp.cmp(refpdb, testpdb, shallow=False)
        assert result is True

