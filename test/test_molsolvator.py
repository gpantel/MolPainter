import pytest
from MolPainter.MolSolvator import solvate_and_export

class TestMolSolvator:

    def test_solvent_construction(self):
        """
        Solvates a system with a set initial seed for the RNG and checks against a reference
        """
        args_tuple = ('./test/PDBfiles/DPPC.pdb',\
        ['./test/PDBfiles/NA+.pdb', './test/PDBfiles/CL-.pdb', './test/PDBfiles/WF.pdb', './test/PDBfiles/W.pdb'],\
        [1, 1, 9, 90],\
        './test/PDBfiles/solvator_test.pdb',\
        4, 4, 8, 4.8, 4.5, 6, 6, 0, 35, 7, False, False, False, 1000, 0)
        solvate_and_export(*args_tuple)
        refpdb = "./test/PDBfiles/solvator_ref.pdb"
        testpdb = "./test/PDBfiles/solvator_test.pdb"
        reffile = open(refpdb, 'r')
        testfile = open(testpdb, 'r')

        with open(refpdb, 'r') as reffile, open(testpdb, 'r') as testfile: 
            for ref, test in zip(reffile, testfile):
                refline = ref.strip()
                testline = test.strip()
                if (refline.split(' ') == testline.split(' ')) == False:
                    result = False
                    break
                else:
                    result = True
        assert result is True

