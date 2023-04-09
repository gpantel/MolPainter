import json
import os
import os.path
import pytest
from MolPainter.commands import Commands
from MolPainter.exporter import SystemExporter
from MolPainter.project import Project
from MolPainter.MolSolvator import solvate_and_export
from stubs import PainterWindowStub

class MockCommands(Commands):
    """
    This class represents the commands used to load and export the system.

    Some functionality is stubbed out so that the commands can be tested in isolation.
    """
    def __init__(self):
        self.project = Project()    # real project data
        self.gui = PainterWindowStub()  # GUI stub

    def fix_mol_import_path(self, newmol):
        relpath = os.path.join(os.path.dirname(__file__), newmol.path)
        if self.file_exists(relpath):
            return relpath
        return None


class TestMolPainterAndSolvator:
    """
    Tests for both MolPainter and MolSolvator
    """
    cmds = MockCommands()
    exporter = SystemExporter()
    proj_infile = os.path.join(os.path.dirname(__file__), 'tutorial.json')
    proj_outfile = os.path.join(os.path.dirname(__file__), 'output.pdb')
    proj_reffile = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'tutorial_output.pdb')
    refpdb = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'solvator_ref.pdb')
    testpdb = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'solvator_test.pdb')

    pdb_DPPC = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'DPPC.pdb')
    pdb_NA = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'NA+.pdb')
    pdb_CL = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'CL-.pdb')
    pdb_WF = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'WF.pdb')
    pdb_W = os.path.join(os.path.dirname(__file__), 'PDBfiles', 'W.pdb')

    def compare_files(self, left, right):
        with open(left, 'r') as lfile, open(right, 'r') as rfile:
            for ref, test in zip(lfile, rfile):
                refline = ref.strip()
                testline = test.strip()
                if (refline.split(' ') == testline.split(' ')) == False:
                    return False
            else:
                return True
    
    @classmethod
    def setup_class(cls):
        if os.path.exists(cls.proj_outfile):
            os.remove(cls.proj_outfile)

    @classmethod
    def teardown_class(cls):
        if os.path.exists(cls.proj_outfile):
            os.remove(cls.proj_outfile)

    def test_load_export_project(self):
        """
        End to end test that loads a project and exports the resulting system. 
        """
        with open(self.proj_infile, "r") as projfile:
            proj = json.load(projfile)
            self.cmds.load_project_from_save_data(proj)
        
        assert self.cmds.project.molecule_count == 7    # Includes "Empty"
        assert self.cmds.project.layer_count == 2

        self.exporter.export_system(self.cmds.project, self.proj_outfile, "Square")

        assert self.compare_files(self.proj_outfile, self.proj_reffile) is True
        
    def test_solvent_construction(self):
        """
        Solvates a system with a set initial seed for the RNG and checks against a reference
        """
        args_tuple = (self.pdb_DPPC,\
        [self.pdb_NA, self.pdb_CL, self.pdb_WF, self.pdb_W],\
        [1, 1, 9, 90],\
        self.testpdb,\
        4, 4, 8, 4.8, 4.5, 6, 6, 0, 35, 7, False, False, False, 1000, 0)
        solvate_and_export(*args_tuple)

        assert self.compare_files(self.refpdb, self.testpdb) is True
