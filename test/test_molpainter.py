import json
import os.path
import pytest
from MolPainter.project import Project
from MolPainter.commands import Commands
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


class TestMolPainter:
    """
    This represents an end to end test that loads a project and exports the resulting system. 
    """
    cmds = MockCommands()

    def test_load_export_project(self):
        path = os.path.join(os.path.dirname(__file__), 'tutorial.json')
        with open(path, "r") as projfile:
            proj = json.load(projfile)
            self.cmds.load_project_from_save_data(proj)

        assert self.cmds.project.molecule_count == 7    # Includes "Empty"
        assert self.cmds.project.layer_count == 2

        # TODO: Test the exporter
        