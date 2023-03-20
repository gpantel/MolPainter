"""
These GUI stubs keep the project logic testable separately from the interface layer.
"""

class PainterWindowStub:
    def __init__(self):
        self.mlist = MoleculeListerStub()
        self.drawarea = DrawareaStub()

    def layer_created(self, layer):
        pass

    def molecule_created(self, molecule):
        pass

    def blend_created(self, blend):
        pass


class MoleculeListerStub:
    def __init__(self):
        pass

class DrawareaStub:
    def forget_all_layers(self):
        pass

