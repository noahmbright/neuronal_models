import neuronal_models as nm
# TODO fix E
# TODO Add other models


class QIF(nm.NeuronalModel):
    def __init__(self):
        super().__init__()
        self.E = 1

    @staticmethod
    def dALLdt(t, V, self, clamp):
        # Define your ordinary differential equation here
        if clamp:
            return 0
        return self.E - V + self.I0 + self.I_inj(t)
