import neuronal_models as nm
import numpy as np
# TODO fix E
# TODO Add other models


class QIF(nm.NeuronalModel):
    def __init__(self):
        super().__init__()
        self.E = 1

    @staticmethod
    def dALLdt(t, X, self, clamp):
        # Define your ordinary differential equation here
        if clamp:
            return np.array([0])

        return np.array([self.E - X.item() + self.I0 + self.I_inj(t)])
