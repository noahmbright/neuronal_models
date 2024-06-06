import neuronal_models as nm


class QIF(nm.NeuronalModel):
    def dVdt(t, V, I, E):
        # Define your ordinary differential equation here
        dVdt = E - V + I(t)
        return dVdt
