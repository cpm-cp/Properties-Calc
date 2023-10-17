import numpy as np

class MixingRules:
    """Calculate the a_ij and b_ij using arrays and using the mixing rules.
    """
    def __init__(self, critical_temperature, critical_pressure, R=83.14) -> None:
        """
        Int the class: MixingRules.
        
        Args:
            critical_temperature (list[float]): critical temperature list.
            critical_pressure (list[float]): critical pressure list.
            R (float, opcional): Gas constant. Defaults to 83.14.
        """
        self.critical_temperature = np.array(critical_temperature)
        self.critical_pressure = np.array(critical_pressure)
        self.R = R

    def calc_a_ij(self):
        """
        Calc the array a_ij for the specific mixing rule.
        
        Return:
            numpy.ndarray: a_ij array.
        """
        a = (27/64) * (self.R * self.critical_temperature)**2 / self.critical_pressure
        a_ij = np.sqrt(np.outer(a, a))
        return a_ij

    def calc_b_ij(self):
        """
        Calc the array b_ij for the specific mixing rule.
        
        Return:
            numpy.ndarray: b_ij array.
        """
        b = (27/64) * (self.R * self.critical_temperature)**2 / self.critical_pressure
        b_ij = (np.outer(b, b) / 2)
        return b_ij

