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
        self._R = R

    def calc_a_ij(self):
        """
        Calc the array a_ij for the specific mixing rule.
        
        Return:
            numpy.ndarray: a_ij array.
        """
        a = (27/64) * (self._R * self.critical_temperature)**2 / self.critical_pressure
        a_ij = np.sqrt(np.outer(a, a))
        return a_ij

    def calc_b_ij(self):
        """
        Calc the array b_ij for the specific mixing rule.
        
        Return:
            numpy.ndarray: b_ij array.
        """
        b = (27/64) * (self._R * self.critical_temperature)**2 / self.critical_pressure
        b_ij = (np.outer(b, b) / 2)
        return b_ij
    
    def calc_a(self, molar_fractions):
        """
        Calc the 'a' value using the molar cmomposition and mixing rules.

        Args:
            molar_fractions (list[float]): molar composition list.

        Return:
            float: a value.
        """
        a_ij = self.calc_a_ij()
        molar_fractions = np.array(molar_fractions)
        a = np.dot(np.dot(molar_fractions, a_ij), molar_fractions.T)
        return a

    def calc_b(self, molar_fractions):
        """
        Calc the 'b' value using the molar cmomposition and mixing rules.

        Args:
            molar_fractions (list[float]): molar composition list.

        Return:
            float: b value.
        """
        b_ij = self.calc_b_ij()
        molar_fractions = np.array(molar_fractions)
        b = np.dot(np.dot(molar_fractions, b_ij), molar_fractions.T)
        return b