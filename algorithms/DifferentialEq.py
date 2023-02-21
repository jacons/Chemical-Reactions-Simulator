import numpy as np
from scipy.integrate import odeint

from utils import OdeAlgorithm


class DifferentialEq(OdeAlgorithm):
    """
    Implementation with Differential Equation. This implementation doesn't handle the events and iterations
    since is not a stochastic approach. The implementation takes as input a list of chemical reactions and
    an initial state (dictionary). Then it builds a function to pass to "odeint" function. Provide an exact solution.
    """

    def __init__(self, reactions: list, initial_state: dict):
        self.reactions: list = reactions
        self.initial_state: dict = initial_state

        self.fun = self.getLambda()

        return

    def solve(self, precision: float, end_time: float):

        time = np.linspace(0, end_time, num=int(end_time / precision))
        result = odeint(self.fun,
                        [v for v in self.initial_state.values()],
                        time)

        return [result], [time]

    def getLambda(self):
        """
        Builds from dictionary the lambda function. It creates a "string version of lambda" then
        it transform the string into a real function by "eval"
        :return:
        """
        molecules = list(self.initial_state.keys())
        m2a = {m: "Z[" + str(idx) + "]" for idx, m in enumerate(self.initial_state.keys())}

        str_lambda = "lambda Z,t : ["
        for m in molecules[:-1]:
            str_lambda += self.getOde(m, m2a) + ","  # foreach molecules concatenate differential equation
        str_lambda += self.getOde(molecules[-1], m2a) + "]"

        # str_lambda is a string
        return eval(str_lambda)  # real lambda function

    def getOde(self, molecule: str, m2a: dict) -> str:
        ode_ = ""
        for r in self.reactions:
            if molecule in r.reactants:
                ode_ += "-" + str(r.reactants[molecule]) + " * " + str(r.kinetic)
                for (k, v) in r.reactants.items():
                    ode_ += " * " + "np.power(" + m2a[k] + "," + str(v) + ")"

            if molecule in r.products:
                ode_ += " + " + str(r.products[molecule]) + " * " + str(r.kinetic)
                for (k, v) in r.reactants.items():
                    ode_ += " * " + "np.power(" + m2a[k] + "," + str(v) + ")"
        return ode_
