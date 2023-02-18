
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

from parsers import json_parser


class DifferentialEq:
    def __init__(self, reactions: list, initial_state: dict):
        self.reactions: list = reactions
        self.state: dict = initial_state

        fun = self.getLambda()

        # initial conditions
        y0 = [v for v in initial_state.values()]

        # values of time
        t = np.linspace(1, 50)

        # solving ODE
        y = odeint(fun, y0, t)

        # plot results
        plt.plot(t, y)
        plt.xlabel("Time")
        plt.ylabel("Y")
        plt.show()

    def getLambda(self):
        molecules = list(self.state.keys())
        m2a = {m: "Z[" + str(idx) + "]" for idx, m in enumerate(self.state.keys())}

        str_lambda = "lambda Z,t : ["
        for m in molecules[:-1]:
            str_lambda += self.getOde(m, m2a) + ","
        str_lambda += self.getOde(molecules[-1], m2a) + "]"

        return eval(str_lambda)

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


reactions, init_state, events = json_parser(path='../sources/source1.json')
ode = DifferentialEq(reactions, init_state)
