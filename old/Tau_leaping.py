import math
import random

from numpy import zeros, arange
from numpy.random import choice, poisson

from StochasticSimulator import StochasticSimulator


class TauLeaping(StochasticSimulator):
    def __init__(self, reactions: list, initial_state: dict, tau: float):
        self.reactions: list = reactions.copy()
        self.state: dict = initial_state.copy()



        self.tau = tau

    def step(self) -> float:
        # perform propensities (instantaneous rate of each reaction)
        n = len(self.reactions)
        a_0, tot = zeros(n), 0

        # calculate the propensities
        for idx, react in enumerate(self.reactions):
            _a = react.kinetic

            for r, l in react.reactants.items():
                _a *= math.comb(self.state[r], l)

            a_0[idx] = _a
            tot += _a

        n_reactions = poisson(a_0 * self.tau)

        for idx, num in enumerate(n_reactions):

            react = self.reactions[idx]
            for _ in range(num):
                # update the state with reaction chosen
                for r, l in react.reactants.items():
                    self.state[r] -= l

                for r, l in react.products.items():
                    self.state[r] += l

        """
        a_0 /= tot

        # time event occurs
        tau = math.log(1 / random.uniform(0, 1)) / tot

        # select the reaction to perform
        react = self.reactions[choice(arange(0, n), p=a_0)]

        # update the state with reaction chosen
        for r, l in react.reactants.items():
            self.state[r] -= l

        for r, l in react.products.items():
            self.state[r] += l

        return tau  # tau
        """
        return self.tau
