import math
import random

from numpy import zeros, arange
from numpy.random import choice

from StochasticSimulator import StochasticSimulator


class Gillespie(StochasticSimulator):
    def __init__(self, reactions: list, initial_state: dict):
        self.reactions: list = reactions
        self.state: dict = initial_state

    def step(self) -> float:
        # perform propensities
        n = len(self.reactions)
        a, tot = zeros(n), 0

        for idx, react in enumerate(self.reactions):
            _a = react.kinetic

            for r, l in react.reactants.items():
                _a *= math.comb(self.state[r], l)

            a[idx] = _a
            tot += _a

        a /= tot

        react = self.reactions[choice(arange(0, n), p=a)]

        for r, l in react.reactants.items():
            self.state[r] -= l

        for r, l in react.products.items():
            self.state[r] += l

        return math.log(1 / random.uniform(0, 1)) / tot  # tau
