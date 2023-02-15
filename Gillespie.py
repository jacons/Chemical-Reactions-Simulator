import math
import random

import numpy as np
from numpy.random import choice


class Gillespie:
    def __init__(self, reactions: list, initial_state: dict):
        self.reactions: list = reactions
        self.state: dict = initial_state

    def step(self) -> float:
        # perform propensities
        a, tot = np.zeros(len(self.reactions)), 0
        for idx, react in enumerate(self.reactions):
            _a = react.kinetic
            for r, l in react.reactants.items():
                _a *= math.comb(self.state[r], l)
            a[idx] = _a
            tot += _a
        a /= tot

        reaction = choice(np.arange(0, len(a)), p=a)
        tau = np.log(1 / random.uniform(0, 1)) / tot

        for r, l in self.reactions[reaction].reactants.items():
            self.state[r] -= l

        for r, l in self.reactions[reaction].products.items():
            self.state[r] += l

        return tau
