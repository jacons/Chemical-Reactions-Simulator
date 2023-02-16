import math
from typing import Union

from numpy import zeros
from numpy.random import poisson

from utils import StochasticAlgorithm


class TauLeaping(StochasticAlgorithm):
    def __init__(self, reactions: list, initial_state: dict, tau: float):
        super().__init__()

        self.reactions: list = reactions
        self.state: dict = initial_state

        self.tau = tau

    def step(self) -> Union[float, None]:
        # perform propensities (instantaneous rate of each reaction)
        n = len(self.reactions)
        a_0 = zeros(n)

        for v in self.state.values():
            if v < 0:
                return None
        # calculate the propensities
        for idx, react in enumerate(self.reactions):
            _a = react.kinetic

            for r, l in react.reactants.items():
                _a *= math.comb(self.state[r], l)

            a_0[idx] = _a

        if a_0.sum() == 0:
            return None

        n_reactions = poisson(a_0 * self.tau)

        for idx, num in enumerate(n_reactions):

            react = self.reactions[idx]
            for _ in range(num):
                # update the state with reaction chosen
                for r, l in react.reactants.items():
                    self.state[r] -= l

                for r, l in react.products.items():
                    self.state[r] += l

        return self.tau
