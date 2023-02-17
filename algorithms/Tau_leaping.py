import math

from numpy import zeros, ndarray
from numpy.random import poisson

from utils import StochasticAlgorithm


class TauLeaping(StochasticAlgorithm):
    def __init__(self, reactions: list, initial_state: dict, tau: float):

        self.reactions: list = reactions
        self.state: dict = initial_state
        self.tau = tau

    def get_state(self):
        return self.state.values()

    def update_molecule(self, molecule: str, qnt: int):
        self.state[molecule] += qnt
        return

    def step(self) -> float:
        # perform propensities (instantaneous rate of each reaction)
        n = len(self.reactions)
        a = zeros(n)

        # calculate the propensities
        for idx, react in enumerate(self.reactions):
            _a = react.kinetic

            for r, l in react.reactants.items():
                # Oss if the molecule quantities is 0 then _a become 0
                _a *= math.comb(self.state[r], l)

            a[idx] = _a

        # if all propensities are zero means that there are no reaction to execute
        if a.sum() == 0:
            return -1

        n_reactions: ndarray = poisson(a * self.tau)

        for idx, num in enumerate(n_reactions):

            react = self.reactions[idx]

            # we fix the number of possible reaction to avoid negative numbers
            d1, d2 = react.reactants, self.state
            available = {k: int(d2[k] / d1[k]) for k in d1.keys() & d2}.values()
            # we take the minimum between: the num of reactions chosen from poisson distribution
            # and the number of max reaction available with the current state.
            num = min(min(available), num)

            for _ in range(num):
                # update the state with reaction chosen
                for r, l in react.reactants.items():
                    self.state[r] -= l

                for r, l in react.products.items():
                    self.state[r] += l

        return self.tau
