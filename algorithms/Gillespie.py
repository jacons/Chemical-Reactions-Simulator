import math
import random

from numpy import zeros, arange
from numpy.random import choice
from scipy.special import binom

from utils import StochasticAlgorithm


class Gillespie(StochasticAlgorithm):
    """
    Implementation with Gillespie method using the dictionaries.
    """
    def __init__(self, reactions: list, initial_state: dict):

        self.reactions: list = reactions
        self.state: dict = initial_state

    def get_state(self):
        return self.state.values()

    def update_molecule(self, molecule: str, qnt: int):
        self.state[molecule] += qnt
        return

    def step(self) -> float:
        # perform propensities (instantaneous rate of each reaction)
        n = len(self.reactions)
        a, a_0 = zeros(n), 0

        # calculate the propensities
        for idx, react in enumerate(self.reactions):
            _a = react.kinetic
            for r, l in react.reactants.items():
                # Oss if the molecule quantities is 0 then _a become 0
                _a *= binom(self.state[r], l)
            a[idx] = _a
            a_0 += _a

        # if all propensities are zero means that there are no reaction to execute
        if a_0 == 0:
            return -1

        a /= a_0
        # time event occurs
        dt = math.log(1 / random.uniform(0, 1)) / a_0

        # select the reaction to perform
        react = self.reactions[choice(arange(0, n), p=a)]

        # update the state with the reaction chosen
        for r, l in react.reactants.items():
            self.state[r] -= l
        for r, l in react.products.items():
            self.state[r] += l

        return dt  # tau
