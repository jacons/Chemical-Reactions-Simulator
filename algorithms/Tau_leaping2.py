import math
from typing import Tuple

from numpy import zeros, ndarray, vectorize, prod, ma

from numpy.random import poisson
from scipy.special import binom

from utils import StochasticAlgorithm


class TauLeaping2(StochasticAlgorithm):
    def __init__(self, initial_state: ndarray, fields: Tuple, tau: float, mol2id: dict):
        self.tau: float = tau
        self.mol2id: dict = mol2id
        self.state: ndarray = initial_state

        self.reactants: ndarray = fields[0]
        self.products: ndarray = fields[1]
        self.kinetic: ndarray = fields[2]

        self.update_weight = self.products - self.reactants

    def get_state(self):
        return self.state

    def update_molecule(self, molecule: str, qnt: int):
        self.state[self.mol2id[molecule]] += qnt
        return

    def step(self) -> float:

        # perform propensities (instantaneous rate of each reaction)
        # Oss if the molecule quantities is 0 then a_i become 0
        a = vectorize(binom)(self.state, self.reactants)
        a = prod(a, axis=1) * self.kinetic

        # if all propensities are zero means that there are no reaction to execute
        if a.sum() == 0:
            return -1

        n_reactions: ndarray = poisson(a * self.tau)

        for idx, num in enumerate(n_reactions):
            # number of reaction available for the current state
            available = vectorize(lambda x, y: x // (y if y != 0 else 1))(self.state, self.reactants[idx, :])
            num = min(min(available[self.reactants[idx, :] > 0]), num)
            self.state += num * self.update_weight[idx, :]

        return self.tau
