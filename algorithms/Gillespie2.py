import math
import random
from typing import Tuple

from numpy import ndarray, vectorize, prod, arange
from numpy.random import choice
from scipy.special import binom

from utils import StochasticAlgorithm


class Gillespie2(StochasticAlgorithm):
    def __init__(self, initial_state: ndarray, fields: Tuple, mol2id: dict):

        self.state: ndarray = initial_state
        self.reactants: ndarray = fields[0]
        self.products: ndarray = fields[1]
        self.kinetic: ndarray = fields[2]
        self.mol2id: dict = mol2id

        self.update_weight = self.products - self.reactants

    def get_state(self):
        return self.state

    def update_molecule(self, molecule: str, qnt: int):
        self.state[self.mol2id[molecule]] += qnt
        return

    def step(self) -> float:
        n = self.kinetic.shape[0]

        # perform propensities (instantaneous rate of each reaction)
        # Oss if the molecule quantities is 0 then a_i become 0
        a = vectorize(binom)(self.state, self.reactants)
        a = prod(a, axis=1) * self.kinetic
        a_0 = a.sum()

        # if all propensities are zero means that there are no reaction to execute
        if a_0 == 0:
            return -1

        a /= a_0

        # time event occurs
        dt = math.log(1 / random.uniform(0, 1)) / a_0

        # update num of molecules based on reaction occurred
        self.state += self.update_weight[choice(arange(0, n), p=a), :]

        return dt  # tau
