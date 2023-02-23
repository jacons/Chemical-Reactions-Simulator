import math
import random
from typing import Tuple

from numpy import ndarray, vectorize, arange
from numpy.random import choice
from scipy.special import binom

from utils import StochasticAlgorithm


class Gillespie2(StochasticAlgorithm):
    """
    Implementation with Gillespie method using matrices.
    """

    def __init__(self, initial_state: ndarray, fields: Tuple, mol2id: dict):
        self.state: ndarray = initial_state

        # since the order of molecules present into matrices is important we need a dictionary
        # that map the molecule to the right index position
        self.mol2id: dict = mol2id

        # in the following matrices the rows represent the chemical reactions,
        # the columns represent a molecules
        self.reactants: ndarray = fields[0]
        self.products: ndarray = fields[1]
        self.kinetic: ndarray = fields[2]

        # this matrix represent the number of molecules to add and subtract that we apply
        # a reaction with index i
        self.update_weight = self.products - self.reactants

    def get_state(self):
        return self.state

    def update_molecule(self, molecule: str, qnt: int):
        self.state[self.mol2id[molecule]] += qnt
        return

    def step(self) -> float:
        n = self.kinetic.shape[0]

        # perform propensities (instantaneous rate of each reaction)
        a = vectorize(binom)(self.state, self.reactants).prod(axis=1) * self.kinetic

        a_0 = a.sum()
        # if all propensities are zero means that there are
        # no reaction to execute
        if a_0 == 0:
            return -1

        a /= a_0
        # time event occurs
        dt = math.log(1 / random.uniform(0, 1)) / a_0

        # update num of molecules based on reaction occurred
        c = choice(arange(0, n), p=a)
        self.state += self.products[c, :] - self.reactants[c, :]

        return dt  # tau
