from typing import Tuple

from numpy import ndarray, vectorize
from numpy.random import poisson
from scipy.special import binom

from utils import StochasticAlgorithm

vect_func = vectorize(lambda x, y: x // (y if y != 0 else 1))


class TauLeaping2(StochasticAlgorithm):
    """
    Implementation with Tau Leaping method using the matrices.
    """
    def __init__(self, initial_state: ndarray, fields: Tuple, tau: float, mol2id: dict):

        self.state: ndarray = initial_state
        self.tau: float = tau

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
        """
        Provide the same solution but the code is more compact and a slightly harder to understand
        :return:
        """

        # perform propensities (instantaneous rate of each reaction)
        a = vectorize(binom)(self.state, self.reactants).prod(axis=1) * self.kinetic

        # if all propensities are zero means that there are no reaction to execute
        if a.sum() == 0:
            return -1

        # predict the number of reaction that may occur in this tau-interval
        # n_reactions is an array where each element represent the number of time that the reaction i
        # "it is happened"
        n_reactions: ndarray = poisson(a * self.tau)

        for idx, num in enumerate(n_reactions):
            # We fix the number of possible reaction to avoid negative numbers
            # we take the minimum between: the num of reactions chosen from poisson distribution
            # and the number of max reaction available with the current state.
            available = vect_func(self.state, self.reactants[idx, :])
            num = min(min(available[self.reactants[idx, :] > 0]), num)
            self.state += num * self.update_weight[idx, :]

        return self.tau
