from abc import abstractmethod
from typing import Tuple

from numpy import array, ma

colors = {0: "b", 1: "g", 2: "r", 3: "c", 4: "m", 5: "y", 6: "k"}


class StochasticAlgorithm:

    @abstractmethod
    def step(self) -> float:
        pass

    @abstractmethod
    def get_state(self):
        pass

    @abstractmethod
    def update_molecule(self, molecule: str, qnt: int):
        pass


class OdeAlgorithm:

    @abstractmethod
    def solve(self, precision: float, end_time: float):
        pass


class Reaction:
    """
    Class that simulates a chemical reaction, each instance contains the reactants,
    the products and a kinetic constant.
    """

    def __init__(self, name: str, reactants: dict, products: dict, kinetic: float):
        self.name = name

        self.reactants: dict = reactants
        self.products: dict = products

        self.kinetic: float = kinetic

    def show_reaction(self) -> str:
        reaction = ""
        for r, l in self.reactants.items():
            reaction += str(l) + str(r) + " "
        reaction += " --> "
        for r, l in self.products.items():
            reaction += str(l) + str(r) + " "

        return reaction + " (" + str(self.kinetic) + ")"


def solveOde(model_: OdeAlgorithm, end_time: float = 100, precision: float = 0.1):
    return model_.solve(precision, end_time)


def simulation(model_: StochasticAlgorithm, events: list, end_time: float = 100):
    current_time: float = 0
    hist, times = [], []

    while current_time < end_time:

        # perform one step of computation
        dt = model_.step()

        if dt == -1:
            current_time = end_time
        else:
            current_time += dt

        # updating the current time and history
        times.append(current_time)
        hist.append(list(model_.get_state()))

        # check if the event must occur
        # the list of events is sorted by the fist element (time) and when a specific event
        # do not occur because is too early then we ignore all the rest because they have higher time values.
        for time, molecule, qnt in events:
            if time <= current_time:
                model_.update_molecule(molecule, qnt)
                events.pop(0)
            else:
                break

    return array(hist), times


def avg_simulations(list_simulations: list, list_times: list) -> Tuple:
    n_molecules = list_simulations[0].shape[1]
    lens = [len(i) for i in list_simulations]

    simulations = ma.empty((max(lens), len(list_simulations), n_molecules))
    times = ma.empty((max(lens), len(list_times)))

    simulations.mask, times.mask = True, True

    for idx, (s, t) in enumerate(zip(list_simulations, list_times)):
        simulations[:len(s), idx] = s
        times[:len(t), idx] = t

    return simulations.mean(axis=1), times.mean(axis=-1)
