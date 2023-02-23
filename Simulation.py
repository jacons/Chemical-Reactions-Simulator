import matplotlib
from matplotlib import pyplot as plt
from numpy import ndarray, zeros

from algorithms.DifferentialEq import DifferentialEq
from algorithms.Gillespie import Gillespie
from algorithms.Gillespie2 import Gillespie2
from algorithms.Tau_leaping import TauLeaping
from algorithms.Tau_leaping2 import TauLeaping2
from parsers import dicts_parser, matrix_parser
from utils import simulation, solveOde, avg_simulations

matplotlib.use('Agg')


def getSimulation_from_path(path: str, algorithm: str, itr: int, end_time: float,
                            precision: float) -> dict:
    """
    Universal simulator from json path
    :param path: json source
    :param algorithm: select which kind of algorithm to execute
    :param itr: Number of simulations
    :param end_time: Duration of simulation
    :param precision: Precision (TauLeaping or Ode)
    :return: dictionary of (mean) final state of simulations
    """

    # I have to separate with two case because some algorithms require different input format
    if algorithm == "Gillespie" or algorithm == "TauLeaping" or algorithm == "ode":

        # Dictionary based algorithms
        reactions, init_state, events = dicts_parser(path=path)
        return getSimulationA(reactions, init_state, events, algorithm, itr, end_time, precision)

    elif algorithm == "Gillespie2" or algorithm == "TauLeaping2":

        # Matrix based algorithms
        state, fields, events, mol2id = matrix_parser(path=path)
        return getSimulationB(state, fields, events, mol2id, algorithm, itr, end_time, precision)


def getSimulationA(reactions: list, init_state: dict, events: list, algorithm: str,
                   itr: int, end_time: float, precision: float = 1, figure: bool = True) -> dict:
    """
    Simulate a dictionary based algorthm.
    :param reactions: list of reactions
    :param init_state: dict of initial state
    :param events: list of events(tuple)
    :param algorithm: name of algorithm
    :param itr: number of simulations
    :param end_time: duration of simulation
    :param precision: precision (only for Ode)
    :param figure: True if we want save the plot
    :return: dict of final state
    """
    hists, times = [], []

    if algorithm == "Gillespie":
        for _ in range(itr):
            simulation_, times_ = simulation(
                model_=Gillespie(reactions=reactions.copy(), initial_state=init_state.copy()),
                end_time=end_time,
                events=events.copy())

            hists.append(simulation_)
            times.append(times_)

    elif algorithm == "TauLeaping":
        for _ in range(itr):
            simulation_, times_ = simulation(
                model_=TauLeaping(reactions=reactions.copy(), initial_state=init_state.copy(), tau=precision),
                end_time=end_time,
                events=events.copy())

            hists.append(simulation_)
            times.append(times_)

    elif algorithm == "ode":
        itr = 1  # it does not make sense more than one iteration
        hists, times = solveOde(
            model_=DifferentialEq(reactions=reactions.copy(), initial_state=init_state.copy()),
            end_time=end_time,
            precision=precision)

    if figure:

        # plot (in background all simulations)
        for sim, times_ in zip(hists, times):
            plt.plot(times_, sim, color="gray", alpha=0.1, linewidth=1)

        # try to build a generalized pattern from all simulation (general behavior)
        x, y = avg_simulations(hists, times)
        y.sort()

        labels = [m.replace("_", "") for m in init_state.keys()]
        plt.plot(y, x, label=labels)
        plt.legend(loc="best")
        plt.xlabel("Time")
        plt.ylabel("N° Molecules")
        plt.grid()
        plt.xlim([-end_time / 100, end_time])
        plt.savefig('static/simulation.png')
        plt.clf()

    # to have a more reliable final state , we take the mean of all simulation performed
    means = zeros(len(init_state))
    for sim in hists:
        means += sim[-1, :]
    means /= itr

    return {m: round(means[idx], 2) for idx, m in enumerate(init_state.keys())}


def getSimulationB(init_state: ndarray, fields: tuple, events: list, mol2id: dict, algorithm: str,
                   itr: int, end_time: float, precision: float) -> dict:
    """
    Simulate a dictionary based algorthm.
    :param init_state: dict of initial state
    :param fields: tuple (reactants,products,kinetic) matrices
    :param events: list of events(tuple)
    :param algorithm: name of algorithm
    :param itr: number of simulations
    :param end_time: duration of simulation
    :param precision: precision (only for Tau leaping)
    :param mol2id: dictionary that map the molecule to an index
    :return: dict of final state
    """
    hists, times = [], []
    if algorithm == "Gillespie2":

        for _ in range(itr):
            simulation_, times_ = simulation(
                model_=Gillespie2(initial_state=init_state.copy(), fields=fields, mol2id=mol2id),
                end_time=end_time,
                events=events.copy())

            hists.append(simulation_)
            times.append(times_)

    elif algorithm == "TauLeaping2":

        for _ in range(itr):
            simulation_, times_ = simulation(
                model_=TauLeaping2(initial_state=init_state.copy(), fields=fields, mol2id=mol2id, tau=precision),
                end_time=end_time,
                events=events.copy())

            hists.append(simulation_)
            times.append(times_)

    # to have a more reliable final state , we take the mean of all simulation performed
    means = zeros(len(init_state))
    for sim in hists:
        means += sim[-1, :]
    means /= itr

    return {k: round(means[v], 2) for k, v in mol2id.items()}
