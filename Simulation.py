import matplotlib
from matplotlib import pyplot as plt
from numpy import ndarray, zeros

from algorithms.DifferentialEq import DifferentialEq
from algorithms.Gillespie import Gillespie
from algorithms.Gillespie2 import Gillespie2
from algorithms.Tau_leaping import TauLeaping
from algorithms.Tau_leaping2 import TauLeaping2
from parsers import json_parser, matrix_parser
from utils import simulation, solveOde, avg_simulations

matplotlib.use('Agg')


def getSimulation_from_path(path: str, algorithm: str, itr: int, end_time: float,
                            precision: float) -> dict:
    if algorithm == "Gillespie" or algorithm == "TauLeaping" or algorithm == "ode":

        reactions, init_state, events = json_parser(path=path)
        return getSimulationA(reactions, init_state, events, algorithm, itr, end_time, precision)

    elif algorithm == "Gillespie2" or algorithm == "TauLeaping2":

        state, fields, events, mol2id = matrix_parser(path=path)
        return getSimulationB(state, fields, events, mol2id, algorithm, itr, end_time, precision)


def getSimulationA(reactions: list, init_state: dict, events: list, algorithm: str,
                   itr: int, end_time: float, precision: float, figure: bool = True) -> dict:
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

        hists, times = solveOde(
            model_=DifferentialEq(reactions=reactions, initial_state=init_state),
            end_time=end_time,
            precision=precision)
    if figure:

        for sim, times_ in zip(hists, times):
            plt.plot(times_, sim, color="gray", alpha=0.1, linewidth=1)

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

    means = zeros(len(init_state))
    for sim in hists:
        means += sim[-1, :]
    means /= itr

    return {m: round(means[idx], 2) for idx, m in enumerate(init_state.keys())}


def getSimulationB(init_state: ndarray, fields: tuple, events: list, mol2id: dict, algorithm: str,
                   itr: int, end_time: float, precision: float, figure: bool = True) -> dict:
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

    if figure:

        for sim, times_ in zip(hists, times):
            plt.plot(times_, sim, color="gray", alpha=0.1, linewidth=1)

        x, y = avg_simulations(hists, times)
        y.sort()

        plt.plot(y, x)
        plt.legend(loc="best")
        plt.xlabel("Time")
        plt.ylabel("N° Molecules")
        plt.grid()
        plt.xlim([-end_time / 100, end_time])
        plt.savefig('static/simulation.png')
        plt.clf()

    means = zeros(len(init_state))
    for sim in hists:
        means += sim[-1, :]
    means /= itr

    return {k: round(means[v], 2) for k, v in mol2id.items()}
