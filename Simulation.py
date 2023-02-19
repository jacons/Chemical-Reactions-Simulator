from matplotlib import pyplot as plt
from numpy import ndarray

from algorithms.DifferentialEq import DifferentialEq
from algorithms.Gillespie import Gillespie
from algorithms.Gillespie2 import Gillespie2
from algorithms.Tau_leaping import TauLeaping
from algorithms.Tau_leaping2 import TauLeaping2
from parsers import json_parser, matrix_parser
from utils import simulation, solveOde, colors


def getSimulation_from_path(path: str, algorithm: str, itr: int, end_time: float,
                            precision: float):
    if algorithm == "Gillespie" or algorithm == "TauLeaping" or algorithm == "ode":

        reactions, init_state, events = json_parser(path=path)
        getSimulationA(reactions, init_state, events, algorithm, itr, end_time, precision)
        return

    elif algorithm == "Gillespie2" or algorithm == "TauLeaping2":

        state, fields, events, mol2id = matrix_parser(path=path)
        getSimulationB(state, fields, events, mol2id, algorithm, itr, end_time, precision)
        return


def getSimulationA(reactions: list, init_state: dict, events: list, algorithm: str,
                   itr: int, end_time: float, precision: float):
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

    for sim, times_ in zip(hists, times):

        for idx, m in enumerate(init_state.keys()):
            plt.plot(times_, sim[:, idx], color=colors[idx], alpha=1, linewidth=0.7)

    plt.savefig('static/simulation.png')
    return


def getSimulationB(state: ndarray, fields: tuple, events: list, mol2id: dict, algorithm: str,
                   itr: int, end_time: float, precision: float):
    hists, times = [], []
    if algorithm == "Gillespie2":

        for _ in range(itr):
            simulation_, times_ = simulation(
                model_=Gillespie2(initial_state=state.copy(), fields=fields, mol2id=mol2id),
                end_time=end_time,
                events=events.copy())

            hists.append(simulation_)
            times.append(times_)

    elif algorithm == "TauLeaping2":

        for _ in range(itr):
            simulation_, times_ = simulation(
                model_=TauLeaping2(initial_state=state.copy(), fields=fields, mol2id=mol2id, tau=precision),
                end_time=end_time,
                events=events.copy())

            hists.append(simulation_)
            times.append(times_)

    for sim, times_ in zip(hists, times):

        for idx in range(state.shape[0]):  # plot chart for each molecule
            plt.plot(times_, sim[:, idx], color=colors[idx], alpha=1, linewidth=0.7)

    plt.savefig('/cache/simulation.png')
    return
