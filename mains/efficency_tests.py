from Simulation import getSimulationA, getSimulationB
from parsers import json_parser, matrix_parser


def run_test(path='../sources/source1.json'):
    itr = 1
    end_time = 50
    precision = 0.01

    reactions, init_state, events = json_parser(path=path)
    getSimulationA(reactions, init_state, events, "Gillespie", itr, end_time, precision, figure=False)
    getSimulationA(reactions, init_state, events, "TauLeaping", itr, end_time, precision, figure=False)
    getSimulationA(reactions, init_state, events, "ode", itr, end_time, precision, figure=False)

    state, fields, events, mol2id = matrix_parser(path=path)
    getSimulationB(state, fields, events, mol2id, "Gillespie2", itr, end_time, precision, figure=False)
    getSimulationB(state, fields, events, mol2id, "TauLeaping2", itr, end_time, precision, figure=False)
    print("ok")


run_test()
