import json
from typing import Tuple

from numpy import array

from Reaction import Reaction


def parser(path: str) -> Tuple:  # 'sources/source1.json'
    """
    Takes the path of json file and parse it into:
    1) list of reactions, list of events and initial state ( number of molecules for each element involved)
    :param path: path of json file
    :return:  list of reactions, list of events and initial state
    """
    with open(path) as f:
        source = json.loads(f.read())

    if source is None:
        return None, None, None

    reactions: list = []

    # Parse reactions
    for reaction in source["reactions"]:
        # build a dict where key is the name of molecule and v n° molecules involved
        reactants = {i["molecule"]: int(i["l"]) for i in reaction["reactants"]}
        products = {i["molecule"]: int(i["l"]) for i in reaction["products"]}

        r = Reaction(reaction["name"], reactants, products, float(reaction["kinetic"]))
        reactions.append(r)

    # Parse initial state
    initial_state = {item['molecule']: int(item["qnt"]) for item in source["initial_state"]}

    # Parse events
    # to do

    return reactions, initial_state, None


def single_simulation(model_, end_time):
    current_time: float = 0
    hist, times = [], []

    while current_time < end_time:
        dt = model_.step()
        current_time += dt
        times.append(current_time)
        hist.append(list(model_.state.values()))
    return array(hist), times
