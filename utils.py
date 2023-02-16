import json
from typing import Tuple, Union
from numpy import array
from Reaction import Reaction
from abc import abstractmethod


class StochasticAlgorithm:
    def __init__(self):
        self.state = None

    @abstractmethod
    def step(self) -> Union[float, None]:
        pass


def json_parser(path: str) -> Tuple:  # 'sources/source1.json'
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
    events = [(event["time"], event["molecule"], event["qnt"]) for event in source["events"]]
    # Observe! we sort the events because it improves the execution (more details belows)
    events.sort()
    return reactions, initial_state, events


def reactions_parser(str_reactions: str):
    tokens = str_reactions.split()

    reactions = []

    reactants, products, steps = {}, {}, 0
    for token in tokens:
        if token == "+":
            continue
        if token == "->":
            steps = 1
            continue
        if token == ",":
            steps = 2
            continue
        if steps == 2:
            reactions.append(Reaction("", reactants, products, float(token)))
            reactants, products, steps = {}, {}, 0
            continue

        molecule = token[token.find("_") + 1:]

        if steps == 0:
            reactants[molecule] = token[:token.find("_")]
        elif steps == 1:
            products[molecule] = token[:token.find("_")]

    return reactions


def simulation(model_: StochasticAlgorithm, end_time: float, events: list):
    current_time: float = 0
    hist, times = [], []

    while current_time < end_time:

        # perform one step of computation
        dt = model_.step()

        # check if the state of molecules is not negative
        if dt is None:
            break

        # updating the current time and history
        current_time += dt
        times.append(current_time)
        hist.append(list(model_.state.values()))

        # check if the event must occur
        # the list of events is sorted by the fist element (time) and when a specific event
        # do not occur because is too early then we ignore all the rest because they have higher time values.
        for time, molecule, qnt in events:
            if time <= current_time:
                model_.state[molecule] += qnt
                events.pop(0)
            else:
                break

    return array(hist), times
