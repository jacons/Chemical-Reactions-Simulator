import json
from typing import Tuple

from numpy import zeros

from utils import Reaction


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


def matrix_parser(path: str):
    with open(path) as f:
        source = json.loads(f.read())

    if source is None:
        return None, None, None

    n_reactions = len(source["reactions"])
    n_molecules = len(source["initial_state"])

    # mapping molecules
    mol2id = {item["molecule"]: idx for idx, item in enumerate(source["initial_state"])}
    id2mol = {idx: item["molecule"] for idx, item in enumerate(source["initial_state"])}

    # state
    state = zeros(n_molecules)
    for item in source["initial_state"]:
        state[mol2id[item["molecule"]]] = item["qnt"]

    # reactants, products and kinetics
    reactants = zeros((n_reactions, n_molecules))
    products = zeros((n_reactions, n_molecules))
    kinetic = zeros(n_reactions)

    for idx, reaction in enumerate(source["reactions"]):

        for item in reaction["reactants"]:
            reactants[idx, mol2id[item["molecule"]]] = item["l"]

        for item in reaction["products"]:
            products[idx, mol2id[item["molecule"]]] = item["l"]

        kinetic[idx] = reaction["kinetic"]

    # Parse events
    events = [(event["time"], event["molecule"], event["qnt"]) for event in source["events"]]
    # Observe! we sort the events because it improves the execution (more details belows)
    events.sort()

    return state, (reactants, products, kinetic), events, mol2id
