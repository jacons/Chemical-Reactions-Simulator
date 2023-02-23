import json
from typing import Tuple

from numpy import zeros

from utils import Reaction


def dicts_parser(path: str) -> Tuple:  # 'sources/source1.json'
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


def reactions_parser(str_reactions: str) -> list:
    """
    Used mostly in the web graphical interface, takes inputs as
    1_A + 2_B -> 3_C , 0.4 ;
    1_D + 2_E -> 3_C + 5_F , 0.4 ;
    where the numbers before the molecules identify the number of molecules involved into reaction.

    The method returns a list of reactions (list of Reaction class)
    """

    # tokenize each word by space
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
            reactants[molecule] = int(token[:token.find("_")])
        elif steps == 1:
            products[molecule] = int(token[:token.find("_")])

    return reactions


def state_parser(str_state: str) -> dict:
    """
    Used mostly in the web graphical interface, takes inputs as
    A : 10 ;
    B : 50 ; ecc,,
    where the numbers before the molecules identify the number of molecules involved into reaction.

    The method returns a dictionary ( molecule -> quantities )
    """
    flag, molecule, result = True, None, {}
    for token in str_state.split():
        if token == ";":
            continue
        if token == ":":
            flag = False if flag else False
        else:
            if flag:
                molecule = token
            else:
                result[molecule] = int(token)
                flag = True
    return result


def events_parser(str_state: str) -> list:
    """
    Used mostly in the web graphical interface, takes inputs as
    0.1 A 4 ;
    0.50 B -56 ;
    Each line identify an event where : the fist number is the time when the event (MORE OR LESS)
    could happen , the letter identify the molecule involved and the last values the quantity to change

    The method returns a list of tuple sorted by time [ (time,mol,qnt),(time,mol,qnt) ...]
    """
    flag, time, molecule, result = 0, 0, None, []
    for token in str_state.split():
        if token == ";":
            continue
        elif flag == 0:
            time = token
            flag = 1
        elif flag == 1:
            molecule = token
            flag = 2
        elif flag == 2:
            result.append((float(time), molecule, int(token)))
            flag = 0
    result.sort()
    return result


def matrix_parser(path: str):
    """
    Takes the path of json file and parse it into:
        - initial state, array(# molecules)
        - kinetic , array(# reaction)
        - reactants , matrix (# reaction, # molecules)
        - products , matrix (# reaction, # molecules)
        - events, list of tuple (as before)
        - dictionary that map the molecule to index (of matrices)

    """
    with open(path) as f:
        source = json.loads(f.read())

    if source is None:
        return None, None, None

    n_reactions = len(source["reactions"])
    n_molecules = len(source["initial_state"])

    # mapping molecules
    mol2id = {item["molecule"]: idx for idx, item in enumerate(source["initial_state"])}

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
