import json

from Reaction import Reaction


def parser(path: str):  # 'sources/source1.json'

    with open(path) as f:
        source = json.loads(f.read())

    if source is None:
        return None

    reactions = []
    for reaction in source["reactions"]:

        reactants = {i["molecule"]: int(i["l"]) for i in reaction["reactants"]}
        products = {i["molecule"]: int(i["l"]) for i in reaction["products"]}

        r = Reaction(reaction["name"], reactants, products, float(reaction["kinetic"]))
        reactions.append(r)

    initial_state = {item['molecule']: int(item["qnt"]) for item in source["initial_state"]}

    return reactions, initial_state
