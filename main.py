import json

from Gillespie import Gillespie
from Reaction import Reaction

with open('sources/source1.json') as f:
    source = json.loads(f.read())

reactions = []
for reaction in source["reactions"]:
    reactants = {i["element"]: int(i["l"]) for i in reaction["reactants"]}
    products = {i["element"]: int(i["l"]) for i in reaction["products"]}

    r = Reaction(reaction["name"], reactants, products, float(reaction["kinetic"]))
    reactions.append(r)

print("Reactions:")
for r in reactions:
    print(r.show_reaction())
print()
print("State:")
initial_state = {item['element']: int(item["quantities"]) for item in source["state"]}

g = Gillespie(reactions, initial_state)

hist_a, hist_b, hist_c = [], [], []
for _ in range(5):
    g.step()
    hist_a.append(g.state["A"])
    hist_b.append(g.state["B"])
    hist_c.append(g.state["c"])

