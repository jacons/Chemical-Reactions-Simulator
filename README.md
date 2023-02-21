# Simulator of chemical reactions with external events

---
Idea: develop (in any programming language) a stochastic simulator
of chemical reactions that takes as input:

* A set of chemical reactions

* An initial configuration (concentrations of molecules)

* A list of events (changes in concentrations), each scheduled at a precise time

Simulation of the chemical reaction should be interrupted every time a scheduled event 
has to be executed. The state has to be updated in  accordance with the event content.
Then simulation can continue

The source code should be modular: it should be easy to add a new  simulation algorithm

---

