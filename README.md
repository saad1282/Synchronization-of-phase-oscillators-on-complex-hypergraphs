# Synchronization-of-phase-oscillators-on-complex-hypergraphs

This file contains codes and necessary files for simulations and theory in the work "Synchronization of phase oscillators on complex hypergraphs."

Generate_edge_list.m: Matlab code with a function to generate list of hyperedges. Required parameters are N - hypergraph size and hyperedge_size - size of hyperedge. Degree distribution should be defined.

Edgesandtriangles.m: Matlab code that generates text files with a list of links, links degree, triangles, triangle degree and natural frequency, using the generate_edge_list function.

edgelist_5000.txt: list of links (hyperedge of size 2)

trianglelist_5000.txt: list of triangles (hyperedge of size 3)

k1_degree_5000.txt: links degree list

k2_degree_5000.txt: triangles degree list

Kuramoto_simulation.m: Matlab code for numerical simulation of Kuramoto model. Calculates the order parameters.

Theory for order parameter.nb: Mathematica code to calculate order parameter predicted by mean field theory.
