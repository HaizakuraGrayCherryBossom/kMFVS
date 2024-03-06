The output_colex_file.py script generates .mod files for use with CPLEX. These files are saved in the folder named random_graph_cplex_mod_files.

For kMFVS problem, the filename format is:
random_graph_seed[{seed}]_nodes[{N}]_m[{m}]_original_edges[{M}]_rate[{rate}]_K[{K}]_maximum_K[{maxK}].mod
The components within the brackets are, in order, the seed used for random graph generation, the number of nodes, the parameter m used to determine the number of edges, the total number of edges M, the rate which is the percentage of edges considered for k, the value of k for the kMFVS, and maxK which indicates the maximum path length in the graph, i.e., the length of the longest path.

For the MFVS (Minimum Feedback Vertex Set) problem, the filename format is:
random_graph_seed[{seed}]_nodes[{N}]_m[{m}]_original_edges[{M}].mod
The meanings of the components are similar to those in the kMFVS filename, representing the seed, the number of nodes, the m parameter defining the number of edges, and the total number of original edges M.

This naming convention facilitates the organization and identification of files corresponding to specific instances of ILP (Integer Linear Programming) problems for solving the kMFVS and MFVS in graph theory, utilizing CPLEX.
