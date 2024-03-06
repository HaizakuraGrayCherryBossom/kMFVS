from k_distance_limited_MFVS import k_distance_limited_MFVS as kMFVS
import matplotlib.pyplot as plt
import networkx as nx
import random
import os

num_of_trials = 15 #not 10 but 15
rand_num = 7000000
rand = random.sample(range(rand_num),rand_num)
dirname_for_random_graph = './random_graph_cplex_mod_files'
if not os.path.exists(dirname_for_random_graph):
    os.makedirs(dirname_for_random_graph)
node_sizes = [100,200,300,400]
rates = [[2,3,4,5],[2,3,4],[2,3,4],[2,3,4]]
m_list = [[2,3,4,5],[2,2.125,2.25],[1.8,2,2.05],[1.2,1.4,1.6]]

node_sizes = [100,200]
rates = [[2,3,4,5]+list(range(10,96,5)),[2,3,4,5]+list(range(10,96,5))]
m_list = [[2,4],[2.25]]

for i, N in enumerate(node_sizes):
    for m in m_list[i]:
        for rate in rates[i]:
            K = int(rate*0.01*N)
            M = int(m*float(N))
            for trial in range(num_of_trials):
                seed = random.sample(rand,1)[0]
                rand.remove(seed)
                dirname_node_size = dirname_for_random_graph + '/node_size[{}]'.format(N)
                if not os.path.isdir(dirname_node_size):
                    os.makedirs(dirname_node_size)

                g = nx.gnm_random_graph(N,M,directed=True,seed=seed)

                maxK = kMFVS.max_path_length_per_node(g)

                filename_kMFVS_without_cycle_detection = dirname_node_size + '/' + \
                f'random_graph_seed[{seed}]_nodes[{N}]_m[{m}]_original_edges[{M}]_rate[{rate}]_K[{K}]_maximum_K[{maxK}].mod'
                file_normal_MFVS_without_cycle_detection = dirname_node_size + '/' + \
                f'random_graph_seed[{seed}]_nodes[{N}]_m[{m}]_original_edges[{M}].mod'

                filename_K_MFVS_with_cycle_detection = filename_kMFVS_without_cycle_detection
                file_normal_MFVS_with_cycle_detection = file_normal_MFVS_without_cycle_detection

                kMFVS.OutputCPLEXFile(g,K,filename_K_MFVS_with_cycle_detection,N_cyclese=6)
                kMFVS.OutputCPLEXFile(g,g.number_of_nodes()-1,file_normal_MFVS_with_cycle_detection,N_cyclese=6)
