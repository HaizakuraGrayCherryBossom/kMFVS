import networkx as nx
import time

class k_distance_limited_MFVS():
    def _detect_cycles_with_time_limit(nodes_with_self_loop,G,N_cyclese):
        tmp = [nodes_with_self_loop]
        n_gon = 1
        while n_gon == 1 or n_gon != N_cyclese:
            n_gon += 1
            exec('all_the_nodes = G.nodes();'\
            +'ll = [G];tmp.append(list(set([tuple(sorted(['+\
            ','.join(['node{}'.format(j) for j in range(1,n_gon+1)])+'])) '\
            +'for graph in ll for node1 in all_the_nodes '\
            + ''.join(['for node{} in [node{} for node{} '.format(k,k,k)\
            +'in set(graph.successors(node{})) if '.format(k-1)+' and '\
            .join(['node{}!=node{}'.format(k,j) for j in range(1,k)])+'] ' \
            if k!=n_gon else 'for node{} in [node{} for node{} '.format(k,k,k)\
            +'in set(graph.successors(node{})) if '.format(k-1)+' and '\
            .join(['node{}!=node{}'.format(k,j) for j in range(1,k)]) \
            for k in range(2,n_gon+1)])\
            +' and (node1 in set(graph.successors(node{})))] if '.format(n_gon)\
            +'all([1 if not set(n_tuple)<set(['\
            +','.join(['node{}'.format(j) for j in range(1,n_gon+1)])\
            +']) else 0 for n_tuple in {}])'.format(sum(tmp,[]))\
            +'])))',{'G':G,'tmp':tmp}\
            )
        return [\
        ' + '.join(['x{} '.format(k) for k in j]) + '>= 1\n' \
        for n_gon in tmp[1:] for j in n_gon\
        ]

    #G: graphs (dict) <ex> {id1:g1,id2:g2,...}, filename: file name (str)
    #Perform ILP formulation from GRAPH data. Then CPLEX file output.
    def OutputCPLEXFile(G,K,filename,N_cyclese=0):
        nodes = G.nodes()
        N = G.number_of_nodes()
        if K >= N:
            print('K,N',K,N)
            print('Warning: If K >= N, the result is the same as K=N-1. In other words, the result is the same as the usual minimum FVS.')
        display = 'display solution variables -'
        obj = ''.join([' x{} +'.format(i) for i in nodes])[1:-2]
        x = obj.replace('+ ','')
        sbj_list = []

        #Subject to: Self loop
        self_loop_edge_candidates = {(i,i) for i in nodes}
        nodes_with_self_loop = list(
        {edge[0] for edge in self_loop_edge_candidates if edge in G.edges()}
        )
        sbj_list += [\
        'x{0} <= 1\n'.format(sl)+'x{0} >= 1\n'.format(sl) \
        for sl in nodes_with_self_loop\
        ]

        #Subject to: Cycle detection
        cycle_detection_time = 0.00
        if N_cyclese != 0 and N_cyclese != 1:
            tmp = k_distance_limited_MFVS._detect_cycles_with_time_limit(nodes_with_self_loop,G,N_cyclese)
            sbj_list += tmp

        #Subject to:conditional formula for l_parameters
        sbj_list += sum(\
        [['l({0}) - l({1}) + {2} x{0} >= 1\n'.format(edge[0],edge[1],N)] for edge in G.edges() if edge[0]!=edge[1]]\
        ,[]\
        )

        #Stringing conditionals
        sbj = ''.join(list(set(sbj_list)))

        #l parameters: CPLEX;generary
        l_parameters = [
        ['l({}) '.format(node),'0 <= l({}) <= {}\n'.format(node,K)]
        for node in nodes if 'l({}) '.format(node) in sbj
        ]
        l, l_b = ''.join([l[0] for l in l_parameters]), ''.join([l[1] for l in l_parameters])

        filename = filename[:-4] + '.mod'
        #Write ILP formula to file and output
        with open(filename,'w') as f:
            f.write('enter {}'.format(filename.split('/')[-1])[:-4]+'\n')#[:-4]:Remove the extension
            f.write("minimize\n")
            f.write(obj+'\n')
            f.write('subject to\n')
            f.write(sbj)
            f.write('bounds\n')
            f.write(l_b)
            f.write('generary\n')
            f.write(l+'\n')
            f.write('binary\n')
            f.write(x+'\n')
            f.write('end\n')
            f.write('optimize\n')
            f.write(display+'\n')
            f.write('quit')

    def max_path_length_per_node(graph):
        Max_path_num = []
        for subgraph_nodes in list(nx.weakly_connected_components(graph)):
            subgraph = graph.subgraph(subgraph_nodes)
            max_path_lengths = {}
            for node in subgraph.nodes:
                paths = nx.single_source_dijkstra_path_length(subgraph,node)
                max_path_lengths[node] = max(paths.values())
            Max_path_num.append(max(max_path_lengths.values()))
        return max(Max_path_num)
