#!/usr/bin/env python3
"""Both correctness and measured small-instance direct-enumeration baselines.
Imports the frozen direct ordinary-matching enumerator, not B/C theorem-based
checking. This is a rerun with new fixed inputs/costs, not an external audit.
"""
import json,sys
from pathlib import Path
import networkx as nx
from controlled_parameters import (ROOT, OUT, measured, incidence_costs,
    block_locality, solve_incidence, subdivide, analyze_matchings,verify_bad,verify_maximal)


def direct_C(K,costs):
    G,look=subdivide(K)
    o=analyze_matchings(G,collect_maximal=True)
    inv={frozenset(e):a for a,e in look.items()}
    arcs=[inv[frozenset(e)] for e in o['edges']]
    vals=[sum(costs[arcs[i]] for i in M) for M in o['maximal']]
    v=min(vals)
    return dict(value=v,optimal_count=vals.count(v),total_count=len(vals),ordinary=o['stats']['ordinary'])

def main():
    net=nx.cycle_graph(3)
    for v in range(3):net.add_edge(v,v+3)
    theta=nx.Graph([(0,1),(0,2),(2,1),(0,3),(3,4),(4,1)])
    lol=nx.cycle_graph(4);lol.add_edge(0,4)
    barb=nx.disjoint_union(nx.cycle_graph(3),nx.cycle_graph(3));barb.add_edge(0,3)
    cases=[('path4',nx.path_graph(4)),('cycle3',nx.cycle_graph(3)),('cycle4',nx.cycle_graph(4)),
           ('theta',theta),('K4',nx.complete_graph(4)),('loop_chain',lol),('two_loops_bridge',barb),
           ('disconnected',nx.disjoint_union(nx.cycle_graph(3),nx.path_graph(2))),('empty',nx.Graph())]
    records=[];inputs=[]
    for i,(name,K) in enumerate(cases):
        costs=incidence_costs(K,'signed',2026091602+i)
        solver,t1=measured(lambda:solve_incidence(K,costs))
        oracle,t2=measured(lambda:direct_C(K,costs))
        assert all(solver[k]==oracle[k] for k in ('value','optimal_count','total_count'))
        G,look=subdivide(K);assert verify_maximal(G,[look[a] for a in solver['matching_arcs']])
        records.append(dict(route='C',case=name,n_half=len(K),n_target=len(G),m_target=G.number_of_edges(),
                            solver=t1,ordinary_enumeration=t2,counts={k:str(oracle[k]) for k in ('optimal_count','total_count')},
                            optimum=oracle['value'],ordinary_matchings=oracle['ordinary'],status='passed'))
        inputs.append(dict(route='C',case=name,nodes=list(K),edges=list(K.edges()),costs=[[u,v,w] for (u,v),w in costs.items()]))
    baseB=[('net',net),('theta_with_leaves',theta.copy()),('K4',nx.complete_graph(4)),
           ('two_nets',nx.disjoint_union(net,net)),('odd_cycle_leaf',nx.cycle_graph(7))]
    baseB[1][1].add_edges_from((v,v+5) for v in range(5))
    baseB[3][1].add_edge(0,6);baseB[4][1].add_edge(0,7)
    for name,G in baseB:
        solver,t1=measured(lambda:block_locality(G))
        oracle,t2=measured(lambda:analyze_matchings(G))
        assert solver['kappa']==oracle['kappa']
        assert solver['kappa']==2 or verify_bad(G,solver['witness'],solver['kappa'])
        records.append(dict(route='B',case=name,n=len(G),m=G.number_of_edges(),solver=t1,
                            ordinary_enumeration=t2,kappa=oracle['kappa'],ordinary_matchings=oracle['stats']['ordinary'],status='passed'))
        inputs.append(dict(route='B',case=name,nodes=list(G),edges=list(G.edges())))
    (OUT/'small_oracle_benchmarks.json').write_text(json.dumps(dict(status='passed',records=records),indent=2))
    (OUT/'small_oracle_inputs.json').write_text(json.dumps(inputs))
    print('PASS',len(records),'measured direct-oracle comparisons')
if __name__=='__main__':main()
