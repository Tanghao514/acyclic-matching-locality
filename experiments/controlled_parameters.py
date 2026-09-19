#!/usr/bin/env python3
"""Controlled v8 B/C measurements; original solvers remain unchanged.
Times: median of 3 sequential calls; tracemalloc: separate fourth call.
Graph construction and validation excluded. Memory is Python allocations, NOT RSS.
"""
from __future__ import annotations
import gc, hashlib, json, os, platform, random, signal, statistics, sys, time, tracemalloc
from pathlib import Path
import networkx as nx
ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/'code'
sys.path.insert(0,str(CODE))
from block_locality import block_locality
from incidence_matching import solve_incidence, subdivide
from direct_checks import analyze_matchings, verify_bad, verify_maximal
OUT=ROOT/'runs/parameters'; OUT.mkdir(parents=True,exist_ok=True)
SEED=2026091601
REPEATS=3
TIMEOUT=90

def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def serial_result(r,route):
    if route=='B':return dict(kappa=r['kappa'], stats=r['stats'])
    return dict(value=r['value'],optimal_count=str(r['optimal_count']),total_count=str(r['total_count']), branch_assignments=r['branch_assignments'])

def measured(fn):
    if not hasattr(signal, 'SIGALRM'):
        raise RuntimeError('These historical timing runners require Linux/WSL; use validation/run_checks.py on Windows.')
    times=[]; output=None
    def alarm(*_): raise TimeoutError('90s per solve time budget exceeded')
    prev=signal.signal(signal.SIGALRM,alarm)
    try:
        for _ in range(REPEATS):
            gc.collect(); signal.alarm(TIMEOUT); t=time.perf_counter()
            got=fn(); elapsed=time.perf_counter()-t; signal.alarm(0)
            times.append(elapsed); output=got
        gc.collect(); tracemalloc.start(); signal.alarm(TIMEOUT)
        got=fn(); _,peak=tracemalloc.get_traced_memory(); signal.alarm(0);tracemalloc.stop()
        return output,dict(seconds=times,median_seconds=statistics.median(times),python_peak_bytes=peak,memory_measurement='separate tracemalloc run')
    finally:
        signal.alarm(0); signal.signal(signal.SIGALRM,prev)
        if tracemalloc.is_tracing():tracemalloc.stop()

def private_theta_chain(blocks):
    G=nx.Graph()
    for j in range(blocks):
        off=10*j
        for u in (0,1):
            for v in (2,3,4):G.add_edge(off+u,off+v)
        for v in range(5):G.add_edge(off+v,off+5+v)
        if j:G.add_edge(off-9,off)
    return G

def chord_cycle(beta):
    G=nx.cycle_graph(30)
    for i in range(30):G.add_edge(i,30+i)
    chords=[(i,j) for i in range(30) for j in range(i+1,30) if j-i not in (1,29)]
    random.Random(SEED+10).shuffle(chords);G.add_edges_from(chords[:beta-1]);return G

def stretched(H,lengths):
    G=nx.Graph();G.add_nodes_from(H);nxt=max(H,default=-1)+1
    for (u,v),L in zip(H.edges(),lengths):
        path=[u]+list(range(nxt,nxt+L-1))+[v];nxt+=L-1;G.add_edges_from(zip(path,path[1:]))
    return G

def fixed_order_half(b,n=64):
    if b==0:return nx.cycle_graph(n)
    if b==2:
        G=nx.Graph();G.add_nodes_from(range(2));nxt=2;lengths=[2,2,2]
        for j in range(n-5):lengths[j%3]+=1
        for L in lengths:
            path=[0]+list(range(nxt,nxt+L-1))+[1];nxt+=L-1;G.add_edges_from(zip(path,path[1:]))
        assert len(G)==n;return G
    H={4:nx.complete_graph(4),6:nx.circular_ladder_graph(3),8:nx.cubical_graph()}[b]
    lengths=[1]*H.number_of_edges()
    for j in range(n-b):lengths[j%len(lengths)]+=1
    return stretched(H,lengths)

def incidence_costs(K,mode,seed):
    rng=random.Random(seed)
    return {(u,v):(1 if mode=='unit' else 0 if mode=='zero' else rng.randint(-5,7)) for u in sorted(K) for v in sorted(K[u])}

def validate_large_C(K,r):
    arcs=r['matching_arcs'];tails=[u for u,v in arcs]
    assert len(tails)==len(set(tails));assert len({frozenset((u,v)) for u,v in arcs})==len(arcs)
    assert all(K.has_edge(u,v) for u,v in arcs)
    F=nx.Graph();F.add_nodes_from(K);F.add_edges_from(arcs);assert not len(F) or nx.is_forest(F)
    comp={v:i for i,C in enumerate(nx.connected_components(F)) for v in C};roots=set(K)-set(tails)
    assert all(comp[r]==comp[v] for r in roots for v in K[r])
    G,lookup=subdivide(K);M=[lookup[a] for a in arcs];verts=[v for e in M for v in e]
    assert len(verts)==len(set(verts));assert not verts or nx.is_forest(G.subgraph(verts))
    return 'explicit target forest + independent base DSU/root-neighbor maximality'

def main():
    records=[];inputs=[];cases=[]
    for r in [1,4,16,64,256,600]:cases.append(('B_blocks',f'r={r}',private_theta_chain(r),{}))
    for beta in [1,2,3,4,5,6,8,10,12]:cases.append(('B_rank',f'beta={beta}',chord_cycle(beta),{}))
    for L in [1,3,10,30,100,300,1000]:
        K=stretched(nx.complete_graph(4),[L]*6);cases.append(('C_chain',f'L={L}',K,dict(mode='signed',seed=SEED+100)))
    for b in [0,2,4,6]:cases.append(('C_branches',f'b={b}',fixed_order_half(b),dict(mode='signed',seed=SEED+200)))
    for i,(family,label,G,opts) in enumerate(cases):
        route=family[0];payload=dict(family=family,label=label,nodes=sorted(G),edges=sorted(sorted(e) for e in G.edges()),options=opts)
        rec=dict(index=i,family=family,label=label,n=G.number_of_nodes(),m=G.number_of_edges())
        if route=='C':
            costs=incidence_costs(G,**opts);payload['costs']=[[u,v,w] for (u,v),w in sorted(costs.items())]
            rec.update(branches=sum(G.degree(v)==3 for v in G),target_n=G.number_of_nodes()+G.number_of_edges(),target_m=2*G.number_of_edges())
            fn=lambda G=G,costs=costs:solve_incidence(G,costs)
        else:fn=lambda G=G:block_locality(G)
        rec['input_sha256']=digest(payload);inputs.append(payload)
        print('START',i,family,label,'n',len(G),'m',G.number_of_edges(),flush=True)
        try:
            result,timing=measured(fn);rec.update(timing);rec['answer']=serial_result(result,route)
            if route=='B':
                assert result['kappa']==2 or verify_bad(G,result['witness'],result['kappa'])
                rec['validation']='returned minimal-bad witness checked by edge removal and induced forest'
            else:
                rec['validation']=validate_large_C(G,result);assert sum(costs[a] for a in result['matching_arcs'])==result['value']
                rec['optimal_count_bits']=result['optimal_count'].bit_length();rec['total_count_bits']=result['total_count'].bit_length()
            rec['status']='passed';print('DONE',round(rec['median_seconds'],6),'s peak',rec['python_peak_bytes'],flush=True)
        except TimeoutError as e:rec.update(status='timeout',error=str(e));print('TIMEOUT',flush=True)
        records.append(rec)
        (OUT/'parameter_measurements.json').write_text(json.dumps(dict(seed=SEED,repeats=REPEATS,records=records),ensure_ascii=False,indent=2))
        (OUT/'parameter_inputs.json').write_text(json.dumps(inputs,separators=(',',':')))
    env=dict(python=sys.version,platform=platform.platform(),networkx=nx.__version__,cpu=next((s.split(':',1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),'unknown'),cpu_count=os.cpu_count(),seed=SEED,repeats=REPEATS,timeout_seconds=TIMEOUT,timing='sequential, construction and validation excluded; no parallel benchmark cases',solver_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in CODE.glob('*.py')})
    (OUT/'experiment_environment.json').write_text(json.dumps(env,indent=2))
if __name__=='__main__':main()
