"""Small, direct checks on actual matching edges. No structural theorem in oracle.
Run: python validation/run_checks.py
Outputs are separate from immutable original_data. Fixed seeds and inputs recorded.
"""
from pathlib import Path
import sys,json,hashlib,itertools,time,random,platform,statistics,math
import networkx as nx
ROOT=Path(__file__).resolve().parents[1]
CODE=ROOT/'code'
sys.path.insert(0,str(CODE))
from incidence_matching import solve_incidence
from block_locality import block_locality
OUT=ROOT/'runs/revision'; OUT.mkdir(parents=True,exist_ok=True)

def graph_data(G):
    return dict(nodes=sorted(G),edges=sorted(tuple(sorted(e)) for e in G.edges()))

def target(K):
    # Build S(K) here, independently of the solver's subdivision helper.
    n=len(K); E=sorted(tuple(sorted(e)) for e in K.edges())
    G=nx.Graph();G.add_nodes_from(range(n+len(E)));arc={}
    for j,(u,v) in enumerate(E):
        for a,b in [(u,v),(v,u)]:G.add_edge(a,n+j);arc[(a,b)]=(a,n+j)
    return G,arc

def arrays(G):
    E=sorted(tuple(sorted(e)) for e in G.edges());bits=[(1<<u)|(1<<v) for u,v in E]
    return E,bits

def induced_forest(n,E,sat):
    # Union/find scans original induced edges; no root or chain characterization.
    p=list(range(n))
    def find(x):
        while p[x]!=x:x=p[x]
        return x
    for u,v in E:
        if (sat>>u)&1 and (sat>>v)&1:
            a,b=find(u),find(v)
            if a==b:return False
            p[a]=b
    return True

def maximal_actual(n,E,bits,sat):
    return induced_forest(n,E,sat) and all(sat&b or not induced_forest(n,E,sat|b) for b in bits)

def matchings(n,E,bits):
    inc=[[] for _ in range(n)]
    for i,(u,v) in enumerate(E):inc[u].append((v,i));inc[v].append((u,i))
    def rec(available,chosen,sat):
        if not available:
            yield chosen,sat;return
        flag=available&-available;u=flag.bit_length()-1;rest=available^flag
        yield from rec(rest,chosen,sat)
        for v,i in inc[u]:
            if rest>>v&1:yield from rec(rest^(1<<v),chosen+(i,),sat|bits[i])
    yield from rec((1<<n)-1,(),0)

def oracle(G):
    E,bits=arrays(G); good=[]; ordinary=[]; bad=[];number=0
    for ids,sat in matchings(len(G),E,bits):
        number+=1
        if all(sat&b for b in bits):ordinary.append(ids)
        if induced_forest(len(G),E,sat):
            if maximal_actual(len(G),E,bits,sat):good.append(ids)
        elif all(induced_forest(len(G),E,sat^bits[i]) for i in ids):bad.append(ids)
    return dict(edges=E,maximal=good,ordinary_maximal=ordinary,bad=bad,
                ordinary_examined=number,kappa=max([2]+[len(m) for m in bad]),
                Q=sorted(set(len(m) for m in bad if len(m)>=3)))

def values(family,edge_costs):
    vals=[sum(edge_costs[i] for i in ids) for ids in family];v=min(vals)
    return dict(value=v,optimal_count=vals.count(v),total_count=len(vals))

def validate_witness(K,costs,res):
    G,arc=target(K);E,bits=arrays(G);arcs=[tuple(a) for a in res['matching_arcs']]
    assert len(arcs)==len(set(arcs))
    ends=[v for a in arcs for v in arc[a]]
    assert len(ends)==len(set(ends))
    sat=sum(1<<v for v in ends)
    assert maximal_actual(len(G),E,bits,sat)
    assert sum(costs[a] for a in arcs)==res['value']
    return True

def parent_baseline(K,costs):
    """Enumerate root/parent choices, then test actual S(K) induced graphs."""
    G,arc=target(K);E,bits=arrays(G);choices=[[None]+sorted(K[u]) for u in range(len(K))]
    best=None; opt=tot=configs=0
    for parents in itertools.product(*choices):
        configs+=1;sat=0;val=0;valid=True
        for u,v in enumerate(parents):
            if v is None:continue
            a,b=arc[(u,v)];flag=(1<<a)|(1<<b)
            if sat&flag:valid=False;break
            sat|=flag;val+=costs[(u,v)]
        if not valid or not maximal_actual(len(G),E,bits,sat):continue
        tot+=1
        if best is None or val<best:best=val;opt=1
        elif val==best:opt+=1
    return dict(value=best,optimal_count=opt,total_count=tot,configurations=configs)

def simple_cycle_lengths(K):
    # DFS all vertex-simple cycles (not cycle basis; chords allowed).
    lens=set()
    for start in K:
        def dfs(path):
            for v in K[path[-1]]:
                if v==start and len(path)>=3:lens.add(len(path))
                elif v>start and v not in path:dfs(path+[v])
        dfs([start])
    return sorted(lens)

def packings(K):
    nodes=list(K);dist=dict(nx.all_pairs_shortest_path_length(K));comps=list(nx.connected_components(K))
    out=[]
    for mask in range(1<<len(K)):
        S={v for v in nodes if mask>>v&1}
        if all(dist[u].get(v,math.inf)>=3 for u,v in itertools.combinations(S,2)):out.append(S)
    return out,[r for r in out if all(r&c for c in comps)]

def theta(length):
    K=nx.Graph();nxt=2
    for _ in range(3):
        path=[0]+list(range(nxt,nxt+length-1))+[1];nxt+=length-1
        K.add_edges_from(zip(path,path[1:]))
    return nx.convert_node_labels_to_integers(K,ordering='sorted')

def main():
    started=time.time();records=[];inputs=[]
    original=json.loads((ROOT/'data/original/small_oracle_inputs.json').read_text(encoding='utf-8-sig'))
    oldresults=json.loads((ROOT/'data/original/small_oracle_benchmarks.json').read_text())['records']
    assert len(original)==len(oldresults)==14
    # Repeat original exact answers, without replacing their historical timings.
    for inp,old in zip(original,oldresults):
        G=nx.Graph();G.add_nodes_from(inp['nodes']);G.add_edges_from(inp['edges'])
        if inp['route']=='C':
            costs={(u,v):w for u,v,w in inp['costs']};S,arc=target(G);o=oracle(S)
            inv={tuple(sorted(e)):a for a,e in arc.items()};ec=[costs[inv[e]] for e in o['edges']]
            expected=values(o['maximal'],ec);res=solve_incidence(G,costs)
            assert all(res[k]==expected[k] for k in expected)
            assert expected['value']==old['optimum'] and str(expected['optimal_count'])==old['counts']['optimal_count'] and str(expected['total_count'])==old['counts']['total_count']
            validate_witness(G,costs,res)
            records.append(dict(suite='original_recheck',case=inp['case'],status='passed',**expected))
        else:
            o=oracle(G);res=block_locality(G);assert o['kappa']==res['kappa']==old['kappa']
            if res['kappa']>2:
                W=nx.Graph();W.add_nodes_from(G);W.add_edges_from(res['witness']);ends=[v for e in res['witness'] for v in e]
                assert len(ends)==len(set(ends));sat=sum(1<<v for v in ends);E,bits=arrays(G)
                assert not induced_forest(len(G),E,sat)
                assert all(induced_forest(len(G),E,sat^((1<<u)|(1<<v))) for u,v in res['witness'])
            records.append(dict(suite='original_recheck',case=inp['case'],status='passed',kappa=o['kappa'],n=len(G),m=G.number_of_edges(),beta_max=res['stats']['beta_max']))
    print('PASS original inputs',len(records),flush=True)
    graphs=[]
    for n in range(5):
        poss=list(itertools.combinations(range(n),2))
        for mask in range(1<<len(poss)):
            K=nx.Graph();K.add_nodes_from(range(n));K.add_edges_from(e for i,e in enumerate(poss) if mask>>i&1)
            graphs.append((f'labelled_n{n}_mask{mask}',K))
    structure=[]
    for name,K in graphs:
        G,arc=target(K);o=oracle(G);packs,admissible=packings(K);rho=max(map(len,packs))
        roots=set()
        for ids in o['maximal']:
            sat=set(v for i in ids for v in o['edges'][i]);roots.add(tuple(v for v in K if v not in sat))
        assert roots=={tuple(sorted(r)) for r in admissible}
        assert min(map(len,o['maximal']))==min(map(len,o['ordinary_maximal']))==len(K)-rho
        assert sorted(set(map(len,o['maximal'])))==list(range(len(K)-rho,len(K)-nx.number_connected_components(K)+1))
        lengths=simple_cycle_lengths(K);assert o['Q']==lengths
        assert block_locality(G)['kappa']==o['kappa']==max([2]+lengths)
        structure.append(dict(case=name,**graph_data(K),rho=rho,Q=o['Q'],kappa=o['kappa'],status='passed'))
    print('PASS structural graph cases',len(structure),flush=True)
    cases=[(i['case'],nx.Graph()) for i in original if i['route']=='C']
    for (name,K),inp in zip(cases,[i for i in original if i['route']=='C']):K.add_nodes_from(inp['nodes']);K.add_edges_from(inp['edges'])
    cases.extend([('isolated',nx.empty_graph(1)),('edge',nx.path_graph(2)),('star3',nx.star_graph(3))])
    c6=None
    for idx,(name,K) in enumerate(cases):
        G,arc=target(K);o=oracle(G);inv={tuple(sorted(e)):a for a,e in arc.items()}
        for scheme in ['unit','zero','negative','signed_asymmetric']:
            seed=2026091900+idx;rng=random.Random(seed)
            costs={a:(1 if scheme=='unit' else 0 if scheme=='zero' else -1 if scheme=='negative' else rng.randint(-5,7)) for a in sorted(arc)}
            expected=values(o['maximal'],[costs[inv[e]] for e in o['edges']]);res=solve_incidence(K,costs);base=parent_baseline(K,costs)
            assert all(res[k]==base[k]==expected[k] for k in expected)
            validate_witness(K,costs,res)
            if scheme=='zero':assert res['optimal_count']==res['total_count']
            rec=dict(suite='boundary_cost',case=name,scheme=scheme,seed=seed,n=len(K),m=K.number_of_edges(),b=sum(K.degree(v)==3 for v in K),**expected,ordinary_examined=o['ordinary_examined'],witness=res['matching_arcs'],witness_checked=True,status='passed')
            records.append(rec);inputs.append(dict(case=name,scheme=scheme,seed=seed,**graph_data(K),costs=[[u,v,w] for (u,v),w in costs.items()]))
        if name=='cycle3':
            c6=dict(ordinary_sizes=sorted(map(len,o['ordinary_maximal'])),acyclic_sizes=sorted(map(len,o['maximal'])))
        print('PASS boundary',name,flush=True)
    # Modest timing baseline. Different host from the original 16 September data.
    benches=[]
    for name,K in [('K4',nx.complete_graph(4))]+[(f'theta_L{ell}',theta(ell)) for ell in (2,3,4)]:
        rng=random.Random(2026091950);costs={(u,v):rng.randint(-5,7) for u in sorted(K) for v in sorted(K[u])}
        times={};outputs={}
        for alg,fn in [('chain',solve_incidence),('parent',parent_baseline)]:
            ts=[]
            for _ in range(3):
                t=time.perf_counter();out=fn(K,costs);ts.append(time.perf_counter()-t)
            times[alg]=dict(seconds=ts,median_seconds=statistics.median(ts));outputs[alg]=out
        a,b=outputs['chain'],outputs['parent'];assert all(a[k]==b[k] for k in ('value','optimal_count','total_count'));validate_witness(K,costs,a)
        rec=dict(case=name,n=len(K),m=K.number_of_edges(),b=sum(K.degree(v)==3 for v in K),parent_configurations=b['configurations'],branch_assignments=a['branch_assignments'],times=times,result={k:a[k] for k in ('value','optimal_count','total_count')},status='passed')
        benches.append(rec);inputs.append(dict(case=name,scheme='baseline_signed',seed=2026091950,**graph_data(K),costs=[[u,v,w] for (u,v),w in costs.items()]))
        print('PASS baseline',name,'configurations',b['configurations'],'seconds',times,flush=True)
    result=dict(status='passed',date=__import__('datetime').date.today().isoformat(),elapsed_seconds=time.time()-started,environment=dict(python=sys.version,networkx=nx.__version__,platform=platform.platform(),processor=platform.processor()),records=records,structural=structure,baselines=benches,c6=c6,solver_sha256=hashlib.sha256((CODE/'incidence_matching.py').read_bytes()).hexdigest(),validation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'revision_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    (OUT/'revision_inputs.json').write_text(json.dumps(inputs,indent=2),encoding='utf-8')
    print('ALL PASSED',len(records),'answer cases;',len(structure),'structural graphs;',len(benches),'timing baselines',flush=True)
if __name__=='__main__':main()
