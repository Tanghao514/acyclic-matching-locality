"""Independent small-instance oracle: enumerate ordinary matchings in G.
No rooted forest, block-locality, path-signature, or interval-order theorem is
used to judge a candidate. The induced-forest predicate is computed by DSU.
"""
from __future__ import annotations
from functools import lru_cache
import networkx as nx


def analyze_matchings(G, collect_maximal=False):
    labels = list(G)
    pos = {v:i for i,v in enumerate(labels)}
    edges = [(pos[u],pos[v]) for u,v in G.edges()]
    emask = [(1<<u)|(1<<v) for u,v in edges]
    incident = [[] for _ in labels]
    for i,(u,v) in enumerate(edges):
        incident[u].append((v,i))
        incident[v].append((u,i))
    @lru_cache(None)
    def forest(support):
        parent=list(range(len(labels)))
        def find(x):
            while parent[x]!=x:
                parent[x]=parent[parent[x]];x=parent[x]
            return x
        for u,v in edges:
            if support>>u&1 and support>>v&1:
                a,b=find(u),find(v)
                if a==b:return False
                parent[a]=b
        return True
    stats=dict(ordinary=0,acyclic=0,maximal=0,minimal_bad=0)
    kappa=2;bad_witness=[];maximal=[]
    def rec(available,support,chosen):
        nonlocal kappa,bad_witness
        if not available:
            stats['ordinary']+=1
            if forest(support):
                stats['acyclic']+=1
                if collect_maximal:
                    if not any(not(mask&support) and forest(support|mask) for mask in emask):
                        stats['maximal']+=1
                        maximal.append(tuple(chosen))
            elif all(forest(support^emask[i]) for i in chosen):
                stats['minimal_bad']+=1
                if len(chosen)>kappa:
                    kappa=len(chosen);bad_witness=list(chosen)
            return
        bit=available&-available;u=bit.bit_length()-1
        rest=available^bit
        rec(rest,support,chosen)
        for v,i in incident[u]:
            if rest>>v&1:
                chosen.append(i)
                rec(rest^(1<<v),support|emask[i],chosen)
                chosen.pop()
    rec((1<<len(labels))-1,0,[])
    original_edges=[(labels[u],labels[v]) for u,v in edges]
    return dict(kappa=kappa, witness=[original_edges[i] for i in bad_witness],
                edges=original_edges,maximal=maximal,stats=stats)


def verify_bad(G,M,expected=None):
    if expected is not None and len(M)!=expected:return False
    verts=[v for e in M for v in e]
    if len(set(verts))!=len(verts) or not all(G.has_edge(*e) for e in M):return False
    if nx.is_forest(G.subgraph(verts)):return False
    for e in M:
        S=set(verts)-set(e)
        if S and not nx.is_forest(G.subgraph(S)):return False
    return True


def verify_maximal(G,M):
    verts=[v for e in M for v in e]
    if len(set(verts))!=len(verts) or not all(G.has_edge(*e) for e in M):return False
    S=set(verts)
    if S and not nx.is_forest(G.subgraph(S)):return False
    for u,v in G.edges():
        if u not in S and v not in S:
            if nx.is_forest(G.subgraph(S|{u,v})):return False
    return True
