"""Exact locality order on arbitrary graphs, parameterized by cycle rank.

No bipartiteness or degree restriction.  The independent oracle in tests does
not use the cycle+private-leaf characterization employed here.
"""
from __future__ import annotations
import networkx as nx


def best_cycle_tiling(cycle,available):
    """Maximize permitted monomers, cover other cycle vertices by dimers.
    Return (monomer_count, monomer_vertices, dimer_edges), or None.
    Each monomer stands for a selected private-leaf matching edge.
    """
    l=len(cycle)
    def linear(lo,hi):
        length=hi-lo;dp=[None]*(length+1);parent=[None]*(length+1);dp[0]=0
        for i in range(length):
            if dp[i] is None:continue
            if cycle[lo+i] in available:
                cand=dp[i]+1
                if dp[i+1] is None or cand>dp[i+1]:dp[i+1]=cand;parent[i+1]=(i,'m')
            if i+2<=length:
                cand=dp[i]
                if dp[i+2] is None or cand>dp[i+2]:dp[i+2]=cand;parent[i+2]=(i,'d')
        if dp[length] is None:return None
        marks=[];dimers=[];j=length
        while j:
            i,kind=parent[j]
            if kind=='m':marks.append(cycle[lo+i])
            else:dimers.append((cycle[lo+i],cycle[lo+i+1]))
            j=i
        return dp[length],marks,dimers
    ans=linear(0,l)
    inside=linear(1,l-1)
    if inside is not None:
        cand=(inside[0],inside[1],inside[2]+[(cycle[-1],cycle[0])])
        if ans is None or cand[0]>ans[0]:ans=cand
    return ans


def parity_run_tiling(cycle,available):
    """Closed form for the maximum monomer count on a permitted cycle.

    Let r be the number of parity runs in the increasing list of available
    positions.  The best count is the largest t<=r with t congruent to |C| mod 2.
    A negative value means no perfect tiling.  Returns a witness, not just t.
    """
    l=len(cycle);selected=[];last=None
    for i,u in enumerate(cycle):
        if u in available and i%2!=last:selected.append(i);last=i%2
    if (len(selected)-l)%2:
        if not selected:return None
        selected.pop()
    if not selected:
        return (0,[],[(cycle[i],cycle[i+1]) for i in range(0,l,2)]) if l%2==0 else None
    dimers=[]
    for i,j in zip(selected,selected[1:]+selected[:1]):
        p=(i+1)%l
        while p!=j:
            nxt=(p+1)%l;assert nxt!=j
            dimers.append((cycle[p],cycle[nxt]));p=(nxt+1)%l
    return len(selected),[cycle[i] for i in selected],dimers


def _ordered_cycle(edges):
    adj={}
    for u,v in edges:adj.setdefault(u,[]).append(v);adj.setdefault(v,[]).append(u)
    if len(adj)<3 or any(len(v)!=2 for v in adj.values()):return None
    start=next(iter(adj));out=[start];seen={start};prev=start;u=adj[start][0]
    while u!=start:
        if u in seen:return None
        out.append(u);seen.add(u);vs=adj[u];v=vs[0] if vs[0]!=prev else vs[1];prev,u=u,v
    return out if len(out)==len(adj) else None


def cactus_locality(G):
    """O(n+m) structural algorithm for arbitrary cactus graphs.
    Raises ValueError outside the class; counts are not returned.
    """
    if G.is_multigraph() or nx.number_of_selfloops(G):raise ValueError('Simple graph required')
    blocks=[]
    for edges_it in nx.biconnected_component_edges(G):
        edges=list(edges_it)
        if len(edges)==1:continue
        cycle=_ordered_cycle(edges)
        if cycle is None:raise ValueError('Not a cactus graph')
        blocks.append(cycle)
    best=2;witness=[];cycle_witness=None;tilings=0
    for cycle in blocks:
        C=set(cycle);available={u for u in cycle if G.degree(u)>2}
        til=parity_run_tiling(cycle,available);tilings+=1
        if til is None:continue
        q=(len(cycle)+til[0])//2
        if q>best:
            leaves={u:next(v for v in G[u] if v not in C) for u in til[1]}
            witness=til[2]+[(u,leaves[u]) for u in til[1]];best=q;cycle_witness=cycle
    return dict(kappa=best,witness=witness,cycle=cycle_witness,
                stats=dict(method='cactus',cycle_blocks=len(blocks),tilings=tilings))


def cycle_rank_locality(G,prefer_cactus=True):
    """Compute kappa and, if kappa>2, an extremal minimal bad matching.

    Worst-case arithmetic/combinatorial time O(4^beta (n+m)), up to polynomial
    factors for the bitset implementation; beta=m-n+c.  Streaming enumeration
    does not store all cycles or all covers.
    """
    if G.is_multigraph() or nx.number_of_selfloops(G):raise ValueError('Simple graph required')
    if prefer_cactus:
        try:return cactus_locality(G)
        except ValueError:pass
    labels=list(G);H=nx.relabel_nodes(G,{v:i for i,v in enumerate(labels)},copy=True)
    stats=dict(method='cycle_space',beta=H.number_of_edges()-len(H)+nx.number_connected_components(H),
               space_vectors=0,induced_cycles=0,leaf_covers=0,max_conflicts=0,tilings=0)
    best=2;witness=[];bestcycle=None
    for component in nx.connected_components(H):
        J=H.subgraph(component);edges=[tuple(sorted(e)) for e in J.edges()]
        index={e:i for i,e in enumerate(edges)}
        basis=[]
        for C in nx.cycle_basis(J):
            v=0
            for u,w in zip(C,C[1:]+C[:1]):v^=1<<index[tuple(sorted((u,w))) ]
            basis.append(v)
        beta=len(basis);vector=0
        for z in range(1,1<<beta):
            # Consecutive binary-reflected Gray vectors differ in bit lsb(z).
            vector^=basis[(z&-z).bit_length()-1];stats['space_vectors']+=1
            ss=vector;chosen=[]
            while ss:
                bit=ss&-ss;ss^=bit;chosen.append(edges[bit.bit_length()-1])
            C=_ordered_cycle(chosen)
            if C is None:continue
            Cset=set(C)
            if any(sum(v in Cset for v in H[u])!=2 for u in C):continue
            stats['induced_cycles']+=1
            parent={}
            for v in component-Cset:
                neigh=[u for u in H[v] if u in Cset]
                if len(neigh)==1:parent[v]=neigh[0]
            # Same-parent leaves cannot both be selected anyway, so their edges
            # do not create an extra conflict between different monomers.
            conflicts=[(u,v) for u,v in J.edges() if u in parent and v in parent and parent[u]!=parent[v]]
            f=len(conflicts);stats['max_conflicts']=max(stats['max_conflicts'],f)
            assert f<=beta-1
            for mask in range(1<<f):
                removed={e[(mask>>i)&1] for i,e in enumerate(conflicts)}
                representative={}
                for v,u in parent.items():
                    if v not in removed:representative.setdefault(u,v)
                til=best_cycle_tiling(C,representative);stats['leaf_covers']+=1;stats['tilings']+=1
                if til is None:continue
                q=(len(C)+til[0])//2
                if q>best:
                    best=q;witness=til[2]+[(u,representative[u]) for u in til[1]];bestcycle=C
    return dict(kappa=best,witness=[(labels[u],labels[v]) for u,v in witness],
                cycle=None if bestcycle is None else [labels[v] for v in bestcycle],stats=stats)
