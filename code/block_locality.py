"""Locality order via pendant-augmented biconnected blocks.

The numerical engine is the explicitly preserved v7 cycle-space algorithm.
The new part is the lossless block reduction, including external leaf flags.
Input: finite simple undirected NetworkX graph with arbitrary hashable labels.
"""
from __future__ import annotations
from typing import Iterator
import networkx as nx
from cycle_rank_locality_v7 import cycle_rank_locality


def augmented_cycle_blocks(G: nx.Graph) -> Iterator[tuple[nx.Graph, dict]]:
    if G.is_directed() or G.is_multigraph() or nx.number_of_selfloops(G):
        raise ValueError("A finite simple undirected graph is required")
    for raw in nx.biconnected_component_edges(G):
        edges = list(raw)
        if len(edges) == 1:
            continue  # No cycle can be supported by a bridge.
        B = nx.Graph()
        B.add_edges_from(edges)
        original = set(B)
        boundary = []
        representatives = {}
        used_external = set()
        for v in list(B):
            inside_degree = B.degree(v)
            if G.degree(v) > inside_degree:
                x = next(x for x in G[v] if x not in original)
                if x in used_external:
                    raise AssertionError("An external component attaches twice to a block")
                representatives[v] = x
                used_external.add(x)
                boundary.append(v)
        beta = B.number_of_edges() - B.number_of_nodes() + 1
        for v, x in representatives.items():
            B.add_edge(v, x)
        yield B, dict(block_vertices=len(original), block_edges=len(edges),
                      beta=beta, boundary=boundary,
                      representatives=representatives)


def block_locality(G: nx.Graph) -> dict:
    """Return exact normalized kappa and an extremal bad matching if >2.

    Combinatorial bound O(n+m + sum_B 4^beta_B (|V(B)|+|E(B)|)).
    The inherited Python engine uses arbitrary-precision bitsets; their bit
    cost is not silently assumed constant.
    """
    stats = dict(n=G.number_of_nodes(), m=G.number_of_edges(),
                 beta_total=G.number_of_edges()-G.number_of_nodes()+nx.number_connected_components(G),
                 beta_max=0, cyclic_blocks=0, boundary_flags=0,
                 augmented_vertices=0, space_vectors=0, induced_cycles=0,
                 leaf_covers=0)
    answer = dict(kappa=2, witness=[], cycle=None, winning_block=None)
    for index, (H, info) in enumerate(augmented_cycle_blocks(G)):
        stats['cyclic_blocks'] += 1
        stats['beta_max'] = max(stats['beta_max'], info['beta'])
        stats['boundary_flags'] += len(info['boundary'])
        stats['augmented_vertices'] += len(H)
        # Keep the same cycle-space engine, rather than silently switching to
        # the old cactus-specific formula in the new wrapper.
        result = cycle_rank_locality(H, prefer_cactus=False)
        for key in ('space_vectors', 'induced_cycles', 'leaf_covers'):
            stats[key] += result['stats'].get(key, 0)
        if result['kappa'] > answer['kappa']:
            answer.update(kappa=result['kappa'], witness=result['witness'],
                          cycle=result['cycle'], winning_block=index)
    answer['stats'] = stats
    return answer
