"""Exact weighted maximal acyclic MATCHING counts on once-subdivided graphs.

Input is a simple subcubic half-graph K.  Cost (u,v) is the independently
specified cost of matching u to the subdivision vertex of uv.  Counts refer
to actual matchings, not just their sets of unmatched original vertices.

Algorithm: constant-state oriented path tables + at most seven states per
cubic branch. Full core paths determine forest connectivity; a dirty single
cut is legal exactly when its two ends reconnect through the full core.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import product
from typing import Any, Mapping
import networkx as nx

# Modes at a path endpoint: ROOT, OUT_ELSEWHERE, OUT_ON_THIS_PATH.
ROOT, ELSE, ON = 0, 1, 2
# Local symbols: root, outgoing left, outgoing right, outgoing off-path.
ZERO, LEFT, RIGHT, OFF = 0, 1, 2, 3

@dataclass(frozen=True)
class Entry:
    cost: int
    optimal_count: int
    total_count: int
    arcs: tuple[tuple[Any, Any], ...]


def plus(entries):
    entries = [e for e in entries if e is not None]
    if not entries:
        return None
    best = min(e.cost for e in entries)
    chosen = next(e for e in entries if e.cost == best)
    return Entry(best, sum(e.optimal_count for e in entries if e.cost == best),
                 sum(e.total_count for e in entries), chosen.arcs)


def path_table(path, costs, left_mode, right_mode):
    """Three exact classes: full, clean cut, demanded reconnection.

    An internal vertex chooses no outgoing edge, left, or right. An edge is
    selected by exactly one endpoint or not selected. Two endpoints selecting
    the same edge is disallowed. A missing edge incident to a root is dirty.
    Two or more missing edges with any dirty one is locally impossible in a
    maximal forest, irrespective of the rest of the graph.
    """
    if len(path) < 2:
        raise ValueError('Path must contain at least one edge')
    first = ZERO if left_mode == ROOT else OFF if left_mode == ELSE else RIGHT
    last = ZERO if right_mode == ROOT else OFF if right_mode == ELSE else LEFT
    dp = {(first, 0, False): (0, 1, 1)}
    parents = [dict()]
    for i in range(1, len(path)):
        u, v = path[i-1], path[i]
        nxt, back = {}, {}
        symbols = (last,) if i == len(path)-1 else (ZERO, LEFT, RIGHT)
        for state, (val, ways, total) in dp.items():
            prev, cuts, dirty = state
            for sym in symbols:
                forward, backward = prev == RIGHT, sym == LEFT
                if forward and backward:
                    continue
                missing = not (forward or backward)
                nc = min(2, cuts + int(missing))
                nd = dirty or (missing and (prev == ZERO or sym == ZERO))
                if nc >= 2 and nd:
                    continue
                extra = costs[(u, v)] if forward else costs[(v, u)] if backward else 0
                ns, nv = (sym, nc, nd), val + extra
                old = nxt.get(ns)
                if old is None:
                    nxt[ns] = (nv, ways, total)
                    back[ns] = state
                else:
                    ov, ow, ot = old
                    if nv < ov:
                        nxt[ns] = (nv, ways, ot+total)
                        back[ns] = state
                    elif nv == ov:
                        nxt[ns] = (ov, ow+ways, ot+total)
                    else:
                        nxt[ns] = (ov, ow, ot+total)
        dp = nxt
        parents.append(back)
    buckets = {'full': [], 'clean': [], 'demand': []}
    for state, (value, ways, total) in dp.items():
        _, cuts, dirty = state
        category = 'full' if cuts == 0 else 'demand' if dirty else 'clean'
        symbols, st = [state[0]], state
        for i in range(len(path)-1, 0, -1):
            st = parents[i][st]
            symbols.append(st[0])
        symbols.reverse()
        arcs = []
        for i in range(len(path)-1):
            if symbols[i] == RIGHT:
                arcs.append((path[i], path[i+1]))
            if symbols[i+1] == LEFT:
                arcs.append((path[i+1], path[i]))
        buckets[category].append(Entry(value, ways, total, tuple(arcs)))
    return {key: plus(value) for key, value in buckets.items()}


def _walk_cycle(H):
    start = next(iter(H))
    path = [start]
    prev, cur = start, next(iter(H[start]))
    while cur != start:
        path.append(cur)
        nxt = next(v for v in H[cur] if v != prev)
        prev, cur = cur, nxt
    path.append(start)
    return path


def _unbranched(H, costs):
    if H.number_of_edges() == 0:
        return Entry(0, 1, 1, ()), dict(branch_assignments=0, paths=0, cycle=False)
    ends = [v for v in H if H.degree(v) == 1]
    candidates = []
    if ends:
        path = [ends[0]]
        prev, cur = None, ends[0]
        while True:
            opts = [x for x in H[cur] if x != prev]
            if not opts:
                break
            prev, cur = cur, opts[0]
            path.append(cur)
        for lm, rm in product((ROOT, ON), repeat=2):
            tab = path_table(path, costs, lm, rm)
            candidates.extend((tab['full'], tab['clean']))
        assignments = 4
    else:
        path = _walk_cycle(H)
        # Anchor is a single vertex appearing at both path ends.
        # It is root, points through its first incidence, or through its last.
        for lm, rm in ((ROOT, ROOT), (ON, ELSE), (ELSE, ON)):
            tab = path_table(path, costs, lm, rm)
            candidates.extend((tab['clean'], tab['demand']))
        assignments = 3
    return plus(candidates), dict(branch_assignments=assignments, paths=1, cycle=not bool(ends))


def _extract_paths(H, branches):
    seen, paths = set(), []
    for u in branches:
        for v in H[u]:
            edge = frozenset((u, v))
            if edge in seen:
                continue
            seen.add(edge)
            path = [u, v]
            prev, cur = u, v
            while H.degree(cur) == 2:
                nxt = next(x for x in H[cur] if x != prev)
                key = frozenset((cur, nxt))
                if key in seen:
                    raise AssertionError('Malformed chain traversal')
                seen.add(key)
                path.append(nxt)
                prev, cur = cur, nxt
            paths.append(path)
    if len(seen) != H.number_of_edges():
        raise AssertionError('An edge was not assigned to a unique branch chain')
    return paths


def _branched(H, costs):
    branches = [v for v in H if H.degree(v) == 3]
    index = {v: i for i, v in enumerate(branches)}
    paths = _extract_paths(H, branches)
    core = [p for p, path in enumerate(paths) if path[-1] in index]
    core_set = set(core)
    pendant = [p for p in range(len(paths)) if p not in core_set]
    tabs, unary = {}, {}
    ports = {v: [] for v in branches}
    for p, path in enumerate(paths):
        u, v = path[0], path[-1]
        ports[u].append((p, 0))
        if p in core_set:
            ports[v].append((p, 1))
            tabs[p] = {(lm, rm): path_table(path, costs, lm, rm)
                       for lm, rm in product(range(3), repeat=2)}
        else:
            unary[p] = {}
            for lm in range(3):
                allowed = []
                for rm in (ROOT, ON):
                    tab = path_table(path, costs, lm, rm)
                    allowed.extend((tab['full'], tab['clean']))
                unary[p][lm] = plus(allowed)
    choices = []
    for v in branches:
        options = [None]  # A true root.
        for p, side in ports[v]:
            if p in core_set:
                options += [(p, side, False), (p, side, True)]
            else:
                options.append((p, side, None))
        assert len(options) <= 7
        choices.append(options)
    best, best_count, all_count, best_entries = None, 0, 0, None
    stats = dict(branches=len(branches), paths=len(paths), core_paths=len(core),
                 pendant_paths=len(pendant), branch_assignments=0,
                 cyclic_full_cores=0, consistent_assignments=0)
    for actions_tuple in product(*choices):
        stats['branch_assignments'] += 1
        actions = dict(zip(branches, actions_tuple))
        full = {a[0] for a in actions_tuple if a is not None and a[2] is True}
        if any(a is not None and a[2] is False and a[0] in full for a in actions_tuple):
            continue
        adj = [[] for _ in branches]
        for p in full:
            u, v = index[paths[p][0]], index[paths[p][-1]]
            adj[u].append(v)
            adj[v].append(u)
        labels, ncomp = [-1]*len(branches), 0
        for i in range(len(branches)):
            if labels[i] >= 0:
                continue
            labels[i] = ncomp
            stack = [i]
            while stack:
                x = stack.pop()
                for y in adj[x]:
                    if labels[y] < 0:
                        labels[y] = ncomp
                        stack.append(y)
            ncomp += 1
        # Valid also for loops and parallel full paths in the compressed core.
        if len(full) != len(branches)-ncomp:
            stats['cyclic_full_cores'] += 1
            continue
        def mode(v, p, side):
            a = actions[v]
            return ROOT if a is None else ON if a[:2] == (p, side) else ELSE
        selected = []
        possible = True
        for p, path in enumerate(paths):
            u, v = path[0], path[-1]
            lm = mode(u, p, 0)
            if p not in core_set:
                entry = unary[p][lm]
            else:
                rm = mode(v, p, 1)
                table = tabs[p][(lm, rm)]
                if p in full:
                    entry = table['full']
                elif labels[index[u]] == labels[index[v]]:
                    entry = plus((table['clean'], table['demand']))
                else:
                    entry = table['clean']
            if entry is None:
                possible = False
                break
            selected.append(entry)
        if not possible:
            continue
        stats['consistent_assignments'] += 1
        val, cnt, tot = 0, 1, 1
        for entry in selected:
            val += entry.cost
            cnt *= entry.optimal_count
            tot *= entry.total_count
        all_count += tot
        if best is None or val < best:
            best, best_count, best_entries = val, cnt, selected
        elif val == best:
            best_count += cnt
    if best is None:
        raise AssertionError('Every finite graph has a maximal acyclic matching')
    arcs = tuple(arc for entry in best_entries for arc in entry.arcs)
    return Entry(best, best_count, all_count, arcs), stats


def solve_incidence(K: nx.Graph, costs: Mapping | None = None) -> dict:
    """Minimize arbitrary signed incidence costs and count actual matchings.

    Returns one matching as directed original edges u->v, meaning the edge
    u--s_{uv} in S(K). Empty graphs have the unique empty solution. Input
    components are optimized independently. Arithmetic-operation complexity
    O(n+m+sum_i 7**b_i*(b_i+1)); big-integer bit costs must be charged separately.
    """
    if K.is_directed() or K.is_multigraph() or nx.number_of_selfloops(K):
        raise ValueError('K must be simple and undirected')
    if any(K.degree(v) > 3 for v in K):
        raise ValueError('Maximum degree of K must be at most three')
    if costs is None:
        costs = {(u, v): 1 for u in K for v in K[u]}
    for u in K:
        for v in K[u]:
            if (u, v) not in costs or not isinstance(costs[(u, v)], int):
                raise ValueError('Provide an exact integer cost for each directed incidence')
    value, optimal_count, total_count, arcs, stats = 0, 1, 1, [], []
    for comp in nx.connected_components(K):
        H = K.subgraph(comp)
        if any(H.degree(v) == 3 for v in H):
            out, st = _branched(H, costs)
        else:
            out, st = _unbranched(H, costs)
        value += out.cost
        optimal_count *= out.optimal_count
        total_count *= out.total_count
        arcs.extend(out.arcs)
        stats.append(st)
    return dict(value=value, optimal_count=optimal_count, total_count=total_count,
                matching_arcs=arcs, component_stats=stats,
                branch_assignments=sum(s['branch_assignments'] for s in stats))


def subdivide(K: nx.Graph):
    """Construct a disjointly tagged S(K) and an incidence-to-edge map."""
    G = nx.Graph()
    G.add_nodes_from(('v', u) for u in K)
    lookup = {}
    for i, (u, v) in enumerate(K.edges()):
        s = ('e', i)
        for a, b in ((u, v), (v, u)):
            edge = (('v', a), s)
            G.add_edge(*edge)
            lookup[(a, b)] = edge
    return G, lookup
