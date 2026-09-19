# Acyclic Matching Locality

Exact algorithms, independent small-instance checks, and recorded experimental
data accompanying the research manuscript **Block Locality and Exact Optimization
and Counting of Maximal Acyclic Matchings in Subdivision Graphs**.

**Hao Tang**, School of Software, Jilin University  
Contact: tanghao5524@mails.jlu.edu.cn

This is a research software release, not a statement of journal acceptance.
The manuscript and third-party literature are not bundled here.

## Install and verify

Use Python 3.11 or later. From the repository root:

```sh
python -m pip install -r requirements-validation.txt
python -u validation/run_checks.py
```

The pinned environment uses NetworkX 3.1. The historical parameter measurements
used Python 3.13.5 / NetworkX 3.6.1 on Linux; the recorded revision checks used
Python 3.11.7 / NetworkX 3.1 on Windows. See the environment records in `data/`.
The main validation script supports Windows and Linux and takes roughly ten
seconds on the development machine; elapsed times depend on the machine.
Do not run Python with `-O`: checks use assertions.

Successful validation ends with:

```text
ALL PASSED 62 answer cases; 76 structural graphs; 4 timing baselines
```

It rechecks 14 historical small instances, tests 48 graph/cost combinations,
enumerates all 76 labelled simple base graphs with at most four vertices, and
compares four instances with a separate root/parent enumeration baseline.
The matching oracle enumerates actual ordinary matchings, checks induced
acyclicity by union/find, and checks all legal single-edge extensions. These
finite checks support implementation correctness; they do not replace proofs.
This oracle is separate code in the same project, not an external audit.

New runs write to ignored `runs/revision/`. Committed data are never overwritten.

## Algorithms

- `code/block_locality.py`: pendant-augmented biconnected block decomposition,
  exact locality order, and an obstruction witness when the order exceeds two.
- `code/cycle_rank_locality_v7.py`: the cycle-space engine used inside each block.
- `code/incidence_matching.py`: minimum cost, number of optimal solutions,
  total number of solutions, and an optimal witness for inclusion-maximal
  acyclic matchings on once-subdivisions of simple subcubic base graphs.
- `code/direct_checks.py`: historical direct matching enumerator used by the
  original small-instance timing runner.

An acyclic matching is one whose saturated vertices induce a forest.
Maximal always means inclusion-maximal. Counting distinguishes actual
matchings, even when their unmatched original vertices coincide.

### Weighted subdivision example

Run this Python code from the repository root:

```python
import sys
sys.path.insert(0, "code")
import networkx as nx
from incidence_matching import solve_incidence, subdivide

K = nx.cycle_graph(3)
costs = {(u, v): 1 for u in K for v in K[u]}
result = solve_incidence(K, costs)
print(result["value"], result["optimal_count"], result["total_count"])
# 2 9 9
S, incidence_edges = subdivide(K)
matching = [incidence_edges[a] for a in result["matching_arcs"]]
print(matching)  # an optimal matching in S(K); tie-breaking may vary
```

The solver input is the **base graph K**, not its subdivision. A key `(u, v)`
means matching original vertex `u` to the subdivision vertex of edge `uv`.
Opposite incidences `(u, v)` and `(v, u)` may have independent signed integer
costs. Supply every incidence explicitly. Integer counts use Python's arbitrary
precision arithmetic. The weighted solver requires maximum degree at most three.

### Block locality example

```python
import sys
sys.path.insert(0, "code")
import networkx as nx
from block_locality import block_locality

G = nx.cycle_graph(3)
G.add_edges_from((v, v + 3) for v in range(3))
result = block_locality(G)
assert result["kappa"] == 3
print(result["witness"])
```

This routine accepts finite simple undirected graphs. The locality order is
normalized to be at least two; an output of two does not assert that a
two-edge obstruction exists. The algorithms are exponential in their relevant
structural parameters; they are not polynomial-time solvers for arbitrary input.

## Reproduce historical experiments

The original timing runners require Linux or WSL (`SIGALRM` and `/proc/cpuinfo`):

```sh
python -u experiments/controlled_parameters.py
python -u experiments/small_oracle_benchmarks.py
```

For the historical dependency version, install `networkx==3.6.1` in a separate
Python 3.13 environment. These commands regenerate 26 controlled parameter
configurations and 14 small-instance comparisons in `runs/parameters/`.
Each timing uses the median of three sequential calls; Python allocation peaks
come from a separate `tracemalloc` call and are not process RSS. The per-call
timeout is 90 seconds. Graph construction and validation are excluded from
solver timings. Compare exact answers, not machine-dependent elapsed times.

## Recorded data and provenance

- `data/original/`: unchanged inputs, measurements, and environment for the
  16 September 2026 experiment series.
- `data/revision/`: unchanged inputs and results for the 19 September 2026
  revision checks and four baseline timing comparisons.
- `docs/provenance.md`: source mapping, adaptations, and validation scope.
- `manifest.json`: SHA-256 hashes of the release files and source provenance.

The four algorithm modules are byte-for-byte copies of the frozen research
implementation. Public runners adapt directory paths, ensure fresh output
directories, and clarify platform requirements. Historical result hashes
refer to the original scripts, not the adapted release runner.

MIT license. See `CITATION.cff` for software citation metadata.
