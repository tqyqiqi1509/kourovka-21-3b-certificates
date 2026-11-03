# Reproduction

Run from the repository root. Every assertion must remain enabled: do not use Python's `-O` option. Use fresh output directories for new graph replays. Generated reports are separate from the frozen input files.

## 1. Payload integrity

```sh
python verify.py integrity --output verification_results/integrity
```

This checks every manifest member. It establishes file identity, not the mathematical conclusion.

## 2. Exhaustive coverage and group models

```sh
python verify.py coverage --output verification_results/coverage
python verify.py models --output verification_results/models
python verify_graph_models.py verification_results/graph_models
```

Coverage verifies all 107 degrees, actual permutations and parities, and all 38,972,586 required pair/sign slots. The compact KC4 decoder checks reconstructed stream hashes and coverage. Models independently reconstructs the envelope partitions, primitive counts and affine-wreath membership of all 703,569 actual generators. The graph-model check regenerates every supplied graph, compares all vertices/colors/edges and extends every actual graph-envelope generator to an automorphism.

## 3. The degree-128 rational calculations

```sh
python verify.py math --output verification_results/math
```

This regenerates the density table, checks all 7,251 primitive involution coefficients, computes all 465 connected density boxes and 267 odd-prime rows with exact fractions, repeats the outward-interval calculation, and independently enumerates small port-matching models. The final maximum must be strictly below 999/1000; the empty-core bound must be strictly below 1/2. Decimal strings never decide acceptance. The analytic infinite tails and monotonicity in n are proved in the manuscript.

## Indexed permutation semantics

A certificate with an explicit `x` uses that image list directly. For an indexed record `(i,j,parity,candidate)`, take the one-based pool entry `u`. If its sign equals `parity`, set `x=u`; otherwise set `x=u(1 2)`. Permutations act on the right, so multiplication by `(1 2)` interchanges **image values** 1 and 2, not positions in the list. The resulting `x` has the requested sign. It is this decoded `x` that labels the leaf identification and whose conjugate intersection is tested.

The current stored pool entries are even. A positive-sign record therefore uses the pool entry unchanged, and a negative-sign record uses its product with `(1 2)`. For degree 40, `(i,j,parity,candidate)=(90,90,-1,1)` selects an even entry and tests its product with `(1 2)`.

This is the existing behavior in `computation/graph_engine.cpp`, in `check` mode immediately after reading each query; the following assertion checks the decoded sign. The source, executable, pools, query records, and intersection acceptance criterion are unchanged in revision 3.

## 4. Full graph intersections, degrees 27–127

On Linux or WSL, mark the supplied executable executable if the ZIP extraction did not preserve permissions:

```sh
chmod +x software_referee/runtime_293/graph_engine
python3 verify.py graph --jobs 4 --output verification_results/graph_full
```

On Windows with WSL configured, `python verify.py graph ...` automatically invokes the Linux executable through WSL. To test a subset first:

```sh
python verify.py graph --degrees 77 127 --output verification_results/graph_test
```

The full run decodes the KC4 query streams for degrees 27–39, then checks all 38,671,525 graph queries with exact Traces automorphism orbits. Each query must have singleton orbits on every original point. A nonzero process status, missing completion line, wrong version, or nonzero failure count rejects the replay. The supplied full records are in `verification_records/graph/`.

For a source build, extract `software_referee/nauty2_9_3.tar.gz` in `software_referee`, run `./configure` and `make nauty.a` in the extracted directory, then run `sh software_referee/runtime_293/build_engine.sh` from the repository root. The official archive provenance and hash are in `software_referee/upstream_manifest.json`.

## 5. Exact GAP intersections, degrees 21–26

From Linux or WSL with GAP 4.12.1, provide the real GAP executable and root containing `lib/init.g`:

```sh
python3 verify.py gap --gap-binary /absolute/path/to/gap \
  --gap-root /absolute/path/to/gap-root \
  --output verification_results/gap_full
```

The GAP replay checks all 301,061 concrete witness intersections, generated group orders, permutations, signs, and complete coverage. The included replay source contains a Unicode-path serialization correction: JSON uses literal Unicode in its GAP file paths. This fixes file access and changes no mathematical acceptance condition. The new completed records are in `verification_records/gap/`.

## Provenance and local changes

The source archive is `Kourovka_21_3b_research_bundle.zip`, SHA-256 `ec682eeb088c478c879ed15c459a0cdc18151100a260e1b86bd010ced3bb63c3`. This public subset retains the latest 128-threshold route and the exact finite inputs used by this manuscript. It omits abandoned routes, duplicate historical data, manuscript files, and private machine/network diagnostics.

The three independent model-audit scripts accept an explicit report destination; the residual script reads the retained `hyperplane_lifting` copy of the same degree-64–127 group files. New public verification entry points compose the existing checks. All input bytes and retained local edits are bound in `MANIFEST.json`.
