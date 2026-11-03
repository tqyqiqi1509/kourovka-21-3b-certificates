# Proof-to-certificate map

| Mathematical step | Inputs | Verification | Expected result |
|---|---|---|---|
| Affine forest envelopes and primitive counts | `structural_reduction/`, field-coordinate JSON files | `verify.py models` | Complete finite partitions; all 703,569 generators recognized |
| Lifting and finite residual classification | `computation/lifting/` and `computation/hyperplane_lifting/` induction and inheritance JSON | `verify.py models`; dependency records in the indices | Degrees 21–127 close by decreasing-degree induction; 5,086 dependencies |
| All pair and sign slots | `CERTIFICATE_INDEX.json`, actual group files and certificates | `verify.py coverage` | 107 degrees; 38,972,428 pairs; 38,972,586 witnesses |
| Exact graph encoding of actual generators | Graph files, word families, coordinates and actual group generators | `verify_graph_models.py` | All supplied graphs reconstructed; each actual generator extends |
| Finite disjoint intersections, n=21–26 | Six GAP group/certificate families | `verify.py gap` | 301,061 exact intersections equal the identity |
| Finite disjoint intersections, n=27–127 | 101 graph families, pools, queries/KC4 certificates | `verify.py graph` | 38,671,525 accepted exact graph checks |
| Atom-density bounds | `prime_order_referee/prime_atom_densities.py` and density table | `verify.py math`, independent density audit in `verify.py models` | Regenerated table agrees exactly |
| Primitive strengthened involution inequality | Exact affine counts for 69 degrees | `verify.py math` | 7,251 coefficient comparisons |
| Conditional connected-involution and odd-prime bounds | `connected_128_certificate.py`, density table | `verify.py math` | 465 boxes, 267 rows, maximum < 999/1000; empty-core < 1/2 |
| Independent outward intervals | `analysis/connected_128_interval_verifier.py` | `verify.py math` | Every interval encloses the exact quantities and meets the strict bounds |
| Connected coefficient interpretation | Small direct port matchings | `verify.py math` | Connected counts agree with the formal logarithm |
| All n >= 128 | Uniform local-lemma bound, coefficientwise inequalities, monotonicity, analytic tails in the manuscript | Mathematical arguments plus the preceding bounded rational calculation | No per-degree enumeration above 127 |

`verification_records/` contains full fresh intersection replay records and independent audit outputs. The manuscript gives the definitions, exhaustive ranges, recurrence formulas, graph acceptance implication, and proofs that these finite checks suffice.
