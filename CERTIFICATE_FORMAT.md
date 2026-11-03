# Indexed-permutation decoding

A certificate with an explicit `x` uses that image list directly. For an indexed record `(i,j,parity,candidate)`, take the one-based pool entry `u`. If its sign equals `parity`, set `x=u`; otherwise set `x=u(1 2)`. Permutations act on the right, so multiplication by `(1 2)` interchanges **image values** 1 and 2, not positions in the list. The resulting `x` has the requested sign. It is this decoded `x` that labels the leaf identification and whose conjugate intersection is tested.

The current stored pool entries are even. A positive-sign record therefore uses the pool entry unchanged, and a negative-sign record uses its product with `(1 2)`. For degree 40, `(i,j,parity,candidate)=(90,90,-1,1)` selects an even entry and tests its product with `(1 2)`.

This is the existing behavior in `computation/graph_engine.cpp`, in `check` mode immediately after reading each query; the following assertion checks the decoded sign. The source, executable, pools, query records, and intersection acceptance criterion are unchanged in revision 3.
