# Kourovka Notebook Problem 21.3(b): computational supplement

This local review package accompanies **Soluble subgroup pairs in symmetric and alternating groups**, revision 3, 3 November 2026. The theorem states: for every integer n >= 21 and soluble H,K <= S_n, each sign class contains x with H ∩ x⁻¹Kx = 1. In particular the second question of Anagnostopoulou-Merkouri and Burness in Kourovka Notebook Problem 21.3 has a positive answer.

The manuscript is distributed separately. Its exact source and PDF hashes are recorded in `PAPER_BINDING.json`; neither manuscript file is included here. This public repository distributes the computational supplement separately from the manuscript.

The finite proof has 107 degree claims, 5,086 strictly decreasing induction dependencies, 38,972,428 unordered pairs, and 38,972,586 concrete sign witnesses. A fresh full replay checked all 301,061 GAP witnesses in degrees 21–26 and all 38,671,525 graph queries in degrees 27–127. All passed. The infinite range n >= 128 follows from the manuscript's uniform estimates and the finite rational computations here, including 7,251 primitive coefficient comparisons, 465 density boxes, and 267 odd-prime rows.

Start with `REPRODUCE.md`. `CERTIFICATE_INDEX.json` identifies every finite input. `PROOF_CERTIFICATE_MAP.md` links the mathematical assertions to the code and data. `MANIFEST.json` binds file sizes and SHA-256 values. `SOURCE_PATH_MAP.json` preserves the original-to-public path correspondence. Completion records describe executed checks; the verification programs repeat the corresponding acceptance computations.

Python 3.11 or newer is recommended; the Python checks use the standard library. Graph replay additionally requires Linux or WSL and the supplied nauty/Traces 2.9.3 executable, or a local build from the included official source. GAP replay requires GAP 4.12.1. See `LICENSING.md` for the supplied third-party software.

Revision 2 expands proof details and notation; the certificate inputs and verification code are unchanged. Its exact degree-128 calculations were repeated successfully.

Revision 3 specifies the existing indexed-pool parity decoding; see `CERTIFICATE_FORMAT.md`. Mathematical data and verifier behavior are unchanged.
