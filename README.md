# CAMBLateDE - CAMB for Late Dark Energy
Authors: Diogo Souza, Vivian Miranda, João Rebouças

If you use this code, please cite arXiv 2408.14628

This is a CAMB modification that implements parametrizations for the equation of state for late dark energy.
The available parametrization are: 
1) Constant w
2) w0wa
3) Constant bin w     - 2, 3, 5 and 10 bins (item 8 generalizes this parametrization)
4) Linear bin w       - 2, 3, 5 and 10 bins
5) Quadratic bin w    - 2, 3 and 5 bins
6) Cubic bin w        - 2, 3 and 5 bins
7) Hyperbolic tangent - Equation 21 in Planck's paper arXiv 1502.01590 
8) Constant bin w     - Arbitrary numbers of bins
9) Chebyshev expansion of w - see arXiv 2405.04216v1 (Obs: no smooth transition of w to C.C. regime)

## Using the LateDE CAMB
Examples of the parametrizations are provided at `docs/LateDE_demo.ipynb`  
Extensive tests of this modified CAMB is given at `docs/LateDE_test.ipynb`