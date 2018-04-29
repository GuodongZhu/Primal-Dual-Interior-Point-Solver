# Primal-Dual Interior-Point Solver for SOCPs

The following code is based on the writings in a PhD thesis by Alexander Domahidi titled "Methods and Tools for Embedded Optimization and Control", ETH Zurich, 2013. It solves general second-order cone programs of the form

```
min  c'x
s.t. Ax = b
     Gx <= h
```

by implementing a modified primal-dual Mehrotra predictor-corrector interior-point algorithm for second-order conic optimization problems, using Nesterov-Todd scaling and self-dual embedding. Where `<=` is a generalized inequality such that `h - Gx` belongs to a cone `K`. A test case is given in `run_test.m`, which provides test data and makes a call to the main optimization algorithm in `PDIP_solver.m`.

## File descriptions:
* `run_test.m`: provides test data and calls PDIP_solver.m in order to find optimal solution to the problem
* `PDIP_solver.m`: main optimization algorithm used to solve generalized second-order cone programs
* `soc_dot.m`: cone product operator for second-order cones (similar to dot product for vectors)
* `soc_dot_inv.m`: inverse to the cone product operator for second-order cones
* `w_soc.m`: computes Nesterov-Todd scaling matrix, which allows for long steps along search direction
