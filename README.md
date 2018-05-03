# Primal-Dual Interior-Point Solver for SOCPs

This is a project that I completed between January and March 2018.

The following code implements techniques described in the PhD thesis by Alexander Domahidi titled "Methods and Tools for Embedded Optimization and Control", ETH Zurich, 2013. It solves general second-order cone programs of the form

```
min  c'x
s.t. Ax = b
     Gx <= h
```

Where `<=` is a generalized inequality such that `h - Gx` is contained in the cone `K`. It implements a modified primal-dual Mehrotra predictor-corrector interior-point algorithm for second-order conic optimization problems, using Nesterov-Todd scaling and self-dual embedding. A test case is given in `run_test.m`, which provides test data and makes a call to the main optimization algorithm in `PDIP_solver.m`.

## File descriptions:
* `run_test.m`: provides test data and calls `PDIP_solver.m` in order to find an optimal solution to the test problem
* `PDIP_solver.m`: main primal-dual interior-point optimization algorithm that is used to solve generalized second-order cone programs
* `soc_dot.m`: cone product operator for second-order cones (similar to dot product for vectors)
* `soc_dot_inv.m`: inverse to the cone product operator for second-order cones
* `w_soc.m`: computes Nesterov-Todd scaling matrix, which allows for long steps along search direction
