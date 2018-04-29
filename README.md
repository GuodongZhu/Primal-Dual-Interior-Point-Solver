# Primal-Dual Interior-Point Solver for Second-Order Conic Optimization Problems

The following code implements a primal-dual interior-point solver to find the optimal solution to general second-order cone programs of the form

min c'x
s.t. Ax = b
Gx <= h

where <= is a generalized inequality such that h - Gx belongs to a cone K. A test case is given in run_test.m, which provides test data and makes a call to conic_primal_dual_interior_point_solver.m

# File descriptions:
* run_test.m: provides test data and calls PDIP_solver.m in order to find optimal solution to the problem
* PDIP_solver.m: main optimization algorithm used to solve generalized second-order cone programs
* soc_dot.m: cone product operator for second-order cones (similar to dot product for vectors)
* soc_dot_inv.m: inverse to the cone product operator for second-order cones
* w_soc.m: computes Nesterov-Todd scaling matrix, which allows for long steps along search direction
