clear
clc
close all

%% Primal-dual optimization problem
%Primal problem
%min c'x
%s.t. Gx + s = h, Ax = b, s >= 0

%Dual problem
%max -h'z - b'y
%s.t. G'z + A'y + c = 0, z >= 0

m = 10; %Numer of rows
n = 10; %Number of columns

%% ECOS Test
dims.l = 0; dims.q = 10;
ECOS_info.exitflag = -999;
while ECOS_info.exitflag ~= 0
    %Test data
    c = rand(m,1);
    b = rand(m,1);
    h = rand(m,1);
    A = rand(m,n);
    G = rand(m,n);

    A = sparse(A);
    G = sparse(G);
    [x_ECOS,y_ECOS,ECOS_info,s_ECOS,z_ECOS] = ecos(c,G,h,dims,A,b);
    optimal_objective_value_ECOS = c'*x_ECOS
    ECOS_info
end
%% End of ECOS Test

c
b
h
A_full = full(A)
G_full = full(G)