%% Primal-Dual Interior-Point solver run using test data example
clear
clc
close all

%% Test data
c = [0.6437;
    0.3012;
    0.1323;
    0.5516;
    0.0240;
    0.8786;
    0.6809;
    0.8723;
    0.1446;
    0.9457];

b = [0.1479;
    0.5280;
    0.8557;
    0.4120;
    0.0798;
    0.5292;
    0.4154;
    0.8352;
    0.6333;
    0.9801];

h = [0.9488;
    0.2785;
    0.9699;
    0.9895;
    0.4555;
    0.7889;
    0.2354;
    0.0887;
    0.0460;
    0.0541];

A = [0.8222    0.7254    0.7230    0.4955    0.1397    0.7388    0.4913    0.3433    0.0438    0.9944;
    0.7423    0.7267    0.8084    0.8450    0.6076    0.6584    0.0214    0.4952    0.7393    0.3896;
    0.0339    0.6516    0.6718    0.4338    0.7693    0.2714    0.6795    0.0388    0.9148    0.3583;
    0.1103    0.7586    0.1770    0.2347    0.3739    0.8107    0.9605    0.3538    0.6187    0.6582;
    0.9254    0.5478    0.6986    0.4196    0.7123    0.9682    0.5862    0.8175    0.4364    0.6855;
    0.3001    0.6269    0.7691    0.0115    0.5303    0.5064    0.0578    0.7933    0.2295    0.4369;
    0.1950    0.6240    0.2869    0.7544    0.2344    0.7539    0.1028    0.0238    0.5510    0.2156;
    0.4000    0.9475    0.7165    0.9831    0.6654    0.6709    0.0267    0.8654    0.3994    0.1801;
    0.7107    0.4688    0.2281    0.7398    0.5863    0.3113    0.7721    0.2210    0.9477    0.0791;
    0.3909    0.3208    0.9991    0.1543    0.7281    0.4812    0.3831    0.4288    0.9458    0.1680];

G = [0.9745    0.2011    0.4182    0.4080    0.6390    0.1863    0.7726    0.0619    0.1233    0.8830;
    0.7253    0.6289    0.1611    0.0118    0.4713    0.6848    0.5106    0.6895    0.0586    0.3879;
    0.6812    0.7978    0.8686    0.5170    0.6874    0.9078    0.1599    0.5284    0.2017    0.2258;
    0.0621    0.8675    0.1529    0.2025    0.4369    0.5784    0.9477    0.1642    0.9744    0.6361;
    0.5628    0.5362    0.1015    0.5413    0.7974    0.8651    0.4431    0.2119    0.9307    0.3565;
    0.4718    0.5652    0.9824    0.2673    0.4797    0.6634    0.4264    0.1737    0.4780    0.8410;
    0.1460    0.7218    0.9563    0.0349    0.1865    0.7092    0.3792    0.2391    0.8580    0.9389;
    0.2667    0.1366    0.4271    0.3942    0.5661    0.1966    0.7276    0.8243    0.8620    0.1726;
    0.2370    0.6274    0.1767    0.8007    0.3157    0.0046    0.3748    0.2518    0.4501    0.3368;
    0.0977    0.9900    0.3493    0.7551    0.9513    0.0688    0.1381    0.2713    0.3232    0.0511];

m = 10; %Numer of rows
n = 10; %Number of columns
SOC_dimension = 10; %Dimension of second-order cone

%Calling PDIP solver
[x_opt,y_opt,s_opt,z_opt,primal_obj,dual_obj,duality_gap] = PDIP_solver(c,b,h,A,G,m,n,SOC_dimension)

% %% ECOS Test
% dims.l = 0; dims.q = SOC_dimension;
% A = sparse(A);
% G = sparse(G);
% [x_ECOS,y_ECOS,ECOS_info,s_ECOS,z_ECOS] = ecos(c,G,h,dims,A,b);
% optimal_objective_value_ECOS = c'*x_ECOS
% %% End of ECOS Test



