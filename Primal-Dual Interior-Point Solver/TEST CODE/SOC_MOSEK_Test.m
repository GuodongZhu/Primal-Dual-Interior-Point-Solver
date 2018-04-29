% clear
% clc
% close all

%options = sdpsettings('solver','mosek','verbose',1);
options = sdpsettings('solver','ecos','verbose',1);

m = 4979;
x = sdpvar(m,1);

% objective = c'*x;
% constraints = [A*x == b];

objective = c'*x;
constraints = [
    A*x == b;
    x(1:2502) >= 0;
    cone(x(2504:2504+2472), x(2503));
    cone(x(4978:4979), x(4977))
    ];

optimize(constraints,objective,options);
x_opt = value(x);