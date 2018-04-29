%% Computes SOC dot product operator
%x and y are elements in a second-order cone Q^n
function xy_dot = soc_dot(x,y)
    xy_dot = [(x'*y); ((x(1)*y(2:end)) + (y(1)*x(2:end)))];
end