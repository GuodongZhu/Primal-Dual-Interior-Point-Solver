%% Computes SOC inverse operator
%v = u\w for u*v = w
function v = soc_dot_inv(u,w)
    rho = (u(1)^2)-(u(2:end)'*u(2:end));
    nu = u(2:end)'*w(2:end);
    v = (1/rho)*[(u(1)*w(1))-nu; ((nu/u(1)-w(1))*u(2:end) + (rho/u(1))*w(2:end))];
end