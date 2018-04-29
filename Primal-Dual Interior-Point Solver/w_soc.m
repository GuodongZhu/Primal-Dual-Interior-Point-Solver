%% Computes SOC Nesterov-Todd scaling matrix
function W = w_soc(s,z,m)
    zbar = z/sqrt((z(1)^2)-(z(2:end)'*z(2:end)));
    sbar = s/sqrt((s(1)^2)-(s(2:end)'*s(2:end)));
    gamma = sqrt((1+(zbar'*sbar))/2);
    wbar = (1/(2*gamma))*(sbar+[zbar(1); -zbar(2:end)]);
    eta = (((s(1)^2)-((s(2:end)')*s(2:end)))/((z(1)^2)-((z(2:end)')*z(2:end))))^(1/4);

    W = eta*[wbar(1), wbar(2:end)'; wbar(2:end), (eye(m-1)+(1+wbar(1))\wbar(2:end)*(wbar(2:end)'))];
end