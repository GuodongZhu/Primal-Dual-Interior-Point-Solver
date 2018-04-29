%% Primal-Dual Interior-Point Solver for Second-Order Cone Programs
%--------------------------------------------------------------------------
% The following algorithm is based on the PhD thesis by Alexander Domahidi
% titled "Methods and Tools for Embedded Optimization and Control", ETH
% Zurich, 2013.
%
% Solves general second-order cone programs by implementing a modified
% primal-dual Mehrotra predictor-corrector interior-point method for
% second-order cones, using Nesterov-Todd scaling and self-dual embedding
%--------------------------------------------------------------------------

%% Primal-dual optimization problem
%The following general primal-dual second-order cone program is solved:

%Primal problem
%min c'x
%s.t. Gx + s = h
%Ax = b
%s >= 0
%where s >= 0 is a second-order cone constraint

%Dual problem
%max -h'z - b'y
%s.t. G'z + A'y + c = 0
%z >= 0
%where z >= 0 is a second-order cone constraint

function [x_opt,y_opt,s_opt,z_opt,primal_obj,dual_obj,duality_gap] = PDIP_solver(c,b,h,A,G,m,n,SOC_dimension)

    e = [1; zeros(SOC_dimension-1,1)];
    SOC_degree = e'*e;

    %% Initializing x
    AG = [zeros(m,n), A', G';
        A, zeros(m,n), zeros(m,n);
        G, zeros(m,n), -eye(m,n)];
    bh = [zeros(m,1); b; h];
    init_xs = AG\bh;
    x = init_xs(1:m);
    s_tilde = -init_xs((2*m+1):(3*m));

    %% Initializing s
    alpha_p = 0.1;
    while s_tilde(1) <= norm(s_tilde(2:m))
        s_tilde(1) = s_tilde(1) + alpha_p;
    end
    s = s_tilde;

    %% Initializing y
    AG = [zeros(m,n), A', G';
        A, zeros(m,n), zeros(m,n);
        G, zeros(m,n), -eye(m,n)];

    cm = [-c; zeros(m,1); zeros(m,1)];

    init_yz = AG\cm;
    y = init_yz((m+1):(2*m));
    z_tilde = init_yz((2*m+1):(3*m));

    %% Initializing z
    alpha_d = 0.1;
    while z_tilde(1) <= norm(z_tilde(2:m))
        z_tilde(1) = z_tilde(1) + alpha_d;
    end
    z = z_tilde;
    
    %% Initializing additional parameters
    kappa = 1;
    tau = 1;
    W = w_soc(s,z,SOC_dimension);
    lambda = W*z;

    max_iter = 100;
    update_history = [];
    mu_history = [];
    gap_history = [];
    equality_constraint_history = [];
    inequality_constraint_history = [];
    dual_constraint_history = [];
    iter_axis = [];

    %Stopping criteria
    FEASTOL = 1e-8;
    ABSTOL = 1e-8;
    RELTOL = 1e-8;
    
    %Optimal, primal infeasible and dual infeasible parameters
    OPT = 0;
    PIF = 0;
    DIF = 0;

    %% Optimization loop
    tic
    for iter = 1:max_iter
        iter
        %% 1) Evaluate residuals, gap and stopping criteria
        %Residuals
        r = [zeros(m,1); zeros(m,1); s; kappa] - [zeros(m,n), A', G', c; -A, zeros(m,n), zeros(m,n), b; -G, zeros(m,n), zeros(m,n), h; -c', -b', -h', 0]*[x; y; z; tau];
        r_x = r(1:length(x));
        r_y = r((length(x)+1):(length(x)+length(y)));
        r_z = r((length(x)+length(y)+1):(length(x)+length(y)+length(z)));
        r_tau = r(end);

        %Gap
        mu = (s'*z+kappa*tau)/(SOC_degree+1);

        %Stopping criteria
        r0x = max(1,norm(c));
        r0y = max(1,norm(b));
        r0z = max(1,norm(h));
        rho = max(-c'*x,(-b'*y-h'*z));

        %Optimality check
        if max([norm(A'*y + G'*z + c)/r0x, norm(A*x - b)/r0y, norm(G*x + s - h)/r0z]) < FEASTOL && (s'*z < ABSTOL || (s'*z)/rho < RELTOL)
            iter
            OPT = 1
            disp('OPTIMAL SOLUTION FOUND')
            break  
        end

        %Infeasibility check
        if (h'*z + b'*y) < 0 && norm(A'*y + G'*z)/r0x < FEASTOL
            PIF = 1
            disp('PRIMAL INFEASIBLE');
            break
        elseif c'*x < 0 && max([norm(A*x)/r0y,norm(G*x + s)/r0z]) < FEASTOL        
            DIF = 1
            disp('DUAL INFEASIBLE');
            break
        end

        %NT-scaling matrix
        W = w_soc(s,z,SOC_dimension);

        %% 2a) Static Newton step
        KKT = [zeros(m,n), A', G'; -A, zeros(m,n), zeros(m,n); -G, zeros(m,n), W'*W];
        b_vec = -[c; b; h];
        x_static = KKT\b_vec;
        x1 = x_static(1:length(x));
        y1 = x_static((length(x)+1):(length(x)+length(y)));
        z1 = x_static((length(x)+length(y)+1):(length(x)+length(y)+length(z)));
        %x_static - [x1; y1; z1] %Test

        %% 2b) Affine-scaling (predictor) Newton step
        d_x = r_x;
        d_y = r_y;
        d_z = r_z;
        d_tau = r_tau;
        d_s = soc_dot(lambda,lambda);
        d_kappa = kappa*tau;

        b_affine = [d_x; d_y; d_z-W'*(soc_dot_inv(lambda,d_s))];
        x_affine = KKT\b_affine;
        x2_affine = x_affine(1:length(x));
        y2_affine = x_affine((length(x)+1):(length(x)+length(y)));
        z2_affine = x_affine((length(x)+length(y)+1):(length(x)+length(y)+length(z)));
        %x_affine - [x2_affine; y2_affine; z2_affine] %Test

        numerator = (d_tau - (d_kappa/tau) + [c', b', h']*x_affine);
        denominator = ((kappa/tau) - [c', b', h']*x_static);
        delta_tau_affine = numerator/denominator;
        delta_x_affine = x2_affine + delta_tau_affine*x1;
        delta_y_affine = y2_affine + delta_tau_affine*y1;
        delta_z_affine = z2_affine + delta_tau_affine*z1;
        delta_s_affine = -W'*(soc_dot_inv(lambda,d_s) + W*delta_z_affine);
        delta_kappa_affine = -(d_kappa + kappa*delta_tau_affine)/tau;

        %% 3) Affine-scaling step-size and centering parameter
        alpha = 0.9;

        %Backtracking linesearch
        while ((s(1) + alpha*delta_s_affine(1)) - (norm(s(2:m)) + norm(alpha*delta_s_affine(2:m))) <= 0) && ((z(1) + alpha*delta_z_affine(1)) - (norm(z(2:m)) + norm(alpha*delta_z_affine(2:m))) <= 0) &&  ((kappa + alpha*delta_kappa_affine) <= 0) && ((tau + alpha*delta_tau_affine) <= 0)
            alpha = alpha - 0.1
        end
        
        %Centering parameter
        sigma = (1 - alpha)^3;

        %% 4) Combined (corrector and centering) direction Newton step
        d_x = (1-sigma)*r_x;
        d_y = (1-sigma)*r_y;
        d_z = (1-sigma)*r_z;
        d_tau = (1-sigma)*r_tau;
        d_s = soc_dot(lambda,lambda) + soc_dot(W'\delta_s_affine,W*delta_z_affine) - sigma*mu*e;
        d_kappa = kappa*tau + delta_kappa_affine*delta_tau_affine - sigma*mu;

        b_combined = [d_x; d_y; d_z-W'*(soc_dot_inv(lambda,d_s))];
        x_combined = KKT\b_combined;
        x2_combined = x_combined(1:length(x));
        y2_combined = x_combined((length(x)+1):(length(x)+length(y)));
        z2_combined = x_combined((length(x)+length(y)+1):(length(x)+length(y)+length(z)));

        numerator = (d_tau - (d_kappa/tau) + [c', b', h']*x_combined);
        denominator = ((kappa/tau) - [c', b', h']*x_static);
        delta_tau_combined = numerator/denominator;

        delta_x_combined = x2_combined + delta_tau_combined*x1;
        delta_y_combined = y2_combined + delta_tau_combined*y1;
        delta_z_combined = z2_combined + delta_tau_combined*z1;
        delta_s_combined = -W'*(soc_dot_inv(lambda,d_s) + W*delta_z_combined);
        delta_kappa_combined = -(d_kappa + kappa*delta_tau_combined)/tau;

        %% 5) Update iterates and scaling matrices
        alpha = 0.99;

        %Backtracking linesearch
        while ((s(1) + (alpha/0.99)*delta_s_combined(1)) - (norm(s(2:m)) + norm((alpha/0.99)*delta_s_combined(2:m))) <= 0) && ((z(1) + (alpha/0.99)*delta_z_combined(1)) - (norm(z(2:m)) + norm((alpha/0.99)*delta_z_combined(2:m))) <= 0) &&  ((kappa + (alpha/0.99)*delta_kappa_combined) <= 0) && ((tau + (alpha/0.99)*delta_tau_combined) <= 0)
           alpha = alpha - 0.1
        end
        
        %Updating optimization variables for use in next iteration
        update = [s; kappa; x; y; z; tau] + alpha*[delta_s_combined; delta_kappa_combined; delta_x_combined; delta_y_combined; delta_z_combined; delta_tau_combined];
        s = update(1:length(s));
        kappa = update((length(s)+1):(length(s)+length(kappa)));
        x = update((length(s)+length(kappa)+1):(length(s)+length(kappa)+length(x)));
        y = update((length(s)+length(kappa)+length(x)+1):(length(s)+length(kappa)+length(x)+length(y)));
        z = update((length(s)+length(kappa)+length(x)+length(y)+1):(length(s)+length(kappa)+length(x)+length(y)+length(z)));
        tau = update(end);

        %Updating NT-scaling matrix, W, and lambda
        W = w_soc(s,z,SOC_dimension);
        lambda = W*z;

        update_history = [update_history, update];
        mu_history = [mu_history, mu];
        gap_history = [gap_history, s'*z];
        equality_constraint_history = [equality_constraint_history, norm(A*x-b)];
        inequality_constraint_history = [inequality_constraint_history, norm(G*x+s-h)];
        dual_constraint_history = [dual_constraint_history, norm(c+A'*y+G'*z)];
        iter_axis = [iter_axis, iter];
    end
    toc

    x_opt = x;
    y_opt = y;
    s_opt = s;
    z_opt = z;
    primal_obj = c'*x; %Primal objective value
    dual_obj =  -h'*z - b'*y; %Dual objective value
    duality_gap = abs(primal_obj - dual_obj); %Duality gap
    
%     %% Plotting output
%     figure
%     stairs(iter_axis,mu_history)
%     grid on
%     title('\mu')
%     
%     figure
%     stairs(iter_axis,gap_history)
%     grid on
%     title('Primal-Dual Gap')
%     
%     figure
%     stairs(iter_axis,equality_constraint_history)
%     grid on
%     title('Equality Constraint Gap')
%     
%     figure
%     stairs(iter_axis,inequality_constraint_history)
%     grid on
%     title('Inequality Constraint Gap')
%     
%     figure
%     stairs(iter_axis,dual_constraint_history)
%     grid on
%     title('Dual Constraint Gap')
end