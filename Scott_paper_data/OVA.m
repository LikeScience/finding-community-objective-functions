function [C_lb, C_ub, objBase, flagBase, sum_eps, C_lb_full, C_ub_full] = OVA(v, lb, ub, S, r_idx, c_fix, c_fix_idx)
    arguments
        v
        lb
        ub
        S
        r_idx
        c_fix = []
        c_fix_idx = []
    end
    x = v';
    m = size(S,1); 
    n = size(S,2); 
    N = size(v,2); 
    
    c_pos = optimvar('c_pos', n, 'LowerBound', 0,'UpperBound',1);
    c_neg = optimvar('c_neg', n, 'LowerBound', 0,'UpperBound',1);
    c = c_pos - c_neg;
    
    c_pos.UpperBound(c_fix_idx) = c_fix;
    c_neg.UpperBound(c_fix_idx) = c_fix;
    c_pos.LowerBound(c_fix_idx) = c_fix;
    c_neg.LowerBound(c_fix_idx) = c_fix;
    
    q1 = optimvar('q1',n,N,'LowerBound',0,'UpperBound',1000);
    q2 = optimvar('q2',n,N,'LowerBound',0,'UpperBound',1000);
    p = optimvar('p',m,N,'LowerBound',-1000,'UpperBound',1000);
    epsilon = optimvar('epsilon',N,'LowerBound',0,'UpperBound',1000);
    
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex','Display', 'none');
    
    prob1 = optimproblem('ObjectiveSense', 'minimize');
    prob1.Constraints.dual = p' * S - q1' + q2' == repmat(c', N, 1);
    prob1.Constraints.bounds = q2' * ub - q1' * lb - epsilon == x*c;
    prob1.Constraints.sum_c = sum(c) == 1;
    prob1.Objective = sum(epsilon);
    [sol1, fval1, exitflag1, ~] = solve(prob1,'Options', options);
    
    objBase = sol1.c_pos - sol1.c_neg;
    flagBase = exitflag1; 
    sum_eps = fval1;

    C_lb = zeros(length(r_idx), 1);
    C_ub = zeros(length(r_idx), 1);
    C_lb_full = zeros(n, length(r_idx));
    C_ub_full = zeros(n, length(r_idx));
    
    prob2 = optimproblem('ObjectiveSense', 'minimize');
    prob2.Constraints = prob1.Constraints;
    prob2.Constraints.eps_lock = epsilon == sol1.epsilon;
    
    for i = 1:length(r_idx)
        r = r_idx(i);
        
        prob2.Objective = c(r) + sum(c_pos + c_neg);
        sol_min = solve(prob2,'Options', options);
        C_lb(i) = sol_min.c_pos(r) - sol_min.c_neg(r);
        C_lb_full(:, i) = sol_min.c_pos - sol_min.c_neg;
        
        prob2.Objective = -c(r) + sum(c_pos + c_neg);
        sol_max = solve(prob2,'Options', options);
        C_ub(i) = sol_max.c_pos(r) - sol_max.c_neg(r);
        C_ub_full(:, i) = sol_max.c_pos - sol_max.c_neg;
    end
end