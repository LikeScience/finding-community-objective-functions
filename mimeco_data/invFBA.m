

function [objL1, flagL1, objL0, flagL0, objBase, flagBase, sum_eps] = invFBA(v, lb, ub, S, l0, useL2, c_fix, c_fix_idx, min_obj, lower_bound_obj)
    arguments
        v
        lb
        ub
        S
        l0 = false
        useL2 = false
        c_fix = []
        c_fix_idx = []
        min_obj = 0
        lower_bound_obj = -1000
    end
    % invFBA: invFBA algorithm to determine an objective vector from flux distribution data. 
    % The algorithm is as described by Zhao et al., 2016, but the
    % implementation uses 
    % 
    % Inputs: 
    % v = flux vector (num_rnxs x num_samples)
    % lb = lower bounds on fluxes (num_rnxs)
    % ub = upper bounds on fluxes (num_rnxs)
    % S = stoichiometric matrix (num_mets x num_rnxs)
    % l0 = boolean indicating if the final l0 regularization should be done (optional)
    % c_fix = values of objective vector elements to fix (optional)
    % c_fix_idx = indices of objective vector elements to fix (optional)
    %
    % Outputs
    % objL1 = Objective vector after linear invFBA and L1 regularization
    % objL0 = Objective vector after linear invFBA, L1 and L0 regularization

    x = v';
    m = size(S,1); %num_mets
    n = size(S,2); %num_rnxs
    N = size(v,2); %num_samples
    
    c_pos = optimvar('c_pos', n, 'LowerBound', 0,'UpperBound',1);
    c_neg = optimvar('c_neg', n, 'LowerBound', 0,'UpperBound',1);
    c = c_pos - c_neg;
    
    c_pos.UpperBound(c_fix_idx) = c_fix;
    c_neg.UpperBound(c_fix_idx) = c_fix;
    c_pos.LowerBound(c_fix_idx) = c_fix;
    c_neg.LowerBound(c_fix_idx) = c_fix;
    if lower_bound_obj ~= -1000
        c_pos.LowerBound(setdiff(1:n, c_fix_idx)) = lower_bound_obj;
        c_neg.LowerBound(setdiff(1:n, c_fix_idx)) = lower_bound_obj;
    end
    
    q1 = optimvar('q1',n,N,'LowerBound',0,'UpperBound',1000);
    q2 = optimvar('q2',n,N,'LowerBound',0,'UpperBound',1000);
    p = optimvar('p',m,N,'LowerBound',-1000,'UpperBound',1000);
    epsilon=optimvar('epsilon',N,'LowerBound',0,'UpperBound',1000);
    
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex','Display', 'none');
    prob1 = optimproblem('ObjectiveSense', 'minimize');
    prob1.Constraints.dual = p' * S - q1' + q2' == repmat(c', N, 1);
    prob1.Constraints.bounds = q2' * ub - q1' * lb - epsilon == x*c;
    prob1.Constraints.sum_c = sum(c) == 1;
    if min_obj > 0
        prob1.Constraints.min_obj = x*c >= repmat(min_obj,N,1);
    end
    prob1.Objective = sum(epsilon);
    [sol1, fval1, exitflag1, output1] = solve(prob1,'Options', options);
    objBase = sol1.c_pos - sol1.c_neg;
    flagBase = exitflag1; 
    sum_eps = fval1

    prob2 = optimproblem('ObjectiveSense', 'minimize');
    prob2.Constraints = prob1.Constraints;
    prob2.Constraints.eps_lock = epsilon == sol1.epsilon;
    if useL2 == true
       prob2.Objective = sum((c_pos - c_neg).^2);
    else
       prob2.Objective = sum(c_pos + c_neg); 
    end
    [sol2, fval2, exitflag2, output2] = solve(prob2,'Options', options);
    disp(output2.message)
    objL1 = sol2.c_pos - sol2.c_neg;
    flagL1 = exitflag2; 
    objL0 = []; flagL0 = [];
    if l0 
        prob3 = optimproblem('ObjectiveSense', 'minimize');
        options = optimoptions('intlinprog', 'Display', 'none');
        prob3.Constraints = prob1.Constraints;
        prob3.Constraints.eps_lock = epsilon == sol1.epsilon;
        z = optimvar('z', n, 'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);
        prob3.Constraints.z_link = c_pos + c_neg <= z;
        prob3.Objective = sum(z);
        init.c_pos = sol2.c_pos;
        init.c_neg = sol2.c_neg;
        init.p = sol2.p;
        init.q1 = sol2.q1;
        init.q2 = sol2.q2;
        init.epsilon = sol2.epsilon;
        init.z = double((sol2.c_pos + sol2.c_neg) > 1e-7);
        [sol3, fval3, exitflag3, output3] = solve(prob3,init,'Options', options);
        objL0 = sol3.c_pos - sol3.c_neg;
        flagL0 = exitflag3;
    end
end