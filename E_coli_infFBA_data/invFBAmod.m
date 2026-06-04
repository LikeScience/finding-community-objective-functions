

function [objBase, sol_c_pos, sol_c_neg, flagBase, sum_eps] = invFBAmod(v, lb, ub, S, l0, useL2, c_fix, c_fix_idx, gurobi_path )
    arguments
        v
        lb
        ub
        S
        l0 = false
        useL2 = false
        c_fix = []
        c_fix_idx = []
        gurobi_path = 'C:\gurobi1301\win64\examples\matlab'
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

    % addpath(gurobi_path);
    % gurobi_setup;

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
    
    q1 = optimvar('q1',n,N,'LowerBound',0,'UpperBound',1000);
    q2 = optimvar('q2',n,N,'LowerBound',0,'UpperBound',1000);
    p = optimvar('p',m,N,'LowerBound',-1000,'UpperBound',1000);
    epsilon=optimvar('epsilon',N,'LowerBound',0,'UpperBound',1000);
    
    options = optimoptions('linprog', 'Algorithm', 'dual-simplex','Display', 'none');
    
    prob1 = optimproblem('ObjectiveSense', 'minimize');
    prob1.Constraints.dual = p' * S - q1' + q2' == repmat(c', N, 1);
    prob1.Constraints.bounds = q2' * ub - q1' * lb - epsilon == x*c;
    prob1.Constraints.sum_c = sum(c_pos)+sum(c_neg) == 1;
    prob1.Objective = sum(epsilon);
    [sol1, fval1, exitflag1, output1] = solve(prob1,'Options', options);
    objBase = sol1.c_pos - sol1.c_neg;
    sol_c_pos = sol1.c_pos;
    sol_c_neg = sol1.c_neg;
    flagBase = exitflag1; 
    sum_eps = fval1

 
end