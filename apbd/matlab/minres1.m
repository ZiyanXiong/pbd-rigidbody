%% Minimum Residual (GNU Octave)
% https://savannah.gnu.org/patch/?9282
% CLEANED
function [x,flag,relres,iter,resvec]  = minres1(A,b,tol,maxit,M,~,x0,~,~,~,mindices,mu)
flag = 1;

n = length(b);
beta_p = 0;
c_pp = 0;
c_p = 0;
s_pp = 0;
s_p = 0;

% Initiation
rM = b - A * x0;
resvec = norm(rM);
MrM = M \ rM;

beta = sqrt(rM' * MrM);
normb = norm(b);
relres = resvec / normb;

v_p = zeros(n, 1);
v_n = rM;
temp2 = beta;
m_p = zeros(n, 1);
m_pp = zeros(n, 1);
Am_p = zeros(n, 1);
Am_pp = zeros(n, 1);
x = x0;

Mv_n = M \ v_n;

% Iteration
for iter = 1 : maxit

    AMv_n = A * Mv_n;

    alpha = Mv_n' * (AMv_n / beta / beta);
    if iter == 1
        v_f = (AMv_n - alpha * v_n) / beta;
    else
        v_f = (AMv_n - alpha * v_n) / beta - v_p * beta / beta_p;
    end

    Mv_f = M \ v_f;

    beta_n = sqrt(v_f' * Mv_f);

    if iter == 1
        epsilon = 0;
        gamma_h = alpha;
        delta = 0;
    elseif iter == 2
        epsilon = 0;
        gamma_h = beta * s_p - alpha * c_p;
        delta = beta * c_p + alpha * s_p;
    else
        epsilon = s_pp * beta;
        gamma_h = - c_pp * beta * s_p - alpha * c_p;
        delta = - c_pp * beta * c_p + alpha * s_p;
    end

    gamma = sqrt(gamma_h ^ 2 + beta_n ^ 2);
    c = gamma_h / gamma;
    s = beta_n / gamma;

    m = (Mv_n / beta - epsilon * m_pp - delta * m_p) / gamma;
    Am = (AMv_n / beta - epsilon * Am_pp - delta * Am_p) / gamma;

    x = x + m * temp2 * c;
    rM = rM - Am * temp2 * c;

    normrM = norm(rM);
    temp2 = temp2 * s;

    relres = normrM / normb;
    resvec(end+1) = normrM; %#ok<AGROW>
    
    % Check convergence
    feasible = true;
    for i = 1:length(x)
        if(mod(mindices(i),3)==1 && x(i) < -1e-6)
            feasible = false;
        end

        if(mod(mindices(i),3)==2 && norm([x(i) x(i+1)]) > mu*x(i-1) + 1e-6)
            feasible = false;
        end
    end

    if(~feasible)
        flag = 1;
        break;
    end

    if (relres <= tol) || (beta_n <= eps)
        flag = 0;
        break
    end

    % Check stagnation
    if norm(resvec(end)-resvec(end-1)) <= eps * norm(resvec(end))
        flag = 3;
        break
    end

    beta_p = beta;
    beta = beta_n;
    c_pp = c_p;
    c_p = c;
    s_pp = s_p;
    s_p = s;
    m_pp = m_p;
    m_p = m;
    Am_pp = Am_p;
    Am_p = Am;
    v_p = v_n;
    v_n = v_f;
    Mv_n = Mv_f;
end

end