%%% STAT 310 Final Project
%%% Aaron Maurer

%%% Caller
setup = 0
prob = [1,1,1,1,1]

%%% Setup
if setup
    cd cvx
    cvx_setup
    cd ..
end

%%% Problem 1
if prob(1)
    % set seed and generate matricies/vectors
    rand('state',0)
    A = rand(300,100)
    c = rand(100,1)-1
    b = .5 * A * ones(100,1)

    % solve the linear programing problem to get optimal y and lower bound l
    cvx_begin
        variable y(100)
        minimize(c'*y)
        subject to
            A*y <= b
            0 <= y <= 1
    cvx_end
    l = cvx_optval

    % For values of tau going from 0 to 1 in increments of .01, round values of y above tau to 0 and above tau to 1 to getxtau (I store all 100 vectors in one matrix). Then, calculate objective value of c'x for each tau
    tau = .01*(1:100)                             % tau comes out as a 1x100 matrix
    xtau = repmat(y,1,100) >= repmat(tau,100,1) % xtau(r,c) is 1 if y(r)>=tau(c), and otherwise 0
    obj = c' * xtau

    % Test inequality constraint - generate the largest entry value of A*x -b for each tau
    con = max(A * xtau - repmat(b,1,100),[],1)

    % Pick the best objective value for which the constraint is satistifed and compare to our initial lower bound
    u = min(obj .* (con<=0))
    diff = u - l

end

