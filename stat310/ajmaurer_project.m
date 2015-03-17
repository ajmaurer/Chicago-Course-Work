%%% STAT 310 Final Project
%%% Aaron Maurer

%%% Caller
setup = 1 
prob = [0,0,0,0,1]

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

    % For values of tau going from 0 to 1 in increments of .01, round values 
    % of y above tau to 0 and above tau to 1 to get xtau (I store all 100 vectors 
    % in one matrix). Then, calculate objective value of c'x for each tau
    tau = .01*(1:100)                     % tau comes out as a 1x100 matrix
    % xtau(r,c) is 1 if y(r)>=tau(c), and otherwise 0
    xtau = repmat(y,1,100) >= repmat(tau,100,1) 
    obj = c' * xtau

    % Test inequality constraint: get the largest entry value of A*x -b for each tau
    con = max(A * xtau - repmat(b,1,100),[],1)

    % Pick the best objective value for which the constraint is satistifed and 
    % compare to our initial lower bound
    % each entry of obj is negative, so zeroing out some entries doesn't matter
    u = min(obj .* (con<=0)) 
    diff = u - l

    % Plot objective value vs constraint
    plot(con,obj,'.k')
    hold on
    plot([0,0],[-45,0],'--r')
    title('Objective Value vs. Constraint')
    xlabel('Maximum Entry in Constraint Vector = max _i (Ax_\tau - b)_i')
    ylabel('Objective Value = c^Tx_\tau')
    legend('Point Corresponding to Particular \tau','Right Extreme of Feasible Region' ...
    ,'location','southwest')
    print('plot_1.pdf')
    hold off
    close
end

%%% Problem 2
if prob(2) 
    % input matricies
    a = [.1;.2;-.05;.1]
    A = a * a'     
    b = [.2;.1;.3;.1]   % known diagonal of X

    % Optimize cvx. The objective is a'Xa = trace(AX), and the constraints specify that
    % X is PSD and has the appropriate entries.
    cvx_begin
        variable X(4,4) symmetric
        minimize(trace(A*X))
        subject to
            X == semidefinite(4)
            diag(X) == b 
            X(1,2) >= 0
            X(1,3) >= 0
            X(2,3) <= 0
            X(2,4) <= 0
            X(3,4) >= 0
    cvx_end
end

%%% Problem 3
if prob(3) 

    % Get Matricies.
    newimage        % Given program

    % Convert matricies so this can be formatted as a QP/LP. Matricies huge but sparse
    % For both problems, the objective is a norm of Cx
    xs = sparse(reshape(S.*X,[m*n,1]))      % elements of X, or 0 if data missing
    SD = sparse(diag(reshape(S,[m*n,1])))   % has elements of S along diagonal
    CT = speye([m*(n-1),m*n]) - [sparse(m*(n-1),m),speye(m*(n-1),m*(n-1))] % Top of C
    CB = speye([(m-1)*n,m*n]) - [sparse((m-1)*n,1),speye((m-1)*n,m*n-1)] % bottom C
    C  = [CT;CB]

    % Show original image
    colormap gray
    imagesc(X)
    print('image_2_orig.pdf')
    close

    % Show sparse image
    colormap gray
    imagesc(abs(X-127*(1-S)))
    print('image_2_miss.pdf')
    close

    %% part a: L2 variation. 

    % optimize the problem, in this case a QP since were minimizing the 2 norm
    % SD*x == xs is asserting x_ij = a_ij when S=1, and is 0=0 when S=0 (always true)
    cvx_begin
        variable x(m*n)
        minimize(norm(C*x,2))
        subject to 
            SD*x == xs 
    cvx_end

    % output the image
    colormap gray
    imagesc(reshape(x,[m,n]))
    print('image_2a.pdf')
    close
    cvx_optval_a = cvx_optval

    %% part b: L1 variation. 

    % optimize the problem, in this case a LP
    cvx_begin
        variable x(m*n)
        minimize(norm(C*x,1))
        subject to 
            SD*x == xs 
    cvx_end
    cvx_optval_b = cvx_optval

    % output the image
    colormap gray
    imagesc(reshape(x,[m,n]))
    print('image_2b.pdf')
    close
end

%%% Problem 4
if prob(4) 

    % get the data.
    censoring           % given program
    AL = A(:,1:k)       % left 25 columns A
    AR = A(:,(k+1):n)   % remainder of A 

    %% part a: fit a constrained QP. 
    % Usual least squares for the first k data points, 
    % constraint that a'x>beta for remainder
    cvx_begin
        variable x_qp(d)
        minimize(norm(AL'*x_qp-b,2))
        subject to
            AR'*x_qp>=beta
    cvx_end

    %% part b: usual least squares
    x_ls = (AL*AL')^-1 * AL*b 

    %% part c: comparison
    re_qp = norm(x_qp - x_true)/norm(x_true)
    re_ls = norm(x_ls - x_true)/norm(x_true)
    
end

%%% Problem 5
if prob(5)

    % get the data
    sparsify        % given program
    Ey = c'*mu      % expectation of c'y
    Vy = c'*Sigma*c % variance of c'y 
    cS = c'*Sigma   % when multiplied by x gives cov(c'a,y'a)

    % Run the QCLP 
    % we are minimizing L1 of x. The constraints are turned into QC using the identities 
    % E(y)^2=E(y^2)+var(y) and MSE(y-hy) = E(y-hy)^2 + var(y-hy)
    cvx_begin
        variable x(n)
        minimize(norm(x,1))
        subject to
            .99*Ey^2 - 2*Ey*x'*mu + (x'*mu)^2 + .99*Vy - 2*cS*x + x'*Sigma*x<=0
    cvx_end

    % now, calculate the number of entries of x above various thresholds:
    xpos=zeros(6,1)
    for i=1:6 
        xpos(i)=sum(x>10^-i) 
    end
end




