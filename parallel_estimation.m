simcount = 100; 
Tmax = 5e5; M = 8; H = 3;
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);

phi1 = @(x) ones(size(x)); phi2 = @(x) -double((x>=0)&(x<1/2))+double((x>=1/2)&(x<=1));
phi3 = @(x) sqrt(2)*phi2(2*x); phi4 = @(x) phi3(x-1/2);
phi5 = @(x) sqrt(2)*phi3(2*x); phi6 = @(x) phi5(x-1/4); phi7 = @(x) phi5(x-1/2); phi8 = @(x) phi5(x-3/4);
phi9 = @(x) sqrt(2)*phi5(2*x); phi10 = @(x) phi9(x-1/8); phi11 = @(x) phi9(x-2/8); phi12 = @(x) phi9(x-3/8);
phi13 = @(x) phi9(x-4/8); phi14 = @(x) phi9(x-5/8); phi15 = @(x) phi9(x-6/8); phi16 = @(x) phi9(x-7/8);
chi1 = @(t) exp(-t)*double(t>=0); chi2 = @(t) exp(-2*t)*double(t>=0); chi3 = @(t) exp(-3*t)*double(t>=0);
phi = {phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8}; chi = {chi1 chi2 chi3};

delete(gcp('nocreate'));
parpool(4);
parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; N = 256; Z = linspace(0,1,N+1); Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 1*Bfun1(x)*Bfun2(y)/2;

    Y = zeros(1,N+1);
    
    while realtime < T
        tau = exprnd(1/Lambda); 
        if ((Lambda-Lambda_0)*exp(-tau)+Lambda_0)/Lambda > rand
            lambda = (lambda - lambda_0)*exp(-tau) + lambda_0; Lambda = mean(lambda);
            U = rand; Yfunction = @(x) double(U<(sum(lambda(1:x))/(sum(lambda))));
            for i = [1:(N+1)]
                Y(i) = Yfunction(i);
            end
            [MM, I] = max(Y); X(end+1) = (I-1)/N; t(end+1) = realtime + tau; realtime = t(end);
            if lambda_0(I)/lambda(I) > rand
                Immigrants(end+1) = 1;
            else
                Immigrants(end+1) = 0;
            end
            for i = 1:(N+1)
                lambda(i) = lambda(i) + W((i-1)/N,X(end));
            end
            Lambda = mean(lambda);
        else
            realtime = realtime + tau;
            lambda = (lambda - lambda_0)*exp(-tau) + lambda_0; Lambda = mean(lambda);
        end
    end
    if t(end) > T
        t(end) = [];
        X(end) = [];
    end

    
    % Estimation 
    length_t = length(t); eta = log(T)/sqrt(T)/4;

    nt = length(Tgrid);
    psi = zeros(1+M*H, nt);
    psi(1,:) = ones(1,nt);
    event_index = 1;
    for i = 1:H
        decay = chi{i}(1/Q);
        for j = 1:M
            B = 0;
            row_index = 1 + (i-1)*M + j;
            for k_index = 1:nt
                B = B*decay; 
                kk = Tgrid(k_index);
                while event_index <= length_t && t(event_index) <= kk
                    B = B + phi{j}(X(event_index))*chi{i}(kk-t(event_index));
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        vartheta(end+1) = sum(phi{j}(X))/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi{i}(X(k))*psi(j,ceil(Q*t(k)));
            end
            vartheta(end+1) = B/T;
        end
    end
    
    G = zeros(M*H+1);
    for i = 1:(M*H+1)
        for j = 1:(M*H+1)
            G(i,j) = mean(psi(i,:).*psi(j,:));
        end
    end
    
    
    %theta = zeros(1,M+M*M+H);
    %alpha = theta(1:M);
    %beta = theta((M+1):(M+M*M));
    %beta = reshape(beta, M, M).';
    %gamma = theta((M+M*M+1):(M+M*M+H));
    %beta_otimes_gamma = zeros(M, M*H);
    %theta_otimes_k = zeros(M, M*H+1);
    %for i = [1:M]
    %    A = [];
    %    for k = [1:H]
    %        A(end+1:end+M) = beta(i,:).*gamma(k);
    %    end
    %    beta_otimes_gamma(i,:) = A;
    %    theta_otimes_k(i,:) = [alpha(i) A];
    %end
    %theta_otimes = [alpha reshape(beta_otimes_gamma.',1,M^2*H)];
    
    alpha0 = zeros(M,1);
    beta0  = zeros(M^2,1);
    gamma0 = ones(H,1)/H;
    u0 = max(gamma0,0);
    v0 = max(-gamma0,0);
    
    theta0 = [alpha0; beta0; u0; v0];
    lb = [-inf(M+M^2,1); zeros(2*H,1)];
    ub = [];
    
    
    Aeq = [zeros(1,M+M^2), ones(1,H), ones(1,H)];
    beq = 1;
    
    %eta = log(T)/sqrt(T)/4;
    obj = @(theta) penalizedLSloss(theta, vartheta, G, eta, M, H);
    nonlcon = @(theta) deal([], norm(theta(M+M*M+1:M+M*M+H) - theta(M+M*M+H+1:M+M*M+2*H),1) - 1);
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'Display', 'off', ...
        'SpecifyObjectiveGradient', false, ...   % true als je de gradient levert
        'MaxFunctionEvaluations', 1e8, ...
        'MaxIterations', 1e6);
    
    
    %theta_hat = fmincon(obj, theta0, [], [], Aeq, beq, lb, ub, ... @(theta) gamma_constraint(theta, M, H), options);
    theta_hat = fmincon(obj, theta0, [], [], [], [], lb, ub, ...
    @(theta) gamma_constraint(theta, M, H), options);
    alpha_hat = theta_hat(1:M);
    beta_hat = theta_hat((M+1):(M+M*M));
    u = theta_hat((M+M*M+1):(M+M*M+H)); v = theta_hat((M+M*M+H+1):(M+M*M+2*H));
    gamma_hat = u - v; signgammahat = sign(gamma_hat(1));
    gamma_hat = gamma_hat * signgammahat;
    norm(gamma_hat,1);   % = 1
    theta_hat = [alpha_hat' beta_hat' gamma_hat'];
    beta_hat = reshape(beta_hat, M, M).'*signgammahat;


    alpha_hat_sims(:,:,p) = alpha_hat;
    beta_hat_sims(:,:,p) = beta_hat;
    gamma_hat_sims(:,:,p) = gamma_hat;
    p
end
save('simulation_results.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');
%alpha_hat_sims
%beta_hat_sims
%gamma_hat_sims
%mean(alpha_hat_sims,3)     std(alpha_hat_sims,0,3)
%mean(beta_hat_sims,3)      std(beta_hat_sims,0,3)
%mean(gamma_hat_sims,3)     std(gamma_hat_sims,0,3)

%beta=zeros(8,8); beta(1,1)=3/4; beta(1,2)=1/4; beta(6,3)=3/20; beta(8,3)=-3/20;
%beta=zeros(8,8); beta(1,1)=4/5; beta(1,4)=1/4; beta(5:8,2)=-1/5;
%beta=zeros(8,8); beta(1,1)=3/4; beta(1,4)=1/6; beta(5:8,2)=-1/6; beta(5,5)=1/12; beta(7,7) = 1/12;
%beta=zeros(8,8); beta(1,1)=3/4; beta(1,2)=1/6; beta(1,3)=-1/6; beta(1,4)=1/6; beta(5:8,2:3)=-1/6; beta(5:8,5:8)=[-1;-1;1;1]*[1 -1 -1 -1]/12;
%MSEvector = zeros(1,simcount);
%for p = 1:simcount
    %A = beta - beta_hat_sims(:,:,p);
    %MSEvector(p) = norm(A, 'fro');
%end
%mean(MSEvector)
%std(MSEvector)
%histogram(MSEvector,10); xlabel('Frobenius norm $\|\hat\beta-\beta\|_F$', 'Interpreter', 'latex'); ylabel('frequency');
%title('Frobenius distance between estimated and true $\beta$', 'Interpreter', 'latex')

%FNRvec = zeros(1,simcount); FPRvec = zeros(1,simcount); F1vec = zeros(1,simcount);
%for p = 1:simcount
    %A = double(beta ~= 0); number1 = sum(A(:)); number0 = 64-number1; B = std(beta_hat_sims,0,3); C = double(abs(beta_hat_sims(:,:,p)) > 2*B); D = ones(8,8) - C; 
    %D0 = D .* A; FN = sum(D0(:)); FNR = FN/number1; C0 = C .* double(beta == 0); FP = sum(C0(:)); FPR = FP/number0;
    %TP = number1 - FN; F1 = 2*TP/(2*TP+FN+FP);
    %FNRvec(p) = FNR; FPRvec(p) = FPR; F1vec(p) = F1; 
%end
%mean(FNRvec); std(FNRvec); 
%mean(FPRvec); std(FPRvec); 
%mean(F1vec); std(F1vec); 



