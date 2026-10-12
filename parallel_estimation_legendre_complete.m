simcount = 100; 
Tmax = 5e5; M = 8; H = 3;
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);

phi = cell(16,1);
for n = 0:15
    phi{n+1} = @(x) sqrt(2*n+1) * legendreP(n, 2*x-1);
end
phi1  = phi{1}; phi2  = phi{2}; phi3  = phi{3}; phi4  = phi{4};phi5  = phi{5}; phi6  = phi{6}; phi7  = phi{7}; phi8  = phi{8};
phi9  = phi{9}; phi10 = phi{10}; phi11 = phi{11}; phi12 = phi{12};phi13 = phi{13}; phi14 = phi{14}; phi15 = phi{15}; phi16 = phi{16};
chi1 = @(t) exp(-t)*double(t>=0); chi2 = @(t) exp(-2*t)*double(t>=0); chi3 = @(t) exp(-3*t)*double(t>=0);
phi = {phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8}; chi = {chi1 chi2 chi3};
N = 1024; Z = linspace(0,1,N+1); Wgrid = zeros(N+1, N+1);
phi1_vec = phi1(Z); phi2_vec = phi2(Z); phi3_vec = phi3(Z); phi4_vec = phi4(Z); phi5_vec = phi5(Z); phi6_vec = phi6(Z); phi7_vec = phi7(Z); phi8_vec = phi8(Z);
phi_mat = [phi1_vec;phi2_vec;phi3_vec;phi4_vec;phi5_vec;phi6_vec;phi7_vec;phi8_vec];
phi9_vec = phi9(Z); phi10_vec = phi10(Z); phi11_vec = phi11(Z); phi12_vec = phi12(Z); phi13_vec = phi13(Z); phi14_vec = phi14(Z); phi15_vec = phi15(Z); phi16_vec = phi16(Z);
phi_mat2 = [phi9_vec;phi10_vec;phi11_vec;phi12_vec;phi13_vec;phi14_vec;phi15_vec;phi16_vec];
for i = 1:(N+1)
    for j = 1:(N+1)
        Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        %Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        %Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        %Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        %Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

%delete(gcp('nocreate'));
%parpool(4);

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_un1.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);
Wgrid = zeros(N+1, N+1); 
for i = 1:(N+1)
    for j = 1:(N+1)
        %Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        %Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        %Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        %Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_un2.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);
Wgrid = zeros(N+1, N+1); 
for i = 1:(N+1)
    for j = 1:(N+1)
        %Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        %Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        %Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        %Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        %Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_un3.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);
Wgrid = zeros(N+1, N+1); 
for i = 1:(N+1)
    for j = 1:(N+1)
        %Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        %Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        %Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        %Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        %Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_un4.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);
Wgrid = zeros(N+1, N+1); 
for i = 1:(N+1)
    for j = 1:(N+1)
        Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        %Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        %Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        %Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_per1.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);
Wgrid = zeros(N+1, N+1); 
for i = 1:(N+1)
    for j = 1:(N+1)
        Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        %Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        %Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        %Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_per2.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_hat_sims = zeros(1,M,simcount);
beta_hat_sims = zeros(M,M,simcount);
gamma_hat_sims = zeros(H,1,simcount);
Wgrid = zeros(N+1, N+1); 
for i = 1:(N+1)
    for j = 1:(N+1)
        Wgrid(i,j) = 3*phi1_vec(i)*phi1_vec(j)/4 + 5/20*phi1_vec(i)*phi2_vec(j) + 3/20*phi6_vec(i)*phi3_vec(j) - 3/20*phi8_vec(i)*phi3_vec(j);
        %Wgrid(i,j) = 7/10*phi1_vec(i)*phi1_vec(j)+1/4*phi1_vec(i)*phi4_vec(j)-(phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/5;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/6*phi1_vec(i)*phi4_vec(j) - (phi2_vec(j)*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/6 + (phi5_vec(i)*phi5_vec(j)+phi7_vec(i)*phi7_vec(j))/12;
        %Wgrid(i,j) = 3/4*phi1_vec(i)*phi1_vec(j) + 1/12*phi1_vec(i)*(phi2_vec(j)-phi3_vec(j)+phi4_vec(j)) - ((phi2_vec(j)+phi3_vec(j))*(phi5_vec(i)+phi6_vec(i)+phi7_vec(i)+phi8_vec(i)))/12 + ((phi5_vec(i)+phi6_vec(i)-phi7_vec(i)-phi8_vec(i))*(-phi5_vec(j)+phi6_vec(j)+phi7_vec(j)+phi8_vec(j)))/12;
       
        % (A) add next basis component
        %Wgrid(i,j) = Wgrid(i,j) + phi9_vec(i)*(phi9_vec(j)+phi10_vec(j)+phi11_vec(j)+phi12_vec(j)+phi13_vec(j)+phi14_vec(j)+phi15_vec(j)+phi16_vec(j))/3;
        % (B) add periodic disturbation 
        %Wgrid(i,j) = Wgrid(i,j) + 2*sin(20*Z(i))*sin(30*Z(j))/3;
        % (C) add rough (Brownian) disturbation 
        dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
        dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
        Wgrid(i,j) = Wgrid(i,j) + Bfun1(Z(i))*Bfun2(Z(j))/3;
    end
end

parfor p = 1:simcount
    rng(74973836 + p, 'philox');
    T=Tmax; Q = 10; Tgrid = linspace(0,T,Q*T+1);

    % simulations 
    X = []; t = []; T = Tmax; eps = 0.05; realtime = 0;
    %Tgrid = linspace(0,T,N*T/10+1);
    lambda_0function = @(x) ones(size(x))/20; Immigrants = []; %lambda_0function = @(x) 0.5*double((x<eps+1/2)&(x>1/2-eps));
    lambda_0 = lambda_0function(Z);
    lambda = lambda_0; Lambda_0 = mean(lambda_0); Lambda = Lambda_0;
    timegrid = linspace(0,T,10*N+1);

    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y);
    %W = @(x,y) 7*phi1(x)*phi1(y)/10 + 1*phi1(x)*phi4(y)/4 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/5 ;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*phi4(y)/6 - phi2(y)*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)*phi5(y)+phi7(x)*phi7(y))/12;
    %W = @(x,y) 3*phi1(x)*phi1(y)/4 + 1*phi1(x)*(phi2(y)-phi3(y)+phi4(y))/6 - (phi2(y)+phi3(y))*(phi5(x)+phi6(x)+phi7(x)+phi8(x))/6 + (phi5(x)+phi6(x)-phi7(x)-phi8(x))*(-phi5(y)+phi6(y)+phi7(y)+phi8(y))/12;

    % (A) add next basis component
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + (phi9(x)*(phi9(y)+phi10(y)+phi11(y)+phi12(y)+phi13(y)+phi14(y)+phi15(y)+phi16(y)))/3;
    
    % (B) add periodic disturbation 
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + 2*sin(20*x)*sin(30*y)/3;
    
    % (C) add rough (Brownian) disturbation 
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun1 = @(t) interp1(Z, B, t, 'linear');
    %dx = 1/N; dB = sqrt(dx) * randn(1,N); B = [0, cumsum(dB)]; Bfun2 = @(t) interp1(Z, B, t, 'linear');
    %W = @(x,y) 3*phi1(x)*phi1(y)/4+5/20*phi1(x)*phi2(y)+3/20*phi6(x)*phi3(y)-3/20*phi8(x)*phi3(y) + Bfun1(x)*Bfun2(y)/6;

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
            for i = [1:(N+1)]
                lambda(i) = lambda(i) + Wgrid(i,X(end)*N+1);
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
                    B = B + phi_mat(j,X(event_index)*N+1)*chi{i}(kk-t(event_index)); 
                    event_index = event_index + 1;
                end
                psi(row_index, k_index) = B;
            end
            event_index = 1;
        end
    end
    
        
    vartheta = []; 
    for j = 1:M
        B = 0;
        for index = 1:length_t
            B = B + phi_mat(j,X(index)*N+1);
        end
        vartheta(end+1) = B/T;
    end
    for i = 1:M
        for j = 2:(M*H+1)
            B = 0;
            for k = 1:length_t
                B = B + phi_mat(i,X(k)*N+1)*psi(j,ceil(Q*t(k))); 
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
save('simulation_results_per3.mat', 'alpha_hat_sims', 'beta_hat_sims', 'gamma_hat_sims');

