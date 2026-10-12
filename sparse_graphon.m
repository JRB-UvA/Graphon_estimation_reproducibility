phi1 = @(x) ones(size(x)); phi2 = @(x) -double((x>=0)&(x<1/2))+double((x>=1/2)&(x<=1));
phi3 = @(x) sqrt(2)*phi2(2*x); phi4 = @(x) phi3(x-1/2);
phi5 = @(x) sqrt(2)*phi3(2*x); phi6 = @(x) phi5(x-1/4); phi7 = @(x) phi5(x-1/2); phi8 = @(x) phi5(x-3/4);
chi1 = @(t) exp(-t)*double(t>=0); chi2 = @(t) exp(-2*t)*double(t>=0); chi3 = @(t) exp(-3*t)*double(t>=0);
phi = {phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8}; chi = {chi1 chi2 chi3};
M = 8; H = 3; 

N = 500; Z = linspace(0,1,N+1); T = 1e4; Q = 25; Tgrid = linspace(0,T,Q*T+1); 
%X = []; t = []; 
length_t = length(t); eta = log(T)/sqrt(T)/4; cutoff = ceil(50*max(length_t/T,1));

clear psi
counter = 0;
psi(1,:) = [phi1(Tgrid)]; 
for i = [1:H]
    for j = [1:M]
        A = [];
        for k = Tgrid
            B = 0;
            counter = max(counter - cutoff,1);
            while counter <= length_t && t(counter) < k
                    B = B + phi{j}(X(counter))*chi{i}(k-t(counter));    
                    counter = counter + 1;
            end
            A(end+1) = B;
        end
        psi(end+1,:) = A;
    end
end

1

vartheta = []; 
for j = [1:M]
    vartheta(end+1) = sum(phi{j}(X))/T;
end
for i = [1:M]
    for j = [2:(M*H+1)]
        B = 0;
        for k = [1:length_t]
            B = B + phi{i}(X(k))*psi(j,ceil(Q*t(k)));
        end
        vartheta(end+1) = B/T;
    end
end

G = zeros(M*H+1);
for i = [1:(M*H+1)]
    for j = [1:(M*H+1)]
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
    'Display', 'iter', ...
    'SpecifyObjectiveGradient', false, ...   % true als je de gradient levert
    'MaxFunctionEvaluations', 1e6, ...
    'MaxIterations', 1e4);


%theta_hat = fmincon(obj, theta0, [], [], Aeq, beq, lb, ub, ... @(theta) gamma_constraint(theta, M, H), options);
theta_hat = fmincon(obj, theta0, [], [], [], [], lb, ub, ...
@(theta) gamma_constraint(theta, M, H), options);
alpha_hat = theta_hat(1:M);
beta_hat = theta_hat((M+1):(M+M*M));
u = theta_hat((M+M*M+1):(M+M*M+H)); v = theta_hat((M+M*M+H+1):(M+M*M+2*H));
gamma_hat = u - v; signgammahat = signsign(gamma_hat(1));
gamma_hat = gamma_hat * signgammahat;
norm(gamma_hat,1);   % = 1
theta_hat = [alpha_hat' beta_hat' gamma_hat'];
beta_hat = reshape(beta_hat, M, M).'*sign(gamma_hat(1));
%gamma = theta((M+M*M+1):(M+M*M+H));


