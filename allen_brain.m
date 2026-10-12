%ecephys session 715093703
load('Allen_neuropixel.mat'); X = (X - 0.003)/.994; %t = t(1:5e5); X = X(1:5e5);
Tmax = max(t); T = Tmax;  

freq = length(t)/Tmax;
M = 16; H = 5;
phi1 = @(x) ones(size(x)); phi2 = @(x) -double((x>=0)&(x<1/2))+double((x>=1/2)&(x<=1));
phi3 = @(x) sqrt(2)*phi2(2*x); phi4 = @(x) phi3(x-1/2);
phi5 = @(x) sqrt(2)*phi3(2*x); phi6 = @(x) phi5(x-1/4); phi7 = @(x) phi5(x-1/2); phi8 = @(x) phi5(x-3/4);
phi9 = @(x) sqrt(2)*phi5(2*x); phi10 = @(x) phi9(x-1/8); phi11 = @(x) phi9(x-2/8); phi12 = @(x) phi9(x-3/8);
phi13 = @(x) phi9(x-4/8); phi14 = @(x) phi9(x-5/8); phi15 = @(x) phi9(x-6/8); phi16 = @(x) phi9(x-7/8);
chi1 = @(t) exp(-freq*t)*double(t>=0); chi2 = @(t) exp(-freq*2*t)*double(t>=0); chi3 = @(t) exp(-freq*3*t)*double(t>=0);
chi4 = @(t) exp(-freq*t/2)*double(t>=0); chi5 = @(t) exp(-freq*t/3)*double(t>=0);
phi = {phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 phi9 phi10 phi11 phi12 phi13 phi14 phi15 phi16}; 
chi = {chi1 chi2 chi3 chi4 chi5};

 
length_t = length(t); eta = log(T)/sqrt(T)/4; res = 2e7; Tgrid = linspace(0,T,res+1); Q = res/length_t;

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


save('data_analysis_penalized.mat', 'alpha_hat', 'beta_hat', 'gamma_hat');




% --------------------------
eta = 0;

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

save('data_analysis_nonpenalized.mat', 'alpha_hat', 'beta_hat', 'gamma_hat');




%------------------------------------------------------------




%baseline = @(x) baseline_fun(x, alpha_hat, phi); baseline_normalized = @(x) baseline(x) / freq;

%graphon = @(x,y) graphon_fun(x,y,beta_hat,phi); graphon_normalized = @(x,y) graphon(x,y) / freq;

%x = linspace(0,1,1000); y = baseline_normalized(x);
%plot(x,y)
%xlabel('Space $x$','Interpreter', 'latex')
%title('Frobenius distance between estimated and true $\beta$', 'Interpreter', 'latex')
%ylabel('Normalized intensity $\nu(x)/f$','Interpreter', 'latex')
%title('Normalized estimated baseline using penalized LS')

%[x,y] = meshgrid(linspace(0,1,110), linspace(0,1,110));
%Z = graphon_normalized(x,y);
%surf(x,y,Z)
%xlabel('Space $x$','Interpreter', 'latex')
%ylabel('Space $y$','Interpreter', 'latex')
%zlabel('Normalized graphon $W(x,y)/f$','Interpreter', 'latex')
%title('Normalized estimated graphon using penalized LS')

