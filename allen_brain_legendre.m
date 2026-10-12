%ecephys session 715093703
load('Allen_neuropixel.mat'); X = (X - 0.003)/.994; %t = t(1:5e5); X = X(1:5e5);
Tmax = max(t); T = Tmax;  

freq = length(t)/Tmax;
M = 16; H = 5;
phi = cell(16,1);
for n = 0:15
    phi{n+1} = @(x) sqrt(2*n+1) * legendreP(n, 2*x-1);
end
phi1  = phi{1}; phi2  = phi{2}; phi3  = phi{3}; phi4  = phi{4};phi5  = phi{5}; phi6  = phi{6}; phi7  = phi{7}; phi8  = phi{8};
phi9  = phi{9}; phi10 = phi{10}; phi11 = phi{11}; phi12 = phi{12};phi13 = phi{13}; phi14 = phi{14}; phi15 = phi{15}; phi16 = phi{16};
chi1 = @(t) exp(-freq*t)*double(t>=0); chi2 = @(t) exp(-freq*2*t)*double(t>=0); chi3 = @(t) exp(-freq*5*t)*double(t>=0);
chi4 = @(t) exp(-freq*t/2)*double(t>=0); chi5 = @(t) exp(-freq*t/5)*double(t>=0);
phi = {phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 phi9 phi10 phi11 phi12 phi13 phi14 phi15 phi16}; chi = {chi1 chi2 chi3 chi4 chi5};
N = 1024; Z = linspace(0,1,N+1); Wgrid = zeros(N+1, N+1);
phi1_vec = phi1(Z); phi2_vec = phi2(Z); phi3_vec = phi3(Z); phi4_vec = phi4(Z); phi5_vec = phi5(Z); phi6_vec = phi6(Z); phi7_vec = phi7(Z); phi8_vec = phi8(Z);
phi9_vec = phi9(Z); phi10_vec = phi10(Z); phi11_vec = phi11(Z); phi12_vec = phi12(Z); phi13_vec = phi13(Z); phi14_vec = phi14(Z); phi15_vec = phi15(Z); phi16_vec = phi16(Z);
phi_mat = [phi1_vec;phi2_vec;phi3_vec;phi4_vec;phi5_vec;phi6_vec;phi7_vec;phi8_vec;phi9_vec;phi10_vec;phi11_vec;phi12_vec;phi13_vec;phi14_vec;phi15_vec;phi16_vec];

 
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
                B = B + phi_mat(j,ceil(X(event_index)*N))*chi{i}(kk-t(event_index)); 
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
        B = B + phi_mat(j,ceil(X(index)*N));
    end
    vartheta(end+1) = B/T;
end
for i = 1:M
    for j = 2:(M*H+1)
        B = 0;
        for k = 1:length_t
            B = B + phi_mat(i,ceil(X(k)*N))*psi(j,ceil(Q*t(k))); 
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


save('data_analysis_penalized_leg.mat', 'alpha_hat', 'beta_hat', 'gamma_hat');




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

save('data_analysis_nonpenalized_leg.mat', 'alpha_hat', 'beta_hat', 'gamma_hat');




%------------------------------------------------------------

eta = log(T)/sqrt(T)/2;

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

save('data_analysis_penalized2_leg.mat', 'alpha_hat', 'beta_hat', 'gamma_hat');




%------------------------------------------------------------

eta = log(T)/sqrt(T)/1;

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

save('data_analysis_penalized3_leg.mat', 'alpha_hat', 'beta_hat', 'gamma_hat');




%------------------------------------------------------------




%baseline = @(x) baseline_fun(x, alpha_hat, phi); baseline_normalized = @(x) baseline(x) / freq;

%graphon = @(x,y) graphon_fun(x,y,beta_hat,phi); graphon_normalized = @(x,y) graphon(x,y) / freq;

%x = linspace(0,1,250); y = baseline_normalized(x);
%plot(x,y)
%xlabel('Space $x$','Interpreter', 'latex')
%ylabel('Normalized intensity $\nu(x)/f$','Interpreter', 'latex')
%title('Normalized estimated baseline using penalized LS')

%[x,y] = meshgrid(linspace(0,1,70), linspace(0,1,70));
%Z = graphon_normalized(x,y);
%surf(x,y,Z)
%xlabel('Space $x$','Interpreter', 'latex')
%ylabel('Space $y$','Interpreter', 'latex')
%zlabel('Normalized graphon $W(x,y)/f$','Interpreter', 'latex')
%title('Normalized estimated graphon using penalized LS')

