function loss = penalizedLSloss(theta, vartheta, G, eta, M, H)
    alpha = theta(1:M);
    beta = theta((M+1):(M+M*M));
    beta = reshape(beta, M, M).';
    %gamma = theta((M+M*M+1):(M+M*M+H));
    u = theta((M+M*M+1):(M+M*M+H)); v = theta((M+M*M+H+1):(M+M*M+2*H));
    gamma = u - v;
    beta_otimes_gamma = zeros(M, M*H);
    theta_otimes_k = zeros(M, M*H+1);
    for i = [1:M]
        A = [];
        for k = [1:H]
            A(end+1:end+M) = beta(i,:).*gamma(k);
        end
        beta_otimes_gamma(i,:) = A;
        theta_otimes_k(i,:) = [alpha(i) A];
    end
    theta_otimes = [alpha' reshape(beta_otimes_gamma.',1,M^2*H)];

    qfG = 0;
    for i = [1:M]
        qfG = qfG + theta_otimes_k(i,:)*G*theta_otimes_k(i,:)';
    end
    loss = -theta_otimes * vartheta' + qfG/2 + eta*sum(abs(beta_otimes_gamma(:)));
end