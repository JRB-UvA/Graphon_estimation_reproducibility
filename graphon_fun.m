function val = graphon_fun(x, y, beta_hat, phi)
    val = zeros(size(x)); 
    for k = 1:size(beta_hat,1)
        phi_k_x = phi{k}(x); 
        for l = 1:size(beta_hat,2)
            phi_l_y = phi{l}(y); 
            val = val + beta_hat(k,l) * (phi_k_x .* phi_l_y);
        end
    end
end