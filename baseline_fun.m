function val = baseline_fun(x, alpha_hat, phi)
    val = zeros(size(x));  % als phi{i}(x) vector teruggeeft
    for i = 1:length(alpha_hat)
        val = val + alpha_hat(i) * phi{i}(x);
    end
end