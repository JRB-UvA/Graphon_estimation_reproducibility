function [c, ceq] = gamma_constraint(theta, M, H)
    u = theta(M+M^2+1 : M+M^2+H);
    v = theta(M+M^2+H+1 : M+M^2+2*H);
    gamma = u - v;
    c = [];           
    ceq = sum(abs(gamma)) - 1; 
end