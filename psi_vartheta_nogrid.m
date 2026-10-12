psi_at_events = zeros(length_t, 1 + M*H);
for i = 1:H
    for j = 1:M
        row_index = 1 + (i-1)*M + j;
        psi_val = 0;
        for k = 1:length_t
            if k == 1
                dt = t(1);
            else
                dt = t(k) - t(k-1);
            end
            psi_val = psi_val * chi{i}(dt);
            psi_at_events(k, row_index) = psi_val;
            psi_val = psi_val + phi{j}(X(k));
        end
    end
end

vartheta = zeros(M + M^2*H, 1);
counter = 1;

for j = 1:M
    vartheta(counter) = sum(phi{j}(X))/T;
    counter = counter + 1;
end
for i = 1:M
    for j = 1:M
        for h = 1:H
            s = 0;
            psi_index = 1 + (h-1)*M + j;
            for k = 1:length_t
                s = s + phi{i}(X(k)) * psi_at_events(k, psi_index);
            end
            vartheta(counter) = s/T;
            counter = counter + 1;
        end
    end
end

G = zeros(1 + M*H);
for i = 1:(1 + M*H)
    for j = 1:(1 + M*H)
        s = 0;
        for k = 1:length_t
            s = s + psi_at_events(k,i) * psi_at_events(k,j);
        end
        G(i,j) = s/(length_t);
    end
end