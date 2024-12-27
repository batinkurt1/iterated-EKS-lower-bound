function Ainv = calculate_Ainv(A,nx)
    n = length(A)/nx;
    Ainv = zeros(n*nx);

    delta = cell(1,n);
    sigma = cell(1,n);

    for i = 1:n
        k = (n-i)+1;

        if i == 1
            D1 = A((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx);
            Dn = A((k-1)*nx+1:k*nx,(k-1)*nx+1:k*nx);

            delta{i} = D1;
            sigma{k} = Dn;
        else
            Di = A((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx);
            Ai = -A((i-1)*nx+1:i*nx,(i-2)*nx+1:(i-1)*nx);

            Dk = A((k-1)*nx+1:k*nx,(k-1)*nx+1:k*nx);
            Akplus1 = -A(k*nx+1:(k+1)*nx,(k-1)*nx+1:k*nx);
            
            delta{i} = Di - Ai / delta{i-1} * Ai';
            sigma{k} = Dk - Akplus1' / sigma{k+1} * Akplus1;
        end
    end

    for i = 1:n
        % calculate Ui.
        Ui = eye(nx);
        for k = 2:i
            Ak = -A((k-1)*nx+1:k*nx,(k-2)*nx+1:(k-1)*nx);
            Ui = inv(Ak)'*delta{k-1}*Ui;
        end
        
        for j = i:n
            % calculate Vj.
            VjT = inv(sigma{1});
            for k = 2:j
                Ak = -A((k-1)*nx+1:k*nx,(k-2)*nx+1:(k-1)*nx);
                VjT = VjT*Ak'/sigma{k};
            end

            Ainv((i-1)*nx+1:i*nx,(j-1)*nx+1:j*nx) = Ui*VjT;
            if i~=j
                Ainv((j-1)*nx+1:j*nx,(i-1)*nx+1:i*nx) = (Ui*VjT)';
            end

        end

    end

end