function H = d_dx(h,x)
    I = eye(length(x));
    h_x =h(x); 
    H = zeros(length(h_x),length(x));
    myeps = 1e-7;
    for k = 1:length(x)
        xtemp = x + myeps*I(:,k);
        H(:,k) = (h(xtemp) - h_x)/myeps;
    end
end