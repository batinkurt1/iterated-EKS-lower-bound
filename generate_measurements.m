function noisy_measurements = generate_measurements(trueTarget,mu,R)
    N=size(trueTarget,2);
    v_k = mvnrnd(mu,R,N)';
    noisy_measurements = trueTarget + v_k;
end