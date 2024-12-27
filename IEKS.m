function [smoothed_states,smoothed_covariances] = IEKS(measurements,A,B,Q,h,R,simulation_length,x0,P0,MC_number,apply_numerical_derivative,dh_dx)

    estimated_states = zeros(4,simulation_length);
    estimated_covariances = cell(1,simulation_length);
    
    smoothed_states = zeros(4,simulation_length);
    smoothed_covariances = cell(1,simulation_length);
    
    state_diff = inf;
    num_of_iterations = 0;
    
    % first perform the standard EKF up to k = N.
    
    for k=1:simulation_length
        if k == 1
            predicted_state = x0;
            predicted_covariance = P0;
            predicted_measurement = h(x0);
        else
            [predicted_state,predicted_covariance,predicted_measurement] = EKF_prediction_update(estimated_states(:,k-1),estimated_covariances{k-1},A,B,h,Q);
        end
        actual_measurement = measurements(:,k);
        [estimated_state,estimated_covariance] = EKF_measurement_update(predicted_state,predicted_covariance,predicted_measurement,actual_measurement,h,R,apply_numerical_derivative,dh_dx);
        estimated_states(:,k) = estimated_state;
        estimated_covariances{k} = estimated_covariance;
    end
    
    while state_diff > 1e-12
        previous_states = estimated_states;
        
        % N given N is known
        smoothed_states(:,simulation_length) = estimated_states(:,simulation_length);
        smoothed_covariances{simulation_length} = estimated_covariances{simulation_length};
    
        % smoothing
        for k = fliplr(1:simulation_length - 1)
            smoother_gain = estimated_covariances{k} * A' / (A * estimated_covariances{k} * A' + B * Q * B');
            smoothed_states(:,k) = estimated_states(:,k) + smoother_gain * (smoothed_states(:,k+1) - A * estimated_states(:,k));
            smoothed_covariances{k} = estimated_covariances{k} + smoother_gain * (smoothed_covariances{k+1} - A * estimated_covariances{k} * A' - B * Q * B') * smoother_gain';
        end
    
        % filtering
        estimated_states = zeros(4,simulation_length);
        estimated_covariances = cell(1,simulation_length);
        for k=1:simulation_length
            if k == 1
                predicted_state = x0;
                predicted_covariance = P0;
                predicted_measurement = h(x0);
            else
                [predicted_state,predicted_covariance,predicted_measurement] = EKF_prediction_update(estimated_states(:,k-1),estimated_covariances{k-1},A,B,h,Q);
            end
            actual_measurement = measurements(:,k);
            % the jacobian location is taken from smoothed states.
            [estimated_state,estimated_covariance] = IEKF_measurement_update(predicted_state,predicted_covariance,predicted_measurement,actual_measurement,h,R,smoothed_states(:,k),apply_numerical_derivative,dh_dx);
            estimated_states(:,k) = estimated_state;
            estimated_covariances{k} = estimated_covariance;
        end
        
        % difference between two consecutive iterations
        state_diff = max(max(abs(previous_states - estimated_states)));
        num_of_iterations = num_of_iterations + 1;
    end
    
    % perform smoothing one last time
    for k = fliplr(1:simulation_length - 1)
        smoother_gain = estimated_covariances{k} * A' / (A * estimated_covariances{k} * A' + B * Q * B');
        smoothed_states(:,k) = estimated_states(:,k) + smoother_gain * (smoothed_states(:,k+1) - A * estimated_states(:,k));
        smoothed_covariances{k} = estimated_covariances{k} + smoother_gain * (smoothed_covariances{k+1} - A * estimated_covariances{k} * A' - B * Q * B') * smoother_gain';
    end

    fprintf("Number of iterations in MC Number %0.5g : %0.5g \n",MC_number,num_of_iterations);
        

end