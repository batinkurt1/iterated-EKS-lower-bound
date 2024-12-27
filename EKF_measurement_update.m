function [estimated_state,estimated_covariance] = EKF_measurement_update(predicted_state,predicted_covariance,predicted_measurement,actual_measurement,h,R)
    H = d_dx(h,predicted_state);
    S = H*predicted_covariance*H' + R;
    K = predicted_covariance*H'/S;
    estimated_state = predicted_state + K*(actual_measurement - predicted_measurement);
    estimated_covariance = predicted_covariance - K*S*K';
    estimated_covariance = 1/2*(estimated_covariance+estimated_covariance');
end