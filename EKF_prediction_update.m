function [predicted_state,predicted_covariance,predicted_measurement] = EKF_prediction_update(prev_state,prev_covariance,A,B,h,Q)
    predicted_state= A*prev_state;
    predicted_covariance = A*prev_covariance*A'+B*Q*B';
    predicted_measurement =h(predicted_state);
end