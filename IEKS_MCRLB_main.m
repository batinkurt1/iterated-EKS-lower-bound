clearvars; clc; 
close all;

plot_open = 0;

Nmc = 30000;
apply_numerical_derivative = 0; % 0 for analytical, 1 for numerical d_dx
apply_tridiagonal_method   = 0; % 0 for built-in inverse, 1 for tridiagonal method
 
% load the true state sequence
load("trueTarget.mat");

t = trueTarget(1, :);
simulation_length = length(t);

% assuming a fixed sampling period
T = t(2) - t(1);

x_true = trueTarget(2, :);
y_true = trueTarget(3, :);

% generate the velocity information
true_states = [x_true;y_true;[(x_true(2)-x_true(1))/T,diff(x_true)/T];[(y_true(2)-y_true(1))/T,diff(y_true)/T]];
state_dimension = size(true_states,1);

if plot_open
    figure;
    plot(x_true,y_true,LineWidth=1.5,Color="#77AC30");
    title("True Target Trajectory");
    ylabel("y position");
    ylim([0,2500]);
    xlabel("x position");
    xlim([500,3000]);
    legend("True Target Trajectory");
    grid on;
end


A = [eye(2),T*eye(2);zeros(2),eye(2)];
B = eye(4);
Q = 10*diag([10,10,1,1]);

f = @(x)(A*x);

f1 = @(x)(A(1,:)*x);
f2 = @(x)(A(2,:)*x);
f3 = @(x)(A(3,:)*x);
f4 = @(x)(A(4,:)*x);

h = @(x)[sqrt(x(1,:).^2+x(2,:).^2);atan2(x(2,:),x(1,:))];

h1 = @(x)(sqrt(x(1,:).^2+x(2,:).^2));
h2 = @(x)(atan2(x(2,:),x(1,:)));

% analytical derivative relations

df_dx = @(x)A;

d2f_dx2 = @(x)zeros(4);

dh_dx = @(x)[x(1)/sqrt(x(1).^2+x(2).^2),x(2)/sqrt(x(1).^2+x(2).^2),0,0;...
            -x(2)/(x(1).^2+x(2).^2),x(1)/(x(1).^2+x(2).^2),0,0];

d2h1_dx2 = @(x)[x(2)^2/sqrt(x(1).^2+x(2).^2)^3,-x(1)*x(2)/sqrt(x(1).^2+x(2).^2)^3,0,0;...
               -x(1)*x(2)/sqrt(x(1).^2+x(2).^2)^3,x(1)^2/sqrt(x(1).^2+x(2).^2)^3,0,0;...
               0,0,0,0;...
               0,0,0,0];
d2h2_dx2 = @(x)[2*x(1)*x(2)/(x(1).^2+x(2).^2)^2,(-x(1)^2+x(2)^2)/(x(1).^2+x(2).^2)^2,0,0;...
               (-x(1)^2+x(2)^2)/(x(1).^2+x(2).^2)^2,-2*x(1)*x(2)/(x(1).^2+x(2).^2)^2,0,0;...
               0,0,0,0;...
               0,0,0,0];

sigma_r = 50;
sigma_theta = 0.5*pi/180;
R = diag([sigma_r^2,sigma_theta^2]);

hbar = @(x)[sqrt(x(1,:).^2+x(2,:).^2);atan2(x(2,:),x(1,:))];
sigma_rbar = 100;
sigma_thetabar = 1.5*pi/180;
Rbar = diag([sigma_rbar^2,sigma_thetabar^2]);

true_measurements_in_polar = hbar([x_true;y_true]);
r_true = true_measurements_in_polar(1,:);
azimuth_true = true_measurements_in_polar(2,:);

x0_mean = [1000; 1000; 0; 0];
P0 = diag([100^2, 100^2, 10^2, 10^2]);
x0 = mvnrnd(x0_mean,P0)';

% determine the pseudotrue state sequence

[pseudotrue_estimated_states,pseudotrue_estimated_covariances] = IEKS(true_measurements_in_polar,A,B,Q,h,R,simulation_length,x0,P0,0,apply_numerical_derivative,dh_dx);

MCRLB = calculate_MCRLB(pseudotrue_estimated_states,f,f1,f2,f3,f4,h,h1,h2,Q,R,Rbar,true_measurements_in_polar,simulation_length,state_dimension,x0,P0,df_dx,dh_dx,d2f_dx2,d2h1_dx2,d2h2_dx2,apply_numerical_derivative,apply_tridiagonal_method);

LB = cell(1,simulation_length);

for k = 1:simulation_length
    row_start = (k-1) * state_dimension + 1;
    row_end = k * state_dimension;
    col_start = row_start;
    col_end = row_end;
    LB{k} = MCRLB(row_start:row_end,col_start:col_end) + (true_states(:,k)-pseudotrue_estimated_states(:,k))*(true_states(:,k)-pseudotrue_estimated_states(:,k))';
end

% perform monte carlo simulation and compare
MSE = cell(1,simulation_length);

% put dummy values in each cell to prevent type error
for k = (t+1)
    MSE{k} = zeros(4);
end

for i = 1:Nmc
    % generate noisy measurements
    noisy_measurements = zeros(2,simulation_length);
    for k = (t+1)
        noisy_measurements(:,k) = generate_measurements(true_measurements_in_polar(:,k),[0;0],Rbar);
    end
    [estimated_states,estimated_covariances] = IEKS(noisy_measurements,A,B,Q,h,R,simulation_length,x0,P0,i,apply_numerical_derivative,dh_dx);
    % (t+1) = 1 to 151, indexing in matlab starts from 1
    for k = (t+1)
        MSE{k} = MSE{k} + (estimated_states(:,k) - true_states(:,k))*(estimated_states(:,k) - true_states(:,k))';
    end
    
end

% normalize with Nmc to compute the average MSE at each time step
for k = (t+1)
    MSE{k} = 1/Nmc*MSE{k};
end

LB_x = zeros(1,simulation_length);
MSE_x = zeros(1,simulation_length);

for k= 1:simulation_length
    current_LB = LB{k};
    current_MSE = MSE{k};

    LB_x(k) = sqrt(current_LB(1,1));
    MSE_x(k) = sqrt(current_MSE(1,1));
end

if plot_open
    figure;
    plot(x_true,y_true,LineWidth=1.5,Color="#77AC30");
    hold on;
    title("True Target Trajectory");
    ylabel("y position");
    ylim([0,2500]);
    xlabel("x position");
    xlim([500,3000]);
    grid on;
    
    plot(estimated_states(1,:),estimated_states(2,:),LineWidth=1.5);
    plot(pseudotrue_estimated_states(1,:),pseudotrue_estimated_states(2,:),LineWidth=1.5);
    scatter(noisy_measurements(1,:).*cos(noisy_measurements(2,:)),noisy_measurements(1,:).*sin(noisy_measurements(2,:)));
    scatter(true_measurements_in_polar(1,:).*cos(true_measurements_in_polar(2,:)),true_measurements_in_polar(1,:).*sin(true_measurements_in_polar(2,:)));
    legend("True Target Trajectory","misspecified state sequence","pseudotrue state sequence","noisy measurements","true measurements");
end

% final result, the lower bound and the MSE of x
figure;
plot(LB_x);
hold on;
plot(MSE_x);
legend('LB_x','MSE_x');

