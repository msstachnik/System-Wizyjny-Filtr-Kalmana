function [xe,P2,K]=KalmanMS(A,B,C,R,Q,P1,u,x,y)
% function to calculate estymation of state using kalman filter.
% Inputs:
% A - A (state) matrix of state
% B - B (control) matrix of state
% C - C (output) matrix of state
% R - matrix of covariance of measurement noise
% Q - matrix of covariance of proces noise
% P1 - covariance of last state
% u - vector of control in time t0
% x - vector of state in time t0
% y - vector of output in time t0 + dt
% Outputs:
% xe - estimated state
% P2 - covariance of estimated state
% K - Kalman Gain
%  Last Update 2015-06-08 21:59 Mateusz Stachnik

% covariance a priori
P1_ = A * P1 * A' + Q; 
% kalman gain
K = P1_ * C' / (C * P1_ * C' + R);
% state a priori
x_ = A * x + B * u;
% estimated state
xe = x_ + K * (y - C * x_);
% covariance of estimated state
P2 = (eye(size(A)) - K * C) * P1_;
end

