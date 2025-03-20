function [x,posCov] = measurement_update(x,y,LTI, cov)
% Measurement Update for Kalman filter
ytilde = y - LTI.C*x;
S = LTI.C*cov.pos*LTI.C' + cov.measurement;
K = cov.pos*LTI.C'/S;
x = x + K*ytilde;
posCov = cov.pos - K*LTI.C*cov.pos;
end
