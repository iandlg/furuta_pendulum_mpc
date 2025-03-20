function [x, posCov] = dynamic_update(x,u,LTI, cov)
% Kalman filter dynamic update
x = LTI.A*x + LTI.B*u;
posCov = LTI.A*cov.pos*LTI.A' + cov.motion;
end
