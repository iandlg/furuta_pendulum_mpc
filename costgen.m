function [H,h] = costgen(cost, LTI,xk)
%costgen Generate quadprog matrices H and h
% J = U'HU + 2hU
H = LTI.Bl'*cost.Ql*LTI.Bl + cost.Rl;
h = LTI.Bl'*LTI.Al*xk;
end
