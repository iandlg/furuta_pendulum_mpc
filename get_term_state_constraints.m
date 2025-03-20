function [Ff, ff] = get_term_state_constraints(LTI, con ,dim)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   Terminal State Constraint : assume the same as state to be ok
LTI.K = -LTI.K;
LTI.Aclosed = LTI.A + LTI.B*LTI.K;

LTI.furutaClosed = LTISystem('A', LTI.Aclosed);
eig(LTI.Aclosed)

LTI.furutaClosed.x.min = con.xmin; 
LTI.furutaClosed.x.max = con.xmax;

con.terminal_state_set = LTI.furutaClosed.invariantSet();
disp(['Terminal state done. Empty : ', num2str(con.terminal_state_set.isEmptySet())])
% tol1 = 0.2; tol2 = 0.2; tol3 = 1;
% con.Ff = con.F;
% con.ff = [tol1; tol1; tol2; tol2; tol3; tol3; tol3; tol3; tol3; tol3];

Ff = con.terminal_state_set.A; ff = con.terminal_state_set.b;
con.F = [eye(dim.nx); -eye(dim.nx)];
con.f = [con.xmax;  -con.xmin];
con.state_set = Polyhedron(con.F,con.f);
if con.state_set.contains(con.terminal_state_set)
    disp('P is fully contained in Q.');
else
    disp('P is NOT fully contained in Q.');
end
end
