% Admissible set using MPT3 ( not tractable with 5 states )
LTI.furutaMPT3 = LTISystem('A', LTI.A, 'B', LTI.B ,'Ts', Ts);

LTI.furutaMPT3.x.min = [-tol1;-tol2;-tol3; -tol3; -tol3]; 
LTI.furutaMPT3.x.max = [ tol1; tol2; tol3; tol3; tol3];

LTI.furutaMPT3.u.min = -10;
LTI.furutaMPT3.u.max = 10;

con.terminal_state_set = LTI.furutaMPT3.invariantSet();
