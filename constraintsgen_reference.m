function [ineqconstraint, eqconstraints]=constraintsgen_reference(LTI,dim,dtilde,yref,con)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C zeros(dim.ny,dim.nu)];
eqconstraints.b=[zeros(dim.nx,1); yref];

ineqconstraint.A = [eye(dim.nx + dim.nu); -eye(dim.nx + dim.nu)];
ineqconstraint.b = [con.xmax; con.umax; -con.xmin; -con.umin];

end