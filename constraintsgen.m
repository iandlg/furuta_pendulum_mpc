function [ineqconstraint, eqconstraints] = constraintsgen(LTI,dim,dtilde,yref,con)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C zeros(dim.ny,dim.nu)];
eqconstraints.b=[LTI.Bdist*dtilde; yref-LTI.Cdist*dtilde];

ineqconstraint.A = [eye(dim.nx + dim.nu); -eye(dim.nx + dim.nu)];
ineqconstraint.b = [con.xmax; con.umax; -con.xmin; -con.umin];
end