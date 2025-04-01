function eqconstraints=eqconstraintsgen(LTI,dim,dtilde,yref)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C zeros(dim.ny,dim.nu)];
eqconstraints.b=[LTI.Bdist*dtilde; yref-LTI.Cdist*dtilde];

end
