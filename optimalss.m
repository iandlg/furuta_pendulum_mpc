function[xref,uref]=optimalss(LTI,dim,weight,constraints,eqconstraints)

H=blkdiag(0,1,0,0,0,eye(dim.nu));
h=zeros(dim.nx+dim.nu,1);


options1 = optimoptions(@quadprog); 
options1.OptimalityTolerance=1e-20;
options1.ConstraintTolerance=1.0000e-15;
options1.Display='off';
xur=quadprog(H,h,constraints.A,constraints.b,eqconstraints.A,eqconstraints.b,[],[],[]);
xref=xur(1:dim.nx);
uref=xur(dim.nx+1:end);

end