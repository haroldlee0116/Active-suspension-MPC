function eqconstraints=eqconstraintsgen(LTI,dim,dtilde)

eqconstraints.A=[eye(dim.nx)-LTI.A -LTI.B; LTI.C LTI.D];
eqconstraints.b=[LTI.Bd*dtilde; LTI.yref-LTI.Cd*dtilde];

end
