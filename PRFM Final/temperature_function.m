function func_T = temperature_function(h,rho,T,Vx,Cp,H_i,dY_i_dx,mL,m,hL,QL,peri,T_prev,drho_dx)

R = 8.31;

summation_term = sum(H_i.*dY_i_dx);
% the_term=(mL*hL-QL)*peri/m
% drhoo=Vx^2*drho_dx/rho
% mlll=-Vx^2*mL*peri/m
dT_dx = Vx^2*(drho_dx/rho - mL*peri/m)/Cp -(summation_term - (mL*hL-mL*(Cp*T+0.5*Vx^2)-QL)*peri/m)/Cp;
func_T = T - T_prev -h*dT_dx;

