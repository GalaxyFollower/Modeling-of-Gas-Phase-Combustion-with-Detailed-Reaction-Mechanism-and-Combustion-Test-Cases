function [func_rho,drho_dx] = density_function(h,CS_area,rho,Vx,Cp,MWmix,M_i,T,H_i,dY_i_dx,mL,m,hL,QL,peri,rho_prev)


R = 8.31;


summation_term_density = sum(dY_i_dx.*(H_i - (1./M_i)*(Cp*T*MWmix)));

tempr_source = Vx^2 * mL*peri/m - peri*(mL*hL-mL*(Cp*T+0.5*Vx^2)-QL)/m;
numerator = rho^2*R*(summation_term_density + tempr_source)/(Cp*MWmix)-2*rho*Vx*mL*peri/CS_area;
denominator = T*rho*R/MWmix*(1+Vx^2/(Cp*T))-rho*Vx^2;

drho_dx = numerator/denominator;

func_rho = rho - rho_prev - h*drho_dx;
