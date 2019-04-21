function [Y_i,dYi_dt]=mass_fraction(Y_i,M_i,w_i,rho_t,dt)
dYi_dt=(M_i.*w_i)./rho_t;
Y_i=Y_i+dYi_dt*dt;
end