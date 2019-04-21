function [func_Y_i,dY_i_dx] = species_con_function(h,Y_i,rho,Vx,M_i,w_i,mL,mL_species,m,peri,Y_i_prev)

source_term = mL_species*peri/m - Y_i*mL*peri/m;
dY_i_dx = (w_i.*M_i)/(rho*Vx) + source_term;
% [r noneed]=size(Y_i(:,1));
% for(i=1:r)
%     
% if(Y_i(i,1)<0)
%    Y_i(i,1) =0;
% end
% end
func_Y_i = Y_i - Y_i_prev - h*dY_i_dx;
