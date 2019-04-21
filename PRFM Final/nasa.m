function [cp_i,H_i,S_i]=nasa(a,temp,T,R)
% [a,temp]=thermodat(spcs);
q=size(a);
% T=1000;
% R=8.314;
cp_i=zeros(q(1),1);
H_i=zeros(q(1),1);
S_i=zeros(q(1),1);

%%%% For Temperature
% T_low=temp(:,1);
% T_high=temp(:,2);
% T_mean=temp(:,3);

for i=1:q(1)
%%ERROR PRINT
if T<temp(i,1)
    disp('Lower temperature out of range')
elseif T>temp(i,2)
    disp('Upper temperature out of range')
elseif T>=temp(i,3) && T <=temp(i,2)
  cp_i(i,1)=(a(i,1)+a(i,2)*T+a(i,3)*T^2+a(i,4)*T^3+a(i,5)*T^4)*R;
  H_i(i,1)=(a(i,1)*T+a(i,2)*T^2/2+a(i,3)*T^3/3+a(i,4)*T^4/4+a(i,5)*T^5/5+a(i,6))*R;
  S_i(i,1)=(a(i,1)*log(T)+a(i,2)*T+a(i,3)*T^2/2+a(i,4)*T^3/3+a(i,5)*T^4/4+a(i,7))*R;
 elseif T<=temp(i,3) && T>= temp(i,1)
  cp_i(i,1)=(a(i,8)+a(i,9)*T+a(i,10)*T^2+a(i,11)*T^3+a(i,12)*T^4)*R;
  H_i(i,1)=(a(i,8)*T+a(i,9)*T^2/2+a(i,10)*T^3/3+a(i,11)*T^4/4+a(i,12)*T^5/5+a(i,13))*R;
  S_i(i,1)=(a(i,8)*log(T)+a(i,9)*T+a(i,10)*T^2/2+a(i,11)*T^3/3+a(i,12)*T^4/4+a(i,14))*R;
end
end
end
