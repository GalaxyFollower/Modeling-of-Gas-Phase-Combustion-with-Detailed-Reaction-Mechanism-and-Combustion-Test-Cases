%% Convert molar concentration in to mole fraction.
function X_it=molfrctn2(X_iconc)
% in_sumfn = diag(1./sum(X_iconc,2));
% X_it=in_sumfn*X_iconc;
[r,c]=size(X_iconc);
X_it=zeros(r,c);
for i= 1: r
X_it(i,:)=X_iconc(i,:)./sum(X_iconc(i,:),2);
end

