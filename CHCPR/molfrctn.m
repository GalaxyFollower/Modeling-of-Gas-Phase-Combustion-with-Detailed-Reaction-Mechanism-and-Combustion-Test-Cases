%% Convert molar concentration in to mole fraction.
function X_it=molfrctn(X_iconc)
X_it=X_iconc./sum(X_iconc);
end
