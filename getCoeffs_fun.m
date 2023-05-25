function normcoeff = getCoeffs_fun(n,fun,eigenf)

% given an initial function f obtain the first n coefficients of its
% expansion in terms of the eigenstates (stationary states)

% initialize coeff vector
coeff = zeros(n,1);

% loop through and integrate: c_n = ∫ conj(eigenf)*fun dx
for i = 1:n
    coeff(i) = integral(@(theta) conj(eigenf(i,theta)).*fun(theta),0,2*pi);
end

% normalize coefficients such that sum |c_n|^2 = 1 as required by
% probability
normcoeff = normalize(coeff,'norm',2);
end