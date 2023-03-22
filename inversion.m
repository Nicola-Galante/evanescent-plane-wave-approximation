%% Description:
%
% Cumulative density function ('cumdensity') inversion through the bise-
% ction method in order to find the values of the evanescence parameter
% 'zeta', according to the Inversion Transform Sampling (ITS) technique.
% The initial samples in [0,1] are stored in 'target'.

%% inversion
%
%  INPUT:
%
% - k:      wavenumber
% - L:      truncation parameter
% - target: target values array
% - tol:    bisection method tolerance
%
%  OUTPUT:
%
% - zeta:   evanescence parameter array

function zeta = inversion(k,L,target,tol)

S=10; lt=length(target);

while cumdensity(k,L,S)<1-1e-03; S=S+10; end
zeta=bisection(@(z)cumdensity(k,L,z),zeros(1,lt),S*ones(1,lt),target,tol);

end