%% Description:
%
% Function defining (the approximation of) the cumulative density function,
% which will then be used to map the initial samples in [0,1] to the
% evanescent parameters 'zeta' through the Inversion Transform Sampling
% (ITS) technique. 

%% cumdensity
%
%  INPUT:
%
% - k:      wavenumber
% - L:      truncation parameter
% - zeta:   evanescence parameter array
%
%  OUTPUT:
%
% - val:    cumulative density function (approximation) values array

function val = cumdensity(k,L,zeta)

N=(L+1)^2; Y=zeros(L+1,length(zeta));
Y_=@(l)(2*l+1)*gammainc(2*k+zeta,2*l+1.5,'upper')/gammainc(2*k,2*l+1.5,'upper');
for l=0:L; Y(l+1,:)=Y_(l); end; val=1-sum(Y,1)/N;

end