%% Description:
%
% Computation of the N-term Christoffel function inverse, i.e.
% '\mu_N^{-1}(zeta)'. It is necessary to treat separately the case in
% which smpl='deterministic', due to the Cartesian product structure of
% the evanescent plane wave parameters.

%% TruncKernel
%
%  INPUT:
%
% - k:      wavenumber
% - L:      truncation parameter
% - zeta:   evanescence parameter array
% - smpl:   parameter sampling strategy
%
%  OUTPUT:
%
% - muN:    N-term Christoffel function inverse, i.e. '\mu_N^{-1}(zeta)'

function muN = TruncKernel(k,L,zeta,smpl)

% Computation of the Associated Legendre polynomials tensor.

Dz=AssociatedLegendre(L,zeta'/(2*k)+1);

% The computation of '\mu_N^{-1}(zeta)' exploits the formula '\mu_N^{-1}
% (zeta) = \sum_{l=0}^L \alpha_l^2*(2*l+1)/(4*pi)*(P_l(zeta/(2*k)+1)^2+2*
% \sum_{m=1}^l*((l-m)!/(l+m)!)*P_l^m(zeta/(2*k)+1)^2)', made possible by
% the fact that '\gamma_l^{-m}*P_l^{-m} = \gamma_l^m*P_l^m'.

G=factorial(max((0:L)-(0:L)',0))./factorial((0:L)+(0:L)'); G=[ones(1,L+1);2*ones(L,L+1)].*G;
muN=sum(alphal(k,0:L).*(2*(0:L)+1)/(4*pi).*reshape(sum(G.*Dz.^2,1),L+1,length(zeta))',2)';

% If smpl='deterministic', due to the Cartesian product structure of the
% evanescent plane wave parameters, we need to repeat the values of
% '\mu_N^{-1}(zeta)' for each Euler angle.

if strcmp(smpl,'deterministic'); muN=repmat(muN,1,length(zeta).^3); end

end

%% Description:
%
% Computation of (the approximation of) the Herglotz densities normali-
% zation constant of degree 'l' and of any order 'm' such that '|m|<=l'.

%% alphal
%
%  INPUT:
%
% - k:      wavenumber
% - l:      degree
%
%  OUTPUT:
%
% - alp:    Herglotz densities normalization constant (approximation)

function alp = alphal(k,l)

alp=(k.^(2*l).*factorial(l))./(2*exp(2*k)*sqrt(pi).*gamma(l+0.5).*igamma(2*l+1.5,2*k));

end