%% Description:
%
% Computation of a solution surrogate. 'M' is a matrix whose lines are
% of the form '[l,m,coeff_l^m]', where 'l' is the degree, 'm' is the order,
% and 'coeff_l^m' is the associated coefficient. The solution surrogate is
% defined as a superposition of spherical waves 'b_l^m', i.e. '\sum_{l,m}
% coeff_l^m*b_l^m'. The solution surrogate is evaluated at the points
% stored in 'x'.

%% sol_surrogate
%
%  INPUT:
%
% - M:      modal matrix
% - k:      wavenumber
% - x:      evaluation points matrix
%
%  OUTPUT:
%
% - val:    surrogate solution values array

function val = sol_surrogate(M,k,x)

[rx,~]=size(x); [rU,~]=size(M); val=zeros(rx,1);

% Computation of the solution surrogate as a superposition of
% spherical waves 'b_l^m'.

for i=1:rU
    b=bl(k,M(i,1),x').'; val=val+b(:,M(i,1)+M(i,2)+1)*M(i,3);
end

end

%% Description:
%
% Computation of the spherical waves 'b_l' of degree 'l'. The output matrix
% 'val' is structured as follows: the spherical waves order 'm' varies
% along the first component, while the points at which the functions 'b_l'
% are evaluated vary along the second component.

%% bl
%
%  INPUT:
%
% - k:      wavenumber
% - l:      degree
% - x:      evaluation points matrix
%
%  OUTPUT:
%
% - val:    spherical waves matrix

function val = bl(k,l,x)

[th2,th1,r]=cart2sph(x(1,:),x(2,:),x(3,:)); th1=pi/2-th1;

Plm=legendre(l,cos(th1));
gamma=sqrt(((2*l+1)*factorial(l-(0:l)'))./(4*pi*factorial(l+(0:l)')));
Plm=gamma.*Plm; Plm=[flipud((-1).^(1:l)'.*Plm(2:l+1,:));Plm];

val=betal(k,l)*sphbes(l,k*r).*Plm.*exp(1i*(-l:l)'*th2);

end

%% Description:
%
% Computation of the normalization constant of the spherical waves of
% degree 'l' and of any order 'm' such that '|m|<=l'.

%% betal
%
%  INPUT:
%
% - k:      wavenumber
% - l:      degree
%
%  OUTPUT:
%
% - nrm:    normalization coefficient

function nrm = betal(k,l)

L2=pi/(4*k)*(besselj(l+0.5,k).^2-besselj(l-0.5,k).*besselj(l+1.5,k));
H1=2*L2+pi/(2*k^3)*(l.*besselj(l+0.5,k).^2-k*besselj(l+0.5,k).*besselj(l+1.5,k));
nrm=1./sqrt(H1);

end

%% Description:
%
% Computation of the spherical Bessel function of order 'l'. The function
% is evaluated at the points stored in 'x'.

%% sphbes
%
%  INPUT:
%
% - l:      degree
% - x:      evaluation points array
%
%  OUTPUT:
%
% - jl:     spherical Bessel function values array

function jl = sphbes(l,x)

jl=sqrt(pi./(2*x)).*besselj(l+.5,x);

if l==0; jl(isnan(jl))=1;
else; jl(isnan(jl))=0;
end

end