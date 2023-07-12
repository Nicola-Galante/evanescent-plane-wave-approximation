%% Description:
%
% This function takes a degree 'l' and an array of evaluation points 'x'
% as inputs and returns a tensor 'T'. For each element 'x_i' in 'x', the
% matrix 'T(:,:,i)' contains all the Associated Legendre polynomials of
% every non-negative order up to degree 'l', evaluated at 'x_i'. Since we
% assume that the Associated Legendre polynomial 'P_l^m = 0 if m > l', the
% matrix 'T(:,:,i)' is upper triangular, with the order varying along the
% first component and the degree varying along the second component.

%% AssociatedLegendre
%
%  INPUT:
%
% - l:      degree
% - x:      evaluation points array
%
%  OUTPUT:
%
% - T:      Associated Legendre polynomials tensor

function T = AssociatedLegendre(l,x)

lx=length(x); if l==0; T=ones(1,1,lx); return; end
T=zeros(l+1,l+1,lx); [Q,V]=LegendreDiag2(l,x);

for m=0:l-2
    W = LegendreLine(l,m,Q(:,m+1),V(:,m+1),x);
    T(m+1,m+1:l+1,:)=reshape(W',1,[],lx);
end

T(l,l,:)=Q(:,l); T(l+1,l+1,:)=Q(:,l+1); T(l,l+1,:)=V(:,l);

end

%% Description:
%
% This function exploit the recursive relation '(l-m)*P_l^m(x) = (2*l-1)*x*
% P_{l-1}^m(x)-(l+m-1)*P_{l-2}^m(x)' in order to compute all the Associated 
% Legendre Polynomials 'P_l^m(x)' for a non-negative fixed order 'm', where
% the initial values are provided by the previous functions
% 'LegendreDiag1', i.e. 'P_l^l(x)', and 'LegendreDiag2', i.e.
% 'P_l^{l-1}(x)'. More precisely, the input parameters 'I' and 'J' are
% arrays containing the initial values, as the evaluation points in 'x'
% vary. The matrix 'Q' first component represents the values of 'x' at
% which the Associated Legendre polynomials 'P_l^m' are evaluated, while
% the second component represents the degree of the polynomials.

%% LegendreLine
%
%  INPUT:
%
% - l:      degree
% - m:      order
% - I:      first initial values array
% - J:      second initial values array
% - x:      evaluation points array
%
%  OUTPUT:
%
% - Q:      Associated Legendre polynomials matrix

function Q = LegendreLine(l,m,I,J,x)

lx=length(x);

if l==m; Q=ones(lx,1).*I;
elseif l==m+1; Q=ones(lx,2).*[I,J];
else
    A=LegendreLine(l-1,m,I,J,x);
    Q=[A(:,1:l-m),((2*l-1)*x.*A(:,l-m)-(l+m-1)*A(:,l-m-1))/(l-m)];
end

end

%% Description:
%
% This function exploit the recursive relation 'P_l^l(x) = (2*l-1)*sqrt(x^2
% -1)*P_{l-1}^{l-1}(x)' in order to compute all the 'first-diagonal'
% Associated Legendre Polynomials 'P_l^l(x)' starting from the initial
% values 'P_0(x) = 1'. The matrix 'Q' first component represents the values
% of 'x' at which the Associated Legendre polynomials 'P_l^l' are evaluated,
% while the second component represents the degree of the polynomials.

%% LegendreDiag1
%
%  INPUT:
%
% - l:      degree
% - x:      evaluation points array
%
%  OUTPUT:
%
% - Q:      Associated Legendre polynomials matrix

function Q = LegendreDiag1(l,x)

if l==0; Q=ones(length(x),1);
else
    A=LegendreDiag1(l-1,x);
    Q=[A(:,1:l),(2*l-1)*sqrt(x.^2-1).*A(:,l)];
end

end

%% Description:
%
% This function exploit the recursive relation 'P_l^{l-1}(x) = x*(2*l-1)*
% P_{l-1}^{l-1}(x)' in order to compute all the 'second-diagonal'
% Associated Legendre Polynomials 'P_l^{l-1}(x)', where the initial values
% are provided by the previous function 'LegendreDiag1', i.e. 'P_l^l(x)'.
% The matrix 'V' first component represents the values of 'x' at which the
% Associated Legendre polynomials 'P_l^{l-1}' are evaluated, while the
% second component represents the degree of the polynomials. 'Q' is the
% matrix inherited from the 'LegendreDiag1' function.

%% LegendreDiag2
%
%  INPUT:
%
% - l:      degree
% - x:      evaluation points array
%
%  OUTPUT:
%
% - Q:      Associated Legendre polynomials matrix
% - V:      Associated Legendre polynomials matrix

function [Q,V] = LegendreDiag2(l,x)

Q=LegendreDiag1(l,x); V=x*(2*(1:l)-1).*Q(:,1:l);

end
