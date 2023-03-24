%% Description:
%
% Implementation of a simple bisection method for non-decreasing functions,
% such as the cumulative density function 'cumdensity'.

%% bisection
%
%  INPUT:
%
% - f:      function
% - a:      left endpoint
% - b:      right endpoint
% - fx:     input value
% - tol:    tolerance parameter
%
%  OUTPUT:
%
% - x:      target    

function x = bisection(f,a,b,fx,tol)

f=@(x)f(x)-fx;

while true
    x=(a+b)/2;
    if ~any(abs(b-x)>tol); break; end
    index=(sign(f(x))==ones(1,length(fx)));
    a(~index)=x(~index); b(index)=x(index);  
end

end