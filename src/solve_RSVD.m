%% Description:
%
% Function able to solve linear system 'A*xi=b' using a regularized SVD.
% The threshold 'e' controls the regularization.  Any singular value
% smaller than or equal to 'Sd(1)*e', where 'Sd(1)' is the largest
% singular value, is approximated by zero when solving the least-squares
% problem. The condition number 'cond' is defined as the ratio of the
% largest over the smallest singular value of the matrix 'A', namely
% 'Sd(1)/Sd(end)'.

%% solve_RSVD
%
%  INPUT:
%
% - A:      rectangular matrix
% - b:      right-hand-side
% - e:      regularization parameter
%
%  OUTPUT:
%
% - xi:     linear system solution
% - cond:   condition number

function [xi,cond] = solve_RSVD(A,b,e)

% SVD computation. The regularization process involves trimming the tail
% of relatively small singular values by setting them to zero.

[U,S,V]=svd(A); Sd=diag(S); Se=spdiags((Sd>Sd(1)*e)./Sd,0,S)';
xi=V*(Se*(U'*b));

% Condition number computation

cond=Sd(1)/Sd(end);

end