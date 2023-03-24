%% Description:
%
% Function that generates a set of plane waves 'set' given their directions
% 'd' and weights 'W'. The plane wave set is evaluated at the points stored
% in 'x'.

%% approx_set
%
%  INPUT:
%
% - k:      wavenumber
% - d:      direction matrix
% - W:      weight array
% - x:      evaluation points matrix
%
%  OUTPUT:
%
% - set:    matrix of evanescent plane waves evaluated at 'x'

function set = approx_set(k,d,W,x)

set=W.*exp(1i*k*(x*d));

end