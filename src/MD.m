%% Description:
%
% The extremal point system of degree 's-1' and related weights are
% retrieved from the folder 'MD' (Maximum Determinant). For this function
% to work properly, it is necessary to download the extremal point sets
% and the related weights from
% https://web.maths.unsw.edu.au/~rsw/Sphere/MaxDet/MaxDet1.html#ExtSys
% and put them in the 'MD' folder.

%% MD
%
%  INPUT:
%
% - s:      sqrt(no. of extremal points), integer
%
%  OUTPUT:
%
% - S:      Extremal points and related weights

function S = MD(s)

% Given an integer 's', the extremal points and weights of degree 's-1' are
% retrived from the 'MD' folder by reconstructing the name of the related
% file.

s_str=num2str(s-1); s__str=s_str; s_digit=numel(s_str);
S_str=num2str(s^2); S__str=S_str; S_digit=numel(S_str);

if s_digit==1; s__str=strcat('00',s_str);
elseif s_digit==2; s__str=strcat('0',s_str);
end

if S_digit==1; S__str=strcat('0000',S_str);
elseif S_digit==2; S__str=strcat('000',S_str);
elseif S_digit==3; S__str=strcat('00',S_str);
elseif S_digit==4; S__str=strcat('0',S_str);
end

try
    S=load(sprintf('MD/md%s.%s',s__str,S__str));
catch
    disp(['File not found. Download the extremal point systems and the related weights from' ...
        ' https://web.maths.unsw.edu.au/~rsw/Sphere/MaxDet/MaxDet1.html#ExtSys' ...
        ' and put them in the MD folder']);
end