%% Description:
%
% Plot of the function 'u' both on inner sections of the unit ball (with
% resolution 'n1') and on the unit sphere (with resolution 'n2'). If you
% are interested in a logarithmic scale, the optional input parameter
% 'Clim' is a 1 by 2 array adjusting the minimum colorbar values for the
% two displayed plots. Otherwise, it can be left blank.

%% plot_sph
%
%  INPUT:
%
% - u:      function
% - n1:     resolution parameter 1
% - n2:     resolution parameter 2
% - Clim:   colorbar minimum value array

function plot_sph(u,n1,n2,Clim)

% If you are not interested in a logarithmic scale, 'Clim' is set blank.

if nargin==3; Clim=[]; end

% Computation of the function values 'U' on inner sections of the unit ball.

if ~mod(n1,2); n1=n1+1; end; n=linspace(-1,1,n1); [X,Y,Z]=meshgrid(n,n,n);
XX=[repmat(n,1,n1);repelem(n,1,n1)]; XX=[repmat(n,1,n1^2);repelem(XX,1,n1)]';
IND=~all(XX,2); UX=u(XX(IND,:)); U=zeros(length(XX),1); U(IND)=UX;
U=reshape(U,n1,n1,n1); U(X.^2+Y.^2+Z.^2>1)=nan;

% Plot of the function 'u' on inner sections of the unit ball.

figure
slice(X,Y,Z,real(U),0,0,0); view([-38,30]);
axis off; axis equal; shading interp;
TT=turbo; colormap(TT(1:128,:)); colorbar;
if ~isempty(Clim)
    clim([Clim(1),max(max(max(real(U))))]);
    set(gca,'ColorScale','log'); colormap(hot);
end

% Computation of the function values 'U' on the unit sphere.

[X,Y,Z]=sphere(n2); n=n2+1; X_=reshape(X,n^2,1); Y_=reshape(Y,n^2,1);
Z_=reshape(Z,n^2,1); U=u([X_,Y_,Z_]); U=reshape(U,n,n);

% Plot of the function 'u' on the unit sphere.

figure
surf(X,Y,Z,real(U)); view([-128,30]);
axis off; axis equal; shading interp;
TT=turbo; colormap(TT(1:128,:)); colorbar;
if ~isempty(Clim)
    clim([Clim(2),max(max(real(U)))]);
    set(gca,'ColorScale','log'); colormap(hot);
end

end