%% Description:
%
% Approximation of the target function 'u' through plane wave expansion.
% The approximated function 'U' is reconstructed by sampling the target
% function 'u' on the boundary of the unit ball (Dirichlet sampling). The
% evanescent plane waves defining the approximation set are determined
% through a sampling strategy 'smpl' of the parameters 'y=(th1,th2,th3,
% zeta)'. Errors and stability estimates are also provided.

%% Dirichlet_sampling
%
%  INPUT:
%
% - k:      wavenumber
% - u:      target function 
% - L:      truncation parameter
% - P:      approximation set dimension
% - e:      regularization parameter
% - smpl:   parameter sampling strategy
%
%  OUTPUT:
%
% - U:      approximated function
% - err:    absolute error function
% - acc:    relative residual
% - stab:   stability measure

function [U,err,acc,stab] = Dirichlet_sampling(k,u,L,P,e,smpl)

% Selection of the sampling strategy and computation of the related
% parameters 'y=(th1,th2,th3,zeta)'. The weights 'W' and the directions 'd'
% of the plane waves associated with these parameters are then computed.
% The approximation set dimension 'P' is possibly resized due to the
% sampling strategy 'smpl' choosen. If smpl='propagative', only propagative
% plane waves are employed.

if strcmp(smpl,'deterministic') || strcmp(smpl,'random') || strcmp(smpl,'sobol')
    if strcmp(smpl,'deterministic')
        Ps=ceil(P^(1/4)); P=Ps^4;
        th1=linspace(0,1,Ps+1); th1=th1(1:Ps); th1=th1+th1(2)/2;
        th2=linspace(0,2*pi,Ps+1); th2=th2(1:Ps);
        th3=linspace(0,2*pi,Ps+1); th3=th3(1:Ps);        
        zeta=linspace(0,1,Ps+1); zeta=zeta(1:Ps); zeta=zeta+zeta(2)/2;
    elseif strcmp(smpl,'random') || strcmp(smpl,'sobol')
        if strcmp(smpl,'random'); points=rand(4,P);
        else; points=net(scramble(sobolset(4),'MatousekAffineOwen'),P)'; end
        th1=points(1,:); th2=2*pi*points(2,:); th3=2*pi*points(3,:);
        zeta=points(4,:);
    end
    th1=acos(1-2*th1); zeta=inversion(k,L,zeta,1e-12);
    W=sqrt(1./(P*TruncKernel(k,L,zeta,smpl)));
    d=direction_set(k,th1,th2,th3,zeta,smpl);
elseif strcmp(smpl,'extremal_random') || strcmp(smpl,'extremal_sobol')
    Ps=ceil(sqrt(P)); P=Ps^2; d=MD(Ps); d=d(:,1:3)';
    [th2,th1,~]=cart2sph(d(1,:),d(2,:),d(3,:)); th1=pi/2-th1;
    if strcmp(smpl,'extremal_random'); points=rand(2,P);
    else; points=net(scramble(sobolset(2),'MatousekAffineOwen'),P)'; end
    th3=2*pi*points(1,:); zeta=points(2,:); zeta=inversion(k,L,zeta,1e-12);
    W=sqrt(1./(P*TruncKernel(k,L,zeta,smpl)));
    d=direction_set(k,th1,th2,th3,zeta,smpl);
elseif strcmp(smpl,'propagative')
    Ps=ceil(sqrt(P)); P=Ps^2; d=MD(Ps); W=sqrt(1/P); d=d(:,1:3)'; 
else; warning('Error in sampling method');
end

% Sampling nodes on the boundary 'X' and related weights 'WX' imported
% from the 'MD folder' thanks to the 'MD function'. The number of sampling
% points on the boundary of the unit ball is given by 'S^2', the smaller
% integer square greater than or equal to '2*P', and therefore the
% (overdetermined) linear system 'A*xi=b' that will be solved is of size
% 'S^2' by 'P'.

S=ceil(sqrt(2*P)); X=MD(S); WX=X(:,4); X=X(:,1:3);

% Computation of the matrix 'A' and the right-hand-side 'b'.

A=sqrt(WX).*approx_set(k,d,W,X); b=sqrt(WX).*u(X);

% Solution of the linear system 'xi' through regularization tuned by the
% input parameter 'e' together with its associated expansion 'U'.

xi=solve_RSVD(A,b,e); U=@(x)approx_set(k,d,W,x)*xi;

% Errors and stability estimate.

err=@(x)abs(U(x)-u(x)); acc=norm(A*xi-b)/norm(b); stab=norm(xi);

end