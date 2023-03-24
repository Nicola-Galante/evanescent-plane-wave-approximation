%% Description:
%
% Main script to compute approximation of Helmholtz solution by plane waves
% within the unit ball.

% 'Dirichlet_sampling' function parameters. The possible input sampling
% strategies are: 'propagative', 'deterministic', 'random', 'sobol',
% 'extremal_random', and 'extremal_sobol'.

k=5;                    % wavenumber
L=5*k;                  % truncation parameter
P=4*(L+1)^2;            % approximation set dimension
e=1e-14;                % regularization parameter
smpl='extremal_sobol';  % parameter sampling strategy

% 'plot_sph' function parameters.

n1=200;                 % resolution parameter 1
n2=500;                 % resolution parameter 2
Clim=[1e-12,1e-9];      % minimum colorbar values

% Construction of the solution surrogate.

M=[];
for l=0:L
    M=[M;l*ones(2*l+1,1),(-l:l)',randn(2*l+1,1)/max(1,l-k)];
end
u=@(x)sol_surrogate(M,k,x);

% Definition of the fundamental solution with singularity in 's'.

%s=[0,-1,1]; u=@(x)exp(1i*k*vecnorm(x-s,2,2))./(4*pi*vecnorm(x-s,2,2));

% Computation of the approximated function 'U'.

[U,err,acc,stab]=Dirichlet_sampling(k,u,L,P,e,smpl);

% Plot of both the target function 'u' and the absolute error 'err'.

disp('plotting target function'); plot_sph(u,n1,n2);
disp('plotting absolute error'); plot_sph(err,n1,n2,Clim);