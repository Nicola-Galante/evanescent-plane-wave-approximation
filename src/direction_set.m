%% Description:
%
% Function defining the direction set given the parameters 'y=(th1,th2,th3,
% zeta)' and the parameter sampling strategy 'smpl'. Fixed an upward
% direction 'dz', dependent on 'zeta', the direction matrix 'd' is built
% through the rotations of 'dz' associated to the Euler angles 'th=(th1,
% th2,th3)'. If smpl='deterministic', the construction procedure varies
% due to the Cartesian product structure of the parameters.

%% direction_set
%
%  INPUT:
%
% - k:      wavenumber
% - th1:    Euler angle 1 array
% - th2:    Euler angle 2 array
% - th3:    Euler angle 3 array
% - zeta:   evanescence parameter array
% - smpl:   parameter sampling strategy
%
%  OUTPUT:
%
% - d:      direction matrix

function d = direction_set(k,th1,th2,th3,zeta,smpl)

% Definition of the rotation matrix associated with the y-axis and the
% z-axis. 

N=length(zeta);
RY=@(s)[cos(s),0,sin(s);0,1,0;-sin(s),0,cos(s)];
RZ=@(s)[cos(s),-sin(s),0;sin(s),cos(s),0;0,0,1];

% Costruction of the direction matrix if smpl='deterministic'.

if strcmp(smpl,'deterministic')
    d=[]; dz=[1i*sqrt((zeta/(2*k)+1).^2-1);zeros(1,N);zeta/(2*k)+1];
    for i=1:N
        R1=RZ(th2(i));
        for j=1:N
            R2=RY(th1(j));
            for h=1:N
                R3=RZ(th3(h));
                R=R1*R2*R3; d=[d,R*dz];
            end
        end
    end

% Costruction of the direction matrix in any other case.

else
    d=zeros(3,N);
    for i=1:N
        dz=[1i*sqrt((zeta(i)/(2*k)+1).^2-1);0;zeta(i)/(2*k)+1];
        R1=RZ(th2(i)); R2=RY(th1(i)); R3=RZ(th3(i));
        R=R1*R2*R3; d(:,i)=R*dz;
    end
end

end