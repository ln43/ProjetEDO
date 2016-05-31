function sol=TravellingWaveKPP2D()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
%


close all;
clear all;

d=0.5;
alpha=0.2;
K=0.3;
u0=0.7;

a = -50;
b = 50;
nel = 200; % number of elements
h = (b-a)/nel; % step size
nv = nel+1;% number of vertices

x = a:h:b; % mesh
y = a:h:b;
dt = 0.01; % time steps
tend = 50; 

e=ones(nv, 1);

% Discretization of the Laplace operator with Neumann bc


% initial guess
u = zeros(201,201);
for i=1:1:201
    for j=1:1:201
        if i>80 & i<120 & j<120 & j>80
            u(i,j)=u0;
        end
    end
end
counter=1;
ind=1;
index1=1;
figure(1)

    % fonction KPP
    function y=KPP(u)
        y=alpha*u.*(1-u./K);
    end    

% explicit Euler scheme
for t=dt:dt:tend
    u = u + dt .* (d.*del2(u,h,h) +  KPP(u));
 
    counter=counter+1;
    index1=index1+1;
    if(mod(counter,1000) == 0)
        subplot(2,3,ind)
        surf(x,y,u,'edgecolor','none');
        ind=ind+1;
    end
    %ustore(counter,:) = u;
    
    
    %ajout animation
    figure(2);
    surf(x,y,u,'edgecolor','none');
    drawnow;
    MOVI(index1) = getframe; % creation de l'animation
    
end
figure(1);
subplot(2,3,ind)
surf(x,y,u,'edgecolor','none');

end