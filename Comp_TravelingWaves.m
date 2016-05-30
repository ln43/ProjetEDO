% This function finds traveling-wave solutions for the
% Pred-Prey model. To use it correctly, one must change
% the parameters in each function defined in this file.
function ppTravelingWaves

% Solver for a diffusive version of the Fitz-Hugh-Nagumo system in 1D using
% finite differences
% 
% u_t - eps du \Delta u = 1/eps (u (1-u)(u-alpha)-v)
% v_t = gamma u - beta v
%

close all;
clear all;

a =-100; %pour enlever avance onde des deux cotés mettre a=0
b = 100;
nel = 200; % number of elements
h = (b-a)/nel; % step size
nv = nel+1;% number of vertices

x = a:h:b; % mesh
dt = 0.1; % time steps
tend = 100; 

e=ones(nv, 1);
%Discretization of the Laplace operator
Delta=1/(h^2)*spdiags([e -2*e e], -1:1, nv, nv);
Delta(1,2) = 2/h^2;
Delta(end,end-1) = -2/h^2;

% parameters of the model
alpha1=0.6;
alpha2=0.6;
K1=0.2;
K2=0.2;
gamma1=0.5;
gamma2=1.5;
d1=1;
d2=0.5;

% initial guess
u = 0.03 * (x < 3).*(x>-3 ); %pour enlever avance onde des deux cotés enlever .*(x>-3 ) dans les deux
v = 0.2 * (x < 3).* (x> -3) ;

counter = 0;
index1=1;
% explicit Euler scheme
for t=dt:dt:tend
    index1 = index1 + 1;
    counter = counter +1 ;
    ustore(counter,:) = u;
    u = u + dt * d1* ((Delta * u')') +  dt*alpha1.*u.*(1-u/K1-gamma1.*v/K1);
    v = v +  dt * d2* ((Delta * v')') + dt * alpha2.*v.*(1-v/K2-gamma2.*u/K2);
    
    if(mod(counter,50) == 0)
        figure(1)
        plot(x,u);
        hold on;
        
        figure(2);
        plot(x,v);
        hold on; 
    end
    
    

    


%figure(3)
%plot(x,u,x,v,'r');
%legend('u','v')
figure(3);
plot(u,'blue','LineWidth',2);
hold on;
plot(v,'green','LineWidth',2);
hold off;
axis([0 200 0 0.2]);                % echelle des axes
drawnow;
MOVI(index1) = getframe; % creation de l'animation
end
end