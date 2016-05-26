function sol=TravellingWaveKPP()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
%


close all;
clear all;

d=0.5;
alpha=0.05;
K=0.3;
u0=0.7;

a = -50;
b = 50;
nel = 200; % number of elements
h = (b-a)/nel; % step size
nv = nel+1;% number of vertices

x = a:h:b; % mesh
dt = 0.01; % time steps
tend = 50; 

e=ones(nv, 1);

% Discretization of the Laplace operator with Neumann bc
Delta=1/(h^2)*spdiags([e -2*e e], -1:1, nv, nv);
Delta(1,2) = 2/h^2;
Delta(end,end-1) = 2/h^2;

% initial guess
u = u0 * (-10<x & x<10 ) ;
ustore(1,:) = u;
counter=1;
figure(2)
plot(x,u);
hold on
legendInfo{counter}=['t=' num2str(counter)];


    % fonction KPP
    function y=KPP(u)
        y=alpha*u.*(1-u./K);
    end    

% explicit Euler scheme
for t=dt:dt:tend
    
    u = u + dt * (d*(Delta * u')' +  KPP(u));
 
    counter=counter+1;
    ustore(counter,:) = u;   
    
    if(mod(counter,1000) == 0)
        figure(2)
        plot(x,u);
        hold on
        legendInfo{counter/1000+1}=['t=' num2str(counter)];
    end
end

figure(1);
plot(linspace(0,1,50), KPP(linspace(0,1,50)),'g',[0 1], [0 0],'--');
title('f(u)')
xlabel('u')
ylabel('f(u)')

figure(2)
title('Fronts d ondes pour differents temps');%tous les milles pas de temps
xlabel('x')
ylabel('u')
ylim([0 max(max(ustore))+0.05])
legend(legendInfo)

figure(3);
indX=[101,121,171,201];
for i=1:1:4
  plot(dt:dt:tend,ustore(2:length(ustore(:,indX(i))),indX(i)))
  hold on
end
hold off
xlabel('t')
ylabel('u')
axis([0 tend+1 0 max(max(ustore))+0.05]);


figure(4)
surf(x,dt:dt:tend,ustore(2:length(ustore(:,1)),:),'edgecolor','none');

xlabel('Distance x')
ylabel('Time t')
zlabel('Specie u')

%%Je ne sais pas ce que �a affiche ???
% [xmesh, tmesh] = meshgrid(x,0:dt:tend);
% zmesh = xmesh-2 * diag(0:dt:tend) * ones(size(xmesh));

% figure(2)
% contour(xmesh,tmesh, ustore);
% hold on
% plot([10 20 20 10],[5 10 5 5],'--');
% text(0,3,'Ref triangle with slope 2');
% xlabel('x');
% ylabel('t');
% title('Contour plot of u(x,t)')
% hold off

%figure(3)
%waterfall(zmesh, tmesh, ustore)
%zlim([0 1])
%xlim([-120,50])
%view (50,30)
end