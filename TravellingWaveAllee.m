
function sol=TravellingWaveAllee()
% Solver for the Allee equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = ku(1-u)(u-A)
%

close all;
clear all;

d=0.005;
A=0.25;
k=4/(1-A)^2;

u0=0.1;

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
u = u0 * (-10<=x & x<=10 ) ;
ustore(1,:) = u;
counter=1;
figure(2)
plot(x,u);
hold on

    % fonction Allee
    function y=Allee(u)
        y=k*u.*(1-u).*(u-A);
    end    

% explicit Euler scheme
for t=dt:dt:tend
    
    u = u + dt * (d*(Delta * u')' +  Allee(u));
 
    counter=counter+1;
    ustore(counter,:) = u;   
    
    if(mod(counter,1000) == 0)
        figure(2)
        plot(x,u);
        hold on
        legendinfo{counter/1000+1}=strcat('t=', num2str(counter));
    end
end

figure(1);
for A=[0.25,0.5,0.75]
    plot(linspace(0,1,50), Allee(linspace(0,1,50)),[0 1], [0 0],'g--');
    hold on
end
title('f(u)')
xlabel('u')
ylabel('f(u)')
legend('A=0.25','A=0.5','A=0.75')
hold off

figure(2)
plot(x,u);
title('Fronts d ondes pour differents temps');%tous les milles pas de temps
xlabel('t');
ylabel('u');
legend(legendinfo)
figure(3);
indX=[101,121,171,201];
for i=1:1:4
  %subplot(2,2,i) 
  plot(dt:dt:tend,ustore(2:length(ustore(:,indX(i))),indX(i)))
  hold on
  title(strcat('u(t), x =', num2str(indX(i)/2-50.5)))
  xlabel('t')
  ylabel('u')
  axis([0 tend+1 0 max(u)+0.1]);
  %axis tight
end


figure(4)
surf(x,dt:dt:tend,ustore(2:length(ustore(:,1)),:),'edgecolor','none');

xlabel('Distance x')
ylabel('Time t')
zlabel('Specie u')

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