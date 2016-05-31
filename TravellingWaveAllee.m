
function sol=TravellingWaveAllee()
% Solver for the Allee equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = ku(1-u)(u-A)
%

close all;
clear all;

d=0.05;
A=0.75;
k=4/(1-A)^2;
u0=0.9;

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
u = u0*exp(-0.01*(x.^2)) ;
ustore(1,:) = u;
counter=1;
index1=0;
figure(2)
plot(x,u);
legendInfo{1}=['t=' num2str(0)];
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
    index1=index1+1;
    if(mod(counter,1000) == 0)
        figure(2)
        plot(x,u);
        hold on
        legendInfo{counter/1000+1}=['t=' num2str(counter)];
    end
%     if(mod(counter,50) == 0)
%         figure(4);
%         plot(x,u,'green','LineWidth',2);
%         axis([a b 0 1])
%         drawnow;
%         MOVI(index1) = getframe; % creation de l'animation
%         hold off;
%     end
    
end
hold off

% figure(1);
% for A=[0.25,0.5,0.75]
%     plot(linspace(0,1,50), Allee(linspace(0,1,50)),[0 1], [0 0],'g--');
%     hold on
% end
% title('f(u)')
% xlabel('u')
% ylabel('f(u)')
% legend('A=0.25','A=0.5','A=0.75')
% hold off

figure(2)
title('Fronts d ondes pour différents temps');%tous les milles pas de temps
xlabel('t');
ylabel('u');
ylim([0 max(max(ustore))+0.05])
legend(legendInfo)

figure(4)
surf(x,dt:dt:tend,ustore(2:length(ustore(:,1)),:),'edgecolor','none');
xlabel('Distance x')
ylabel('Time t')
zlabel('Specie u')
end