function sol=TravellingWaveKPP()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
%


close all;
clear all;

d=0.5;
alpha=0.5;
K=0.5;
u0=0.1;

a = -50;
b = 50;
nel = 200; % number of elements
h = (b-a)/nel; % step size
nv = nel+1;% number of vertices

x = a:h:b; % mesh
dt = 0.01; % time steps
tend = 100; 

e=ones(nv, 1);

% Discretization of the Laplace operator with Neumann bc
Delta=1/(h^2)*spdiags([e -2*e e], -1:1, nv, nv);
Delta(1,2) = 2/h^2;
Delta(end,end-1) = 2/h^2;

% initial guess
 u = u0*exp(-0.01*(x.^2)) ; %
%u= u0 * ( x<-10 ) ;  
ustore(1,:) = u;
counter=1;
index1=0;
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
    index1=index1+1;
    ustore(counter,:) = u;   
    
    if(mod(counter,1000) == 0)
        figure(2)
        plot(x,u);
        hold on
        legendInfo{counter/1000+1}=['t=' num2str(t+0.01)];
    end
    
%     if(mod(counter,10) == 0)
%         figure(5);
%         plot(x,u,'green','LineWidth',2);
%         title(['Propagation de l''onde avec une croissance logistique K =', num2str(K),', \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)])
%         xlabel('x')
%         ylabel('u')
%         axis([-50 50 0 1])
%         drawnow;
%         MOVI(counter) = getframe; % creation de l'animation
%         hold off;
%     end    

    
end

figure(1);
plot(linspace(0,1,50), KPP(linspace(0,1,50)),'g',[0 1], [0 0],'--');
ylabel('f(u)')
xlabel('u')
title(['f(u) avec une croissance logistique K =', num2str(K),', \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)]);%tous les milles pas de temps

figure(2)
title(['Fronts d''ondes avec une croissance logistique K =', num2str(K),', \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)]);%tous les milles pas de temps
xlabel('x')
ylabel('u')
ylim([0 max(max(ustore))+0.05])
legend(legendInfo)

figure(4)
surf(x,dt:dt:tend,ustore(2:length(ustore(:,1)),:),'edgecolor','none');
xlabel('Distance x')
ylabel('Time t')
zlabel('Specie u')
title(['Fronts d''onde avec une croissance logistique K =', num2str(K),', \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)]);%tous les milles pas de temps

end