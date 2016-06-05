function sol=TravellingWaveKPP2D()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
%


close all;
clear all;

d=0.5;
alpha=0.02;
K=0.5;
u0=0.01;

a = -50;
b = 50;
nel = 200; % number of elements
h = (b-a)/nel; % step size
nv = nel+1;% number of vertices

x = a:h:b; % mesh
y = a:h:b;
dt = 0.01; % time steps
tend = 100; 

% initial guess
u = zeros(201,201);
for i=1:1:201
    for j=1:1:201
        r= (i-100).^2 + (j-100).^2 ;
        if (r<10) 
            u(i,j)=u0* exp(-0.01*r);
        end   
    end
end
counter=1;
% ind=1;
index1=1;
% figure(1)
% surf(x,y,u,'edgecolor','none');
% axis([a b a b 0 1])
% xlabel('x')
% ylabel('y')
% zlabel('u')
% saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
% ind=ind+1;

    % fonction KPP
    function y=KPP(u)
        y=alpha*u.*(1-u./K);
    end    

% explicit Euler scheme
for t=dt:dt:tend
    u = u + dt .* (d.*del2(u,h,h) +  KPP(u));
    counter=counter+1;
    index1=index1+1;
%     if(mod(counter,1000) == 0)
%         figure(1)
%         subplot(2,3,ind)
%         surf(x,y,u,'edgecolor','none');
%         axis([a b a b 0 1])
%         xlabel('x')
%         ylabel('y')
%         zlabel('u')
%         saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
%         ind=ind+1;
%     end    
    if(mod(counter,100) == 0)
        figure(2);
        surf(x,y,u,'edgecolor','none');
        axis([a b a b 0 1]);
        title(['Propagation de l''onde avec une croissance logistique K =', num2str(K),', \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)])
        drawnow;
        MOVI(index1) = getframe; % creation de l'animation
    end    
    
end

figure(1)
suplabel(['Propagation de l''onde avec une croissance logistique K =', num2str(K),', \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)])

% 
% surf(x,y,u,'edgecolor','none');
% axis([a b a b 0 1])
% xlabel('x')
% ylabel('y')
% zlabel('u')
% saveas(figure(1),strcat('F',num2str(ind)),'jpeg')

end