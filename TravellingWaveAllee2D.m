function sol=TravellingWaveAllee2D()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
%

close all;
clear all;

d=0.5;
A=0.25;
k=4/(1-A)^2;
u0=0.5;

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
index1=1;
ind=1;
figure(1)
surf(x,y,u,'edgecolor','none');
axis([a b a b 0 1])
xlabel('x')
ylabel('y')
zlabel('u')
saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
ind=ind+1;

    % fonction Allee
    function y=Allee(u)
        y=k*u.*(1-u).*(u-A);
    end   

% explicit Euler scheme

%Version pour A<0.5 : tous les 500 pas de temps
for t=dt:dt:tend
    u = u + dt .* (d.*del2(u,h,h) +  Allee(u));
    counter=counter+1;
    index1=index1+1;
    %if(mod(counter,10) == 0 && ind <=12)
    if(mod(counter,500) == 0)
        surf(x,y,u,'edgecolor','none');
        axis([a b a b 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('u')
        saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
        ind=ind+1;
    end 
    if(mod(counter,10) == 0)
        %ajout animation je sais pas trop si ca marche vraiment !!
        figure(2);
        surf(x,y,u,'edgecolor','none');
        axis([a b a b 0 1]);
        drawnow;
        MOVI(index1) = getframe; % creation de l'animation
    end    
end

surf(x,y,u,'edgecolor','none');
axis([a b a b 0 1])
xlabel('x')
ylabel('y')
zlabel('u')
saveas(figure(1),strcat('F',num2str(ind)),'jpeg')

end