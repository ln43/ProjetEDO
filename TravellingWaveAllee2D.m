function sol=TravellingWaveAllee2D()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
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
figure(1)
subplot(4,3,ind);
surf(x,y,u,'edgecolor','none');
ind=ind+1;

    % fonction Allee
    function y=Allee(u)
        y=k*u.*(1-u).*(u-A);
    end   

% explicit Euler scheme

%Version pour A>0.5 : tous les 2 pas de temps, au bout d'environ 25 déjà
%plus rien
for t=dt:dt:tend
    u = u + dt .* (d.*del2(u,h,h) +  Allee(u));
    counter=counter+1;
    if(mod(counter,2) == 0 && ind <=12)
        subplot(4,3,ind)
        %figure(ind)
        surf(x,y,u,'edgecolor','none');
        ind=ind+1;
    end
end

%Version pour A<0.5 : tous les 500 pas de temps
% for t=dt:dt:tend
%     u = u + dt .* (d.*del2(u,h,h) +  Allee(u));
%     counter=counter+1;
%     if(mod(counter,500) == 0)
%         subplot(4,3,ind)
%         %figure(ind)
%         surf(x,y,u,'edgecolor','none');
%         ind=ind+1;
%     end
% end

%figure(ind)
%subplot(4,3,ind)
%surf(x,y,u,'edgecolor','none');

end