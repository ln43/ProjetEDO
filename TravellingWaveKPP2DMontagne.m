function sol=TravellingWaveKPP2DMontagne()
% Solver for the Fisher/KPP equation with homogeneous Neumann boundary conditions
% in 1D using finite differences
% 
% u_t - \Delta u = u (1-u)
%


close all;
clear all;

d=0.5;
alpha=0.5;
K=0.8;
u0=0.1;

a = -100;
b = 100;
nel = 200; % number of elements
h = (b-a)/nel; % step size
nv = nel+1;% number of vertices

x = a:h:b; % mesh
y = a:h:b;
dt = 0.01; % time steps
tend = 1000; 

e=ones(nv, 1);

% Discretization of the Laplace operator with Neumann bc


% initial guess
u = zeros(201,201);
for i=1:1:201
    for j=1:1:201
        r= (i-110).^2 + (j-120).^2 ;
        if (r<10) 
            u(i,j)=u0* exp(-0.01*r);%*((i-100).^2 + (j-100).^2 ) ;
        end   
%         if i>80 & i<120 & j<120 & j>80
%             u(i,j)=u0*(i.^2 + j.^2 ) ;
%         end
    end
end
counter=1;
ind=1;
index1=1;
figure(1)
subplot(2,3,ind)
surf(x,y,u,'edgecolor','none');
axis([a b a b 0 1])
xlabel('x')
ylabel('y')
zlabel('u')
title('Condition initiale');
saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
ind=ind+1;

    % fonction KPP
    function y=KPP(u)
        mk=zeros(201,201);
        mk=mk+1/0.2;
        for l=0:5
            mk=abs(mk+(1/0.001)*diag(ones(1,201-l),l));
        end
        for l=5:10
            mk=abs(mk+(1/0.01)*diag(ones(1,201-l),l));
        end
         for l=10:15
            mk=abs(mk+(1/0.1)*diag(ones(1,201-l),l));
         end
         for l=15:20
            mk=abs(mk+(1/0.2)*diag(ones(1,201-l),l));
         end
        for l=20:25
            mk=abs(mk+(1/0.3)*diag(ones(1,201-l),l));
        end
        for l=25:30
            mk=abs(mk+(1/0.4)*diag(ones(1,201-l),l));
        end       
        for l=0:201
            mk=abs(mk+(1/0.000001)*diag(ones(1,201-l),-l));
        end   
%         mk(90:110,90:110)=(1/0.55);
%         mk(95:105,95:105)=(1/0.6);
%         mk(35:65,35:65)=1/0.45;
%         mk(40:60,40:60)=1/0.4;
%         mk(160:180,160:180)=1/0.53;
%         mk(40:60, 160:180)=1/0.48;
%         mk(160:180, 40:60)=1/0.57;
        
%         mk(50:150,50:150)=(1/K/0.2);
        %mk(75:125,75:125)=1/K;
        y=alpha*u.*(1-u.*mk);
    end    

% explicit Euler scheme
for t=dt:dt:tend
    u = u + dt .* (d.*del2(u,h,h) +  KPP(u));
   
    
    counter=counter+1;
    %if(mod(counter,10) == 0 && ind <=12)
    index1=index1+1;
    if(mod(counter,2000) ==0 )
        figure(1)
        subplot(2,4,ind)
       
        surf(x,y,u,'edgecolor','none');
        axis([a b a b 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('u')
        title(['t=' num2str(t+0.01)]);
        saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
        ind=ind+1;
    end    
    if(mod(counter,100) == 0)
        %ajout animation je sais pas trop si ca marche vraiment !!
        figure(2);
        surf(x,y,u,'edgecolor','none');
        xlabel('x')
        ylabel('y')
        zlabel('u')
        axis([a b a b 0 1]);
        title(['Propagation de l''onde avec une croissance logistique K variable, \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)])
        drawnow;
        MOVI(index1) = getframe; % creation de l'animation
    end    
    
end

figure(1)
suplabel(['Propagation de l''onde avec une croissance logistique K variable, \alpha =', num2str(alpha),', u_0 =', num2str(u0),' et d =', num2str(d)])

% 
% surf(x,y,u,'edgecolor','none');
% axis([a b a b 0 1])
% xlabel('x')
% ylabel('y')
% zlabel('u')
% saveas(figure(1),strcat('F',num2str(ind)),'jpeg')
end