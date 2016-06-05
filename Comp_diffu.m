
%%%% Vérifier la conditions initiale %%%%% Dans le modèle on trouve : N --> u    et   u---> psy  


%%%%%%% Important %%%%%%  Chercher l'équilibre et le mettre dans le bord

clear all;
clc;

% parameters of the model
alpha1=0.6;
alpha2=0.6;
K1=0.5;
K2=0.5;
gama1=0.5;
gama2=1.5;
d1=1;
d2=0.5;

% discritisation
L=80; T=20;
dx=0.1;dt=0.05;
x=-L:dx:L;M_max=length(x);
t=0:dt:T;N_max=length(t);

% Equilibres
 equ=(K1-gama1*K2)/(1-gama1*gama2);
 eqv=(K2-gama2*K1)/(1-gama1*gama2);
 
% initiation 
u=zeros(N_max,M_max);
v=zeros(N_max,M_max);
pqy=zeros(N_max,M_max);

%equilibre aux bords
u(:,M_max)=0;
v(:,M_max)=0;


% Condition initial pour u et v 

for i=1:M_max
    u(1,i)=exp(-0.01*x(i)^2);
%     if x(i)>-30 & x(i)<30
%        u(1,i)=equ;
%     elseif x(i)<=-30 | x(i)>=30
%         u(1,i)=0;
%     else u(1,i)=(equ/60)*x(i)+equ/2;
%     end
end

for i=1:M_max
    v(1,i)=exp(-0.01*x(i)^2);
%     if x(i)>-30 & x(i)<30
%        v(1,i)=eqv;
%     elseif x(i)<=-30 | x(i)>=30
%         v(1,i)=0;
%     else v(1,i)=(eqv/60)*x(i)+eqv/2;
%     end
end

% tracer la condition initiale de u ou v
figure(2);
plot(x,u(1,:))
hold on;
plot(x,v(1,:))

% schema central (dérivée seconde)

 a=zeros(M_max-2,M_max-2);
     for i=1:M_max-2
     for j=1:M_max-2      
              if (i==j) 
                  a1(i,j)=1+2*d1*dt/power(dx,2);
                 else
                  if i==j+1||j==i+1 
                      a1(i,j)=-d1*dt/power(dx,2) ;
                 
                  end
              end
     end
     end
     
  a2=zeros(M_max-2,M_max-2);
     for i=1:M_max-2
     for j=1:M_max-2      
              if (i==j) 
                  a2(i,j)=1+2*d2*dt/power(dx,2);
                 else
                  if i==j+1||j==i+1 
                      a2(i,j)=-d2*dt/power(dx,2) ;
                 
                  end
              end
     end
     end    



for n=2:N_max

for i=1:M_max-2

    % non linearité 
    
 z1=alpha1*u(n-1,i+1)*(1-u(n-1,i+1)/K1-gama1*v(n-1,i+1)/K1);
 z2=alpha2*v(n-1,i+1)*(1-v(n-1,i+1)/K2-gama2*u(n-1,i+1)/K2);
       
          if i==1   
              b1(i)=dt*z1+u(n-1,i+1)+d1*dt/power(dx,2)*0 ;
          else if i==M_max-2 
                  b1(i)=dt*z1+u(n-1,i+1)+d1*dt/power(dx,2)*u(n,M_max);
          else b1(i)=z1*dt+u(n-1,i+1);
              end
          end
          
           if i==1   
              b2(i)=dt*z2+v(n-1,i+1)+d2*dt/power(dx,2)*0 ;
          else if i==M_max-2 
                  b2(i)=dt*z2+v(n-1,i+1)+d2*dt/power(dx,2)*v(n,M_max);
          else b2(i)=z2*dt+v(n-1,i+1);
              end
          end
end
uuu=a1\b1';
vvv=a2\b2';
%  for i=2:M_max-1
%    u(n,i)=uuu(i-1);
%  end
u(n,2:M_max-1)= uuu;
v(n,2:M_max-1)= vvv;
end


wu=u(:,1:M_max-100);
wv=v(:,1:M_max-100);


subplot(1,2,1);mesh(x(1:M_max-100),t,wu);
ylabel('t');
xlabel('Distance x')
zlabel('Densité espèce Homo Sapiens')
title(' Homo Sapiens  ');

subplot(1,2,2);mesh(x(1:M_max-100),t,wv);
xlabel('Distance x');
ylabel('t');
zlabel('Densité espèce Néandertal ');
title(' Néandertal ');

suplabel(['Propagation de l''onde avec compétition  \alpha_1=', num2str(alpha1),',  \alpha_2=', num2str(alpha2),',  \gamma_1=', num2str(gama1),',  \gamma_2=', num2str(gama2),',  K_1=', num2str(K1),',  K_2=', num2str(K2),',  d_1=', num2str(d1),' et  d_2=', num2str(d2)])

