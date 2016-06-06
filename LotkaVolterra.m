function LotkaVolterra
clear; close all;

L = 100; maxt = 100;
P(1) = 0.6; %alpha1
P(2) = 0.6; %alpha2
P(3) = 0.2; %K1
P(4) = 0.2; %K2
P(5) = 0.5; %gamma1
P(6) = 0.5; %gamma2
P(7) = 0.5; %d1
P(8) = 0.5; %d2

m = 0;
t = linspace(0,maxt,200);
x = linspace(0,L,200);

sol = pdepe(m,@LVfun,@ICfun,@BCfun,x,t,[],P);
u = sol(:,:,1);
v = sol(:,:,2);

subplot(1,2,1)
surf(x,t,u,'edgecolor','none')
axis([0 L 0 maxt 0 max(max(u))])
ylabel('t');
xlabel('Distance x')
zlabel('Densité espèce Homo Sapiens')
title('Homo Sapiens')

subplot(1,2,2)
surf(x,t,v,'edgecolor','none')
axis([0 L 0 maxt 0 max(max(v))])
xlabel('Distance x');
ylabel('t');
zlabel('Densité espèce Néandertal ');
title('Homo Neanderthalensis')


    function [c,f,s] = LVfun(x,t,u,dudx,P)
        c = [1;1];
        f = [P(7);P(8)].*dudx;
        s = [P(1)*u(1)*(1-u(1)/P(3)-P(5)/P(3)*u(2));
        P(2)*u(2)*(1-u(2)/P(4)-P(6)/P(4)*u(1))];      
    end

    function u0 = ICfun(x,P)
        ui=0.03;
        vi=0.2;
        u0 = [ui;vi];
    end

    function [pl,ql,pr,qr] = BCfun(xl,ul,xr,ur,t,P)
        pl = [0;0]; ql = [1;1];
        pr = [0;0]; qr = [1;1];
    end

end