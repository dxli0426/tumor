randn('state',64) 
T=3000;
dt=0.1;
N=T/dt;
gw=randn(N,1);
x=zeros(N,1);
t=0:dt:(T-dt);
a=1;b=2;
miua=0.9;
miub=0.6;
r=zeros(1,N);
rzero=a;
p1=miub/(miua+miub)+(miua/(miua+miub))*exp(-(miua+miub)*dt);
p2=miub/(miua+miub)-(miub/(miua+miub))*exp(-(miua+miub)*dt);
R=unifrnd(0,1,1,N);
if p1>R(1)
   r(1)=a;
else r(1)=b;
end
for i =2:N
    if r(i-1)==a
       if p1>R(i)
          r(i)=a;
       else r(i)=b;
       end
    else  r(i-1)=b;
          if p2>R(i)
             r(i)=a;
          else r(i)=b;
          end
    end
end
x(1)=0.1;
theta(1)=0.81;theta(2)=0.81;
beta(1)=0.95;beta(2)=0.95;
for j=2:N
    sigma=0.9;
     x(j)=x(j-1)+dt*(x(j-1)*(1-x(j-1)*theta(r(j)))-(beta(r(j))*x(j-1))/(1+x(j-1)))...
        +sqrt(dt)*sigma*x(j-1)^2*gw(j-1)+0.5*sigma*x(j-1)^4*sqrt(dt)*(gw(j-1)^2-1);
end
hold on
plot(t,x,'r-','linewidth',2.3)
box on;
set(gca,'YTick',0:0.05:0.9)
axis([0 30 0 0.9])
legend('\theta_{\rho(t)}=0.81,\beta_{\rho(t)}=0.95')
xlabel('time')
title('Weak persistence of tumor cells')
ylabel('Normalized tumor cell number')

