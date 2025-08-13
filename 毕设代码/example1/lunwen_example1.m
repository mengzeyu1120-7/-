  
warning off all; clear
tic
ni=58; nb=2; nc=60; nt=20;%N=60， 边界2， 内部58
tau=0.01;c=1;
u0= @(x) sin(pi*x);
mu0= @(x) sin(pi*x);
ureal= @(x,t) exp(t).*sin(pi*x);
mureal= @(x,t) exp(t).*sin(pi*x);
F_1= @(x,t) exp(2*t).*(sin(pi*x).^2)./(exp(t).*sin(pi*x)+1);
F_2= @(x,t) exp(2*t).*(sin(pi*x).^2)./(exp(t).*sin(pi*x)+1);
phi = @(r,c) sqrt(c.^2.+r.^2); 
xi=linspace(0,1,ni+2)';
xi=xi(2:ni+1);%内部点Xi
xb=[0,1]';%边界点Xb
xy=[xb;xi];
xc=linspace(0,1,nc)';%中心点y
xt=linspace(0,1,nc)';
DM1=DistanceMatrix(xi,xc);%计算内部r
DM2=DistanceMatrix(xb,xc);%计算边界r
DM3=DistanceMatrix(xt,xc);
m1=repmat(xi,1,nc);m2=repmat(xc',ni,1);%size：58*60 
m3=repmat(xt,1,nc);m4=repmat(xc',nc,1);%size：60*60
m5=repmat(xb,1,nc);m6=repmat(xc',nb,1);%size：2*60
mi1=m1-m2;%（xi-xc'）
mc1=m3-m4;%（xt-xc'）
mb1=m5-m6;
Phi=phi(DM1,c);
un=zeros(nc,10000);%解的矩阵
mun=zeros(nc,10000);
un(:,1)= u0(xt);
uni(:,1)= un((2:ni+1),1);
mun(:,1)= mu0(xt);
Ab=phi(DM2,c);        Ob=zeros(nb,nc); Obv=zeros(nb,1);
A1=(2-tau)*phi(DM1,c);       A2=tau*(uni(:,1)./(uni(:,1)+1)).*phi(DM1,c);
A3=zeros(ni,nc);             A4=(2-tau+tau*(uni(:,1)./(uni(:,1)+1))).*phi(DM1,c);
B1v=(tau+2)*u0(xi);          B2v=-tau*(uni(:,1)./(uni(:,1)+1)).*mu0(xi);
B3v=zeros(ni,1);             B4v=(-tau*(uni(:,1)./(uni(:,1)+1))+tau+2).*mu0(xi);
t_0=0*ones(ni);
t_1=tau*ones(ni);
F1=(F_1(xi,t_0)+F_1(xi,t_1))/2;
F2=(F_2(xi,t_0)+F_2(xi,t_1))/2;
%F1=F_1(xi,0);
%F2=F_2(xi,tau);
bn=[zeros(nb,1);zeros(nb,1);2*tau*F1(:,1);2*tau*F2(:,1)];
A=[Ab Ob;Ob Ab;A1 A2;A3 A4];%A
tol=det(A);
Bv=[Obv;Obv;B1v+B2v;B3v+B4v];%B0
v=zeros(2*nc,5000);
% v=A\(Bv+bn);
v(:,2)=pinv(A)*(Bv+bn);%求伪逆效果更好 得到v1
u_0=u0(xt);
u1=phi(DM3,c)*v((1:nc),2);
un(:,2)= phi(DM3,c)*v((1:nc),2);
uni(:,2)= un(2:ni+1,2);
mun(:,2)= phi(DM3,c)*v((nc+1:2*nc),2);
%un(:,2)= (un(:,2)+mun(:,2))/2;
%mun(:,2)= (un(:,2)+mun(:,2))/2;

tend=1;
for i=2:floor(tend/tau)
%un(:,i)=(un(:,i-1)+un(:,i))/2;
%uni(:,2)= un(2:ni+1,2);
%mun(:,i)=(mun(:,i-1)+mun(:,i))/2;

Ab=phi(DM2,c);        Ob=zeros(nb,nc); 
A1=(2-tau)*phi(DM1,c);       A2=tau*(uni(:,i)./(uni(:,i)+1)).*phi(DM1,c);
A3=zeros(ni,nc);             A4=(2-tau+tau*(uni(:,i)./(uni(:,i)+1))).*phi(DM1,c);
B1=(tau+2)*phi(DM1,c);       B2=-tau*(uni(:,i)./(uni(:,i)+1)).*phi(DM1,c);
B3=zeros(ni,nc);             B4=(-tau*(uni(:,i)./(uni(:,i)+1))+tau+2).*phi(DM1,c);
A=[Ab Ob;Ob Ab;A1 A2;A3 A4];
B=[Ob Ob;Ob Ob;B1 B2;B3 B4];

t_n=tau*(i-1)*ones(ni);
t_n1=tau*i*ones(ni);
F1=(F_1(xi,t_n)+F_1(xi,t_n1))/2;
F2=(F_2(xi,t_n)+F_2(xi,t_n1))/2;
%F1=F_1(xi,t_n);
%F2=F_2(xi,t_n);
bn=[zeros(nb,1);zeros(nb,1);2*tau*F1(:,1);2*tau*F2(:,1)];
v(:,i+1)=pinv(A)*(B*v(:,i)+bn);
un(:,i+1)=phi(DM3,c)*v((1:nc),i+1);
mun(:,i+1)=phi(DM3,c)*v((nc+1:2*nc),i+1);
uni(:,i+1)=un(2:ni+1,i+1);
end

unreal= zeros(nc,10000);
munreal= zeros(nc,10000);
for i=1:floor(tend/tau)
unreal(:,i)= ureal(xt,tau*(i-1));
munreal(:,i)= mureal(xt,tau*(i-1));
end

E_1= zeros(floor(tend/tau));
E_2= zeros(floor(tend/tau));
for i=2:floor(tend/tau)
E_1(i)= norm(un(:,i)-unreal(:,i),inf);
E_2(i)= norm(mun(:,i)-munreal(:,i),inf);
end

fprintf('Eu1 = %8.3e\n', E_1(2));
fprintf('Emu1 = %8.3e\n', E_2(2));
fprintf('Eun = %8.3e\n', E_1(floor(tend/tau)));
fprintf('Emun = %8.3e\n', E_2(floor(tend/tau)));
Eu_1=norm(E_1(1:floor(tend/tau)),inf);
fprintf('Eu = %8.3e\n', Eu_1);
Eu_2=norm(E_2(1:floor(tend/tau)),inf);
fprintf('Emu = %8.3e\n', Eu_2);
Ec=norm((un(:,(floor(tend/tau)))- un(:,(floor(tend/tau)-1))),inf)/norm(un(:,(floor(tend/tau))),inf);
fprintf('Ec = %8.3e\n', Ec);
plot(xt,un(:,1),'-ok');
hold on
plot(xt,un(:,(floor(0.25/tau))),'-sr');
hold on
plot(xt,un(:,(floor(0.5/tau))),'-g.');
hold on
plot(xt,un(:,(floor(0.75/tau))),'-+m');
hold on
plot(xt,un(:,(floor(tend/tau))),'-dy');
legend('t=0','t=0.25','t=0.5','t=0.75','t=1')
figure;

plot(xt,mun(:,1),'-ok');
hold on
plot(xt,mun(:,(floor(0.25/tau))),'-sr');
hold on
plot(xt,mun(:,(floor(0.5/tau))),'-g.');
hold on
plot(xt,mun(:,(floor(0.75/tau))),'-+m');
hold on
plot(xt,mun(:,(floor(tend/tau))),'-dy');
legend('t=0','t=0.25','t=0.5','t=0.75','t=1')
figure;

plot(xt,unreal(:,1),'-ok');
hold on
plot(xt,unreal(:,(floor(0.25/tau))),'-sr');
hold on
plot(xt,unreal(:,(floor(0.5/tau))),'-g.');
hold on
plot(xt,unreal(:,(floor(0.75/tau))),'-+m');
hold on
plot(xt,unreal(:,(floor(tend/tau))),'-dy');
legend('t=0','t=0.25','t=0.5','t=0.75','t=1')
figure;

plot(xt,munreal(:,1),'-ok');
hold on
plot(xt,munreal(:,(floor(0.25/tau))),'-sr');
hold on
plot(xt,munreal(:,(floor(0.5/tau))),'-g.');
hold on
plot(xt,munreal(:,(floor(0.75/tau))),'-+m');
hold on
plot(xt,munreal(:,(floor(tend/tau))),'-dy');
legend('t=0','t=0.25','t=0.5','t=0.75','t=1')
figure;

plot(xt,un(:,2),'-ok');
hold on
plot(xt,mun(:,2),'-g');
hold on
plot(xt,unreal(:,2),'-sr');
legend('u','mu','real')
figure;

plot(xt,un(:,(floor(tend/tau))),'-ok');
hold on
plot(xt,mun(:,(floor(tend/tau))),'-g');
hold on
plot(xt,unreal(:,(floor(tend/tau))),'-sr');
legend('u','mu','real')
figure;

N= linspace(1,100,100)';
plot(N,E_1);
figure;
plot(N,E_2)
