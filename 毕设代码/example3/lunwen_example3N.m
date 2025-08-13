warning off all; clear
tic
Eu= zeros(5);
Emu= zeros(5);
NN= [2 3 16 18 23];
for j=1:5
ni=NN(j)^2; nb=4*(NN(j)+1); nc=ni+nb;
tau=0.01;c=1;
u0= @(x,y) sin(pi*x).*sin(pi*y);
u0_1= @(x,y) pi.*cos(pi*x).*sin(pi*y)+pi.*sin(pi*x).*cos(pi*y);
u0_2= @(x,y) -2*(pi^2).*sin(pi*x).*sin(pi*y)+2*(pi^2).*cos(pi*x).*cos(pi*y);
mu0= @(x,y) sin(pi*x).*sin(pi*y);
mu0_1= @(x,y) pi.*cos(pi*x).*sin(pi*y)+pi.*sin(pi*x).*cos(pi*y);
mu0_2= @(x,y) -2*(pi^2).*sin(pi*x).*sin(pi*y)+2*(pi^2).*cos(pi*x).*cos(pi*y);
ureal= @(x,y,t) exp(t).*sin(pi*x).*sin(pi*y);
mureal= @(x,y,t) exp(t).*sin(pi*x).*sin(pi*y);
F_1= @(x,y,t) exp(2*t).*(sin(pi*x).^2).*(sin(pi*y).^2)./(exp(t).*sin(pi*x).*sin(pi*y)+1)+(2*pi^2+1).*exp(t).*sin(pi*x).*sin(pi*y)-(2*pi^2).*exp(t).*cos(pi*x).*cos(pi*y);
F_2= @(x,y,t) exp(2*t).*(sin(pi*x).^2).*(sin(pi*y).^2)./(exp(t).*sin(pi*x).*sin(pi*y)+1)+(2*pi^2+1).*exp(t).*sin(pi*x).*sin(pi*y)-(2*pi^2).*exp(t).*cos(pi*x).*cos(pi*y);
phi = @(r,c) (c.^2.+r.^2).^(1/2); 
phi_r = @(r,c) ((c.*c+r.*r).^(-1/2));%部分
%phi_rr= @(r,c) c^2.*((c.*c+r.*r).^(-3/2));
phi_r1=@(r,c) ((c.*c+r.*r).^(-1/2));
phi_r2=@(r,c) ((c.*c+r.*r).^(-3/2));

% 边界点
c1=zeros(nb/4,1);c2=linspace(0,1,nb/4)';
c3=ones(nb/4,1);c4=linspace(0,1,nb/4)';
r1=linspace(0,1,nb/4)';r2=zeros(nb/4,1);
r3=linspace(0,1,nb/4)';r4=ones(nb/4,1);
xb=[c1;c2;c3;c4];
yb=[r1;r2;r3;r4];

%内部点
nig=sqrt(ni);
xi1=linspace(0,1,nig);yi11=linspace(0,1,nig)';%yi12=linspace(1,0,nt)';
xi2=repmat(xi1,nig,1);yi2=repmat(yi11,1,nig);
xi=reshape(xi2,[ni,1]);yi=reshape(yi2,[ni,1]);

%中心点
xc=[xb;xi];yc=[yb;yi];

DM1=DistanceMatrix([xi yi],[xc yc]);% 内点与中心点的距离
DM2=DistanceMatrix([xb yb],[xc yc]);% 边界点与中心点的距离
DM3=DistanceMatrix([xc yc],[xc yc]);% 测试点与中心点的距离 

m11=repmat(xi,1,nc);m12=repmat(xc',ni,1);
m21=repmat(yi,1,nc);m22=repmat(yc',ni,1);
m31=repmat(xc,1,nc);m32=repmat(xc',nc,1);
m41=repmat(yc,1,nc);m42=repmat(yc',nc,1);
m51=repmat(xb,1,nc);m52=repmat(xc',nb,1);
m61=repmat(yb,1,nc);m62=repmat(yc',nb,1);
mi1=m11-m12;mi2=m21-m22;
mc1=m31-m32;mc2=m41-m42;
mb1=m51-m52;mb2=m61-m62;

Phi=phi(DM1,c);
Phi_rr=(2.*phi_r1(DM1,c)-(mi1.^2).*phi_r2(DM1,c)-(mi2.^2).*phi_r2(DM1,c)-2*(mi1.*mi2).*phi_r2(DM1,c));
un=zeros(nc,10000);%解的矩阵
mun=zeros(nc,10000);
uni=zeros(ni,10000);
un(:,1)= u0(xc,yc);
uni(:,1)= un((nb+1:nc),1);
mun(:,1)= mu0(xc,yc);
Ab=phi(DM2,c);        Ob=zeros(nb,nc); Obv=zeros(nb,1);
A1=2*phi(DM1,c)-tau*Phi_rr;              A2=tau*(uni(:,1)./(uni(:,1)+1)).*phi(DM1,c);
A3=zeros(ni,nc);                         A4=(2+tau*(uni(:,1)./(uni(:,1)+1))).*phi(DM1,c)-tau*Phi_rr;
B1v=2*u0(xi,yi)+tau*u0_2(xi,yi);         B2v=-tau*(uni(:,1)./(uni(:,1)+1)).*mu0(xi,yi);
B3v=zeros(ni,1);                         B4v=(-tau*(uni(:,1)./(uni(:,1)+1))+2).*mu0(xi,yi)+tau*mu0_2(xi,yi);
t_0=0;
t_1=tau;
F1=(F_1(xi,yi,t_0)+F_1(xi,yi,t_1))/2;
F2=(F_2(xi,yi,t_0)+F_2(xi,yi,t_1))/2;
%F1=F_1(xi,0);
%F2=F_2(xi,0);
bn=[zeros(nb,1);zeros(nb,1);2*tau*F1(:,1);2*tau*F2(:,1)];
A=[Ab Ob;Ob Ab;A1 A2;A3 A4];%A
tol=det(A);
Bv=[Obv;Obv;B1v+B2v;B3v+B4v];%B0
v=zeros(2*nc,5000);
% v=A\(Bv+bn);
v(:,2)=pinv(A)*(Bv+bn);%求伪逆效果更好 得到v1
u_0=u0(xc,yc);
u1=phi(DM3,c)*v((1:nc),2);
un(:,2)= phi(DM3,c)*v((1:nc),2);
uni(:,2)= un(nb+1:nc,2);
mun(:,2)= phi(DM3,c)*v((nc+1:2*nc),2);


tend=1;
for i=2:floor(tend/tau)
%un(:,i)=(un(:,i-1)+un(:,i))/2;
%uni(:,2)= un(2:ni+1,2);
%mun(:,i)=(mun(:,i-1)+mun(:,i))/2;

Ab=phi(DM2,c);        Ob=zeros(nb,nc); 
A1=2*phi(DM1,c)-tau*Phi_rr;              A2=tau*(uni(:,i)./(uni(:,i)+1)).*phi(DM1,c);
A3=zeros(ni,nc);                         A4=(2+tau*(uni(:,i)./(uni(:,i)+1))).*phi(DM1,c)-tau*Phi_rr;
B1=2*phi(DM1,c)+tau*Phi_rr;              B2=-tau*(uni(:,i)./(uni(:,i)+1)).*phi(DM1,c);
B3=zeros(ni,nc);                         B4=(-tau*(uni(:,i)./(uni(:,i)+1))+2).*phi(DM1,c)+tau*Phi_rr;
A=[Ab Ob;Ob Ab;A1 A2;A3 A4];
B=[Ob Ob;Ob Ob;B1 B2;B3 B4];

t_n=tau*(i-1);
t_n1=tau*i;
F1=(F_1(xi,yi,t_n)+F_1(xi,yi,t_n1))/2;
F2=(F_2(xi,yi,t_n)+F_2(xi,yi,t_n1))/2;
%F1=F_1(xi,t_n);
%F2=F_2(xi,t_n);
bn=[zeros(nb,1);zeros(nb,1);2*tau*F1(:,1);2*tau*F2(:,1)];
v(:,i+1)=pinv(A)*(B*v(:,i)+bn);
un(:,i+1)=phi(DM3,c)*v((1:nc),i+1);
mun(:,i+1)=phi(DM3,c)*v((nc+1:2*nc),i+1);
uni(:,i+1)=un(nb+1:nc,i+1);
end

unreal= zeros(nc,10000);
munreal= zeros(nc,10000);
uni=zeros(ni,10000);
for i=1:floor(tend/tau)
unreal(:,i)= ureal(xc,yc,tau*(i-1));
munreal(:,i)= mureal(xc,yc,tau*(i-1));
end

E_1= zeros(floor(tend/tau));
E_2= zeros(floor(tend/tau));
for i=2:floor(tend/tau)
E_1(i)= norm(un(:,i)-unreal(:,i),inf);
E_2(i)= norm(mun(:,i)-munreal(:,i),inf);
end

Eu(j)=norm(E_1(1:floor(tend/tau)),inf);
Emu(j)=norm(E_2(1:floor(tend/tau)),inf);
fprintf('Eu = %8.3e\n', Eu(j));
fprintf('Emu = %8.3e\n', Emu(j));
end

x= [14 25 324 400 625];
plot(x,Eu);
figure;

plot(x,Emu);