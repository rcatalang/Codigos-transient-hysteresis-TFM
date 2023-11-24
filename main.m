%% Cargamos las variables

%Condición inicial IPTG=10^-7
load('ci7_p1_i1_500.mat') %i1=logspace(-6,-5,20);
load('ci7_p2_i1_500.mat') %i1=logspace(-6,-5,20);

load('ci7_p1_i2_500.mat') %i2=logspace(-5,-4,20);
load('ci7_p2_i2_500.mat') %i2=logspace(-5,-4,20);

load('ci7_p1_i3_500.mat') %i3=logspace(-4,-3,20);
load('ci7_p2_i3_500.mat') %i3=logspace(-4,-3,20);

% condición inicial IPTG=10^-2
load('ci2_p1_i1_500.mat') %i1=logspace(-6,-5,20);
load('ci2_p2_i1_500.mat') %i1=logspace(-6,-5,20);

load('ci2_p1_i2_500.mat') %i2=logspace(-5,-4,20);
load('ci2_p2_i2_500.mat') %i2=logspace(-5,-4,20);

load('ci2_p1_i3_500.mat') %i3=logspace(-4,-3,20);
load('ci2_p2_i3_500.mat') %i3=logspace(-4,-3,20);

load('ci2_p1_2000.mat')
load('ci2_p2_2000.mat')

%% Representación de gráficas

% Creamos los vectores de señal y tiempo
%ia=logspace(-7,-6,10);
i1=logspace(-6,-5,20);
i2=logspace(-5,-4,20);
i3=logspace(-4,-3,20);
iy=i2(18:20);

%ib=logspace(-3,-2,10);
sig=[i1 i2 i3];
% sig=[ia i1 i2 i3 ib];

time=[25,50,100,150,200,400,500];

% Pasamos los vectores a formato matricial para representar los diagramas
% en 3D
s=zeros(length(time),length(sig));
for i=1:length(sig)
    s(:,i)=sig(i);
end

t=zeros(length(time),length(sig));
for i=1:length(time)
    t(i,:)=time(i);
end

%% Ploteamos

ci2_p1=[ci2_p1_i1_500;ci2_p1_i2_500;ci2_p1_i3_500];
ci2_p2=[ci2_p2_i1_500;ci2_p2_i2_500;ci2_p2_i3_500];
ci7_p1=[ci7_p1_i1_500;ci7_p1_i2_500;ci7_p1_i3_500];
ci7_p2=[ci7_p2_i1_500;ci7_p2_i2_500;ci7_p2_i3_500];

figure(1)
%Cond inicial 2
plot(sig,ci2_p1(:,5),'r:','DisplayName','t=5,ci 1e-2')
hold on
plot(sig,ci2_p1(:,50),'r--','DisplayName','t=50,ci 1e-2')
plot(sig,ci2_p1(:,150),'r-','DisplayName','t=150,ci 1e-2')
plot(sig,ci2_p1(:,501),'r.-','DisplayName','t=500,ci 1e-2')

%Cond inicial 1
plot(sig,ci7_p1(:,5),'b:','DisplayName','t=5,ci 1e-7')
plot(sig,ci7_p1(:,50),'b--','DisplayName','t=50,ci 1e-7')
plot(sig,ci7_p1(:,150),'b-','DisplayName','t=150,ci 1e-7')
plot(sig,ci7_p1(:,501),'b.-','DisplayName','t=500,ci 1e-7')

%plot(iy,ci2_p1_2000(:,2001),'k','DisplayName','t=2000,ci 1e-2')
%plot(iy,ci2_p1_2000(:,1001),'g','DisplayName','t=1000,ci 1e-2')

ax=gca;
ax.XAxis.Limits = [1e-6 1e-3];
ax.XAxis.Scale = 'log';
xlabel('IPTG')
ylabel('x mean')
title('Transient hysteresis in Protein 1')
legend

figure(2)
%Cond inicial 2
plot(sig,ci2_p2(:,5),'r:','DisplayName','t=5,ci 1e-2')
hold on
plot(sig,ci2_p2(:,50),'r--','DisplayName','t=50,ci 1e-2')
plot(sig,ci2_p2(:,150),'r-','DisplayName','t=150,ci 1e-2')
plot(sig,ci2_p2(:,501),'r.-','DisplayName','t=500,ci 1e-2')

%Cond inicial 1
plot(sig,ci7_p2(:,5),'b:','DisplayName','t=5,ci 1e-7')
plot(sig,ci7_p2(:,50),'b--','DisplayName','t=50,ci 1e-7')
plot(sig,ci7_p2(:,150),'b-','DisplayName','t=150,ci 1e-7')
plot(sig,ci7_p2(:,501),'b.-','DisplayName','t=500,ci 1e-7')


%plot(iy,ci2_p2_2000(:,2001),'k','DisplayName','t=2000,ci 1e-2')
%plot(iy,ci2_p2_2000(:,1001),'g','DisplayName','t=1000,ci 1e-2')

ax=gca;
ax.XAxis.Limits = [1e-6 1e-3];
ax.XAxis.Scale = 'log';
xlabel('IPTG')
ylabel('x mean')
title('Transient hysteresis in Protein 2')
legend

%% Plot en 3D
time=[25,50,100,150,200,300,400,500];

copia_p1=ci2_p1_i2_500(:,501);
copia_p1(18:20)=ci2_p1_2000(:,2001);
ci2_p1_i2000=[ci2_p1_i1_500(:,501);copia_p1;ci2_p1_i3_500(:,501)];
copia_p1(18:20)=ci2_p1_2000(:,1001);
ci2_p1_i1000=[ci2_p1_i1_500(:,501);copia_p1;ci2_p1_i3_500(:,501)];

copia_p2=ci2_p2_i2_500(:,501);
copia_p2(18:20)=ci2_p2_2000(:,2001);
ci2_p2_i2000=[ci2_p2_i1_500(:,501);copia_p2;ci2_p2_i3_500(:,501)];
copia_p2(18:20)=ci2_p2_2000(:,1001);
ci2_p2_i1000=[ci2_p2_i1_500(:,501);copia_p2;ci2_p2_i3_500(:,501)];


figure(3)
cero=ones(length(t),length(s));
% plot3(s',cero(1,:)'*0,ci2_p1(:,1),'r.','DisplayName','t=0')
% hold on
% plot3(s',cero(1,:)'*0,ci7_p1(:,1),'b.','DisplayName','t=0')
plot3(s',cero(1,:)'*5,ci2_p1(:,26),'r:','DisplayName','t=25')
hold on
plot3(s',cero(1,:)'*5,ci7_p1(:,26),'b:','DisplayName','t=25')
plot3(s',cero(1,:)'*10,ci2_p1(:,51),'r--','DisplayName','t=50')
plot3(s',cero(1,:)'*10,ci7_p1(:,51),'b--','DisplayName','t=50')
plot3(s',cero(1,:)'*15,ci2_p1(:,101),'r.-','DisplayName','t=100')
plot3(s',cero(1,:)'*15,ci7_p1(:,101),'b.-','DisplayName','t=100')
plot3(s',cero(1,:)'*20,ci2_p1(:,151),'r-','DisplayName','t=150')
plot3(s',cero(1,:)'*20,ci7_p1(:,151),'b-','DisplayName','t=150')
plot3(s',cero(1,:)'*25,ci2_p1(:,201),'r:','DisplayName','t=200')
plot3(s',cero(1,:)'*25,ci7_p1(:,201),'b:','DisplayName','t=200')
plot3(s',cero(1,:)'*30,ci2_p1(:,301),'r:','DisplayName','t=300')
plot3(s',cero(1,:)'*30,ci7_p1(:,301),'b:','DisplayName','t=300')
plot3(s',cero(1,:)'*35,ci2_p1(:,401),'r--','DisplayName','t=400')
plot3(s',cero(1,:)'*35,ci7_p1(:,401),'b--','DisplayName','t=400')
plot3(s',cero(1,:)'*40,ci2_p1(:,501),'r-.','DisplayName','t=500')
plot3(s',cero(1,:)'*40,ci7_p1(:,501),'b-.','DisplayName','t=500')
%plot3(s',cero(1,:)'*45,ci2_p1_i1000,'r-','DisplayName','t=1000')
%plot3(s',cero(1,:)'*45,ci7_p1(:,501),'b-.','DisplayName','t=1000')
%plot3(s',cero(1,:)'*50,ci2_p1_i2000,'r-','DisplayName','t=2000')
%plot3(s',cero(1,:)'*50,ci7_p1(:,501),'b-.','DisplayName','t=2000')
ax=gca;
ax.XAxis.Limits = [1e-6 1e-3];
ax.XAxis.Scale = 'log';
xlabel('IPTG')
ylabel('Time')
zlabel('x mean')
title('Transient hysteresis in Protein 1')


figure(4)
cero=ones(length(t),length(s));
% plot3(s',cero(1,:)'*0,ci2_p2(:,1),'r.','DisplayName','t=0')
% hold on
% plot3(s',cero(1,:)'*0,ci7_p2(:,1),'b.','DisplayName','t=0')
plot3(s',cero(1,:)'*5,ci2_p2(:,26),'r:','DisplayName','t=25')
hold on
plot3(s',cero(1,:)'*5,ci7_p2(:,26),'b:','DisplayName','t=25')
plot3(s',cero(1,:)'*10,ci2_p2(:,51),'r--','DisplayName','t=50')
plot3(s',cero(1,:)'*10,ci7_p2(:,51),'b--','DisplayName','t=50')
plot3(s',cero(1,:)'*15,ci2_p2(:,101),'r.-','DisplayName','t=100')
plot3(s',cero(1,:)'*15,ci7_p2(:,101),'b.-','DisplayName','t=100')
plot3(s',cero(1,:)'*20,ci2_p2(:,151),'r-','DisplayName','t=150')
plot3(s',cero(1,:)'*20,ci7_p2(:,151),'b-','DisplayName','t=150')
plot3(s',cero(1,:)'*25,ci2_p2(:,201),'r:','DisplayName','t=200')
plot3(s',cero(1,:)'*25,ci7_p2(:,201),'b:','DisplayName','t=200')
plot3(s',cero(1,:)'*30,ci2_p2(:,301),'r:','DisplayName','t=300')
plot3(s',cero(1,:)'*30,ci7_p2(:,301),'b:','DisplayName','t=300')
plot3(s',cero(1,:)'*35,ci2_p2(:,401),'r--','DisplayName','t=400')
plot3(s',cero(1,:)'*35,ci7_p2(:,401),'b--','DisplayName','t=400')
plot3(s',cero(1,:)'*40,ci2_p2(:,501),'r-.','DisplayName','t=500')
plot3(s',cero(1,:)'*40,ci7_p2(:,501),'b-.','DisplayName','t=500')
%plot3(s',cero(1,:)'*45,ci2_p2_i1000,'r-','DisplayName','t=1000')
%plot3(s',cero(1,:)'*45,ci7_p2(:,501),'b-.','DisplayName','t=1000')
%plot3(s',cero(1,:)'*50,ci2_p2_i2000,'r-','DisplayName','t=2000')
%plot3(s',cero(1,:)'*50,ci7_p2(:,501),'b-.','DisplayName','t=2000')
ax=gca;
ax.XAxis.Limits = [1e-6 1e-3];
ax.XAxis.Scale = 'log';
xlabel('IPTG')
ylabel('Time')
zlabel('x mean')
title('Transient hysteresis in Protein 2')

% %% Surfs
% 
% solt=linspace(0,500,501);
% 
% ci2_p1_r(:,1)=ci2_p1(:,1);
% ci2_p2_r(:,1)=ci2_p2(:,1);
% ci7_p1_r(:,1)=ci7_p1(:,1);
% ci7_p2_r(:,1)=ci7_p2(:,1);
% for i=1:(length(ci2_p1)/10)
%     ci2_p1_r(:,i+1)=ci2_p1(:,i+10);
%     ci2_p2_r(:,i+1)=ci2_p2(:,i+10);
%     ci7_p1_r(:,i+1)=ci7_p1(:,i+10);
%     ci7_p2_r(:,i+1)=ci7_p2(:,i+10);
% end
% 
% figure(5)
% surf(sig,solt,ci2_p1')
% ax=gca;
% ax.XAxis.Limits = [1e-6 1e-3];
% ax.XAxis.Scale = 'log';
% xlabel('IPTG')
% ylabel('Time')
% zlabel('x mean')
% colormap('jet')
% colorbar
% title('Protein 1 with initial condition IPTG=1e-2')
% 
% figure(6)
% surf(sig,solt,ci2_p2')
% ax=gca;
% ax.XAxis.Limits = [1e-6 1e-3];
% ax.XAxis.Scale = 'log';
% xlabel('IPTG')
% ylabel('Time')
% zlabel('x mean')
% colormap('jet')
% colorbar
% title('Protein 2 with initial condition IPTG=1e-2')
% % % 
% % 
% % % Gráfica de la evolución con el tiempo y la señal para la coincición
% % % inicial 2
% 
% figure(7)
% surf(sig,solt,ci7_p1')
% ax=gca;
% ax.XAxis.Limits = [1e-6 1e-3];
% ax.XAxis.Scale = 'log';
% xlabel('IPTG')
% ylabel('Time')
% zlabel('x mean')
% colormap('jet')
% colorbar
% title('Protein 1 with initial condition IPTG=1e-7')
% % 
% figure(8)
% surf(sig,solt,ci7_p2')
% ax=gca;
% ax.XAxis.Limits = [1e-6 1e-3];
% ax.XAxis.Scale = 'log';
% xlabel('IPTG')
% ylabel('Time')
% zlabel('x mean')
% colormap('jet')
% colorbar
% title('Protein 2 with initial condition IPTG=1e-7')
