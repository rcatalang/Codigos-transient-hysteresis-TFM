%Diagrama tal que:
% eje x: señal(var IPTG/b)
% eje y: valor de x_mean (sale de la integral de la marginal)
% eje x: tiempo (solution.T)
% 

clear x_mean1;
clear x_mean2;
%s1=logspace(-100,-8,15);
%s2=logspace(-7,-6,30);
%i1=logspace(-6,-5,20);
i2=logspace(-5,-4,20);
%ix=i2(17:20);
iy=i2(18:20);
%i3=logspace(-4,-3,20);
%s4=logspace(-3,0,40);
%s5=logspace(0,100,10);
sig=[iy];
%sig=logspace(-7,-2,10);
% sig=1e-6;

for j=1:length(sig)
    display(['Inicio simulacion:',num2str(sig(j))])
    tic
    sol(j)=SELANSI_Solve('Simulacion3',sig(j));
    solt(:,j)=sol(j).T;
    solx1(j,:)=sol(j).x{1};
    solx2(j,:)=sol(j).x{2};

    for i=1:length(sol(j).PTX)
        PX=sol(j).PTX{i};
        y1=sum(PX,2);
        %x_mean1(j,1)=sig(j);
        x_mean1(j,i)=trapz(solx1(j,:),solx1(j,:).*y1');
        y2=sum(PX,1);
        %x_mean2(j,1)=sig(j);
        x_mean2(j,i)=trapz(solx2(j,:),solx2(j,:).*y2);
    end
    toc
end


% ci2_p1_i1_500=x_mean1;
% ci2_p2_i1_500=x_mean2;
% save('ci2_p1_i1_500.mat','ci2_p1_i1_500');
% save('ci2_p2_i1_500.mat','ci2_p2_i1_500');




