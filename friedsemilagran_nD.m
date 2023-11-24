function [x,Xgrid,TT,PX_sol]=friedsemilagran_nD(path_forder_DR,dato,SLdato,IPTG)
%%%
% [x,Xgrid,TT,PX_sol]=friedsemilagran_nD(path_forder_DR,dato,SLdato)
% Numerical solution of the Friedman equation in general dimension with 
% the semilagrangian method
%%%

PathCurrent_SL = pwd;

% Dimensionless parameters (gamma1 = 0.01 and gamma2 = 4e-4)
b=cell(dato.n_gene,1);
for i=1:dato.n_gene
    b{i}=dato.R_constants(i,2)/dato.R_constants(i,3);
end      

% Spatial discretization
iN=cell(dato.n_gene,1);
x=cell(dato.n_gene,1);
for i=1:dato.n_gene
    iN{i} = SLdato.Prot_mesh(i,3) + 1;
    x{i} = linspace(SLdato.Prot_mesh(i,1),SLdato.Prot_mesh(i,2), iN{i});
%     fprintf('\n The discretization step in the dimension %d is: %g \n',i,x{i}(2)-x{i}(1));
end

% Protein "spatial" Mesh
Xgrid=cell(dato.n_gene,1);
[Xgrid{1:dato.n_gene}] = ndgrid(x{1:dato.n_gene});

% Time definition
t0     = SLdato.Time_mesh(1);
tmax   = SLdato.Time_mesh(2);
nt     = SLdato.Time_mesh(3)*SLdato.Time_mesh(4) + 1;
deltat = (tmax-t0)/nt;
tl     = linspace(t0, tmax, nt);
% fprintf('\n The time discretization is: %g \n',tl(2)-tl(1));

% Initial conditions (Gaussian density function)
cd(path_forder_DR);

%PX0un=IC_Function(x,Xgrid);

n_gene = length(Xgrid);

% % Parameters of the Gaussian density funcion
% x0g=[50 50];  % Mean of the Gaussian density
% sigmag=[5 5]; % Standard deviation of the Gaussian density 
% 
% gausker=0;
% for i=1:n_gene
%     gausker = gausker+((Xgrid{i}-x0g(i))/sigmag(i)).^2;
% end
% PX0un=exp(-gausker/2);

%Introducimos la matriz de la distribución de probabilidad del estacionario
load(fullfile(path_forder_DR,'distIPTG2.mat'))
load(fullfile(path_forder_DR,'distIPTG7.mat'))
load(fullfile(path_forder_DR,'distIPTG6.mat'))

PX0un=distIPTG2;

% Norm of the initial condition
auxnor0=PX0un;
for i=1:n_gene
    auxnor0 = trapz(x{i},auxnor0);
end

% Normalization of the initial condition to be a density function
IC=PX0un/auxnor0;

PX0un=IC;

cd(PathCurrent_SL);

% Normalized initial condition
auxnor0=PX0un;
for i=1:dato.n_gene
    auxnor0 = trapz(x{i},auxnor0);
end
% initial condition normalized
PX0=PX0un/auxnor0;

% Computation of characteristics curves
xbar=cell(dato.n_gene,1);
xbarlim=cell(dato.n_gene,1);
for i=1:dato.n_gene
    xbar{i}=x{i}*exp(deltat*dato.R_constants(i,4));
    xbarlim{i} = find(xbar{i}>=x{i}(end));
%     fprintf('\nThe length of X%g is: %g and there are %g points of Xbar biggest than Xmax \n',i,length(x{i}),length(xbarlim{i}));
end
% Caracteristics grid
Xbargrid=cell(dato.n_gene,1);
[Xbargrid{1:dato.n_gene}] = ndgrid(xbar{1:dato.n_gene});

% Initialization 
PX      = PX0;

% Time Independent functions
e_x=cell(dato.n_gene,1);
e_lx=cell(dato.n_gene,1);
for i=1:dato.n_gene
    e_x{i}=exp(Xgrid{i}/b{i});
    e_lx{i}=exp(-Xgrid{i}/b{i});
    e_x{i}(isinf(e_x{i})==1)=realmax;
    e_lx{i}(isinf(e_lx{i})==1)=realmax;
end

% Input function, c_i(x), construction:
% if strcmp(dato.IF_Type,'Hill')==1
%     cx=cell(dato.n_gene,1);
%     for i=1:dato.n_gene
%         cx{i}=IF_FeedbackMechanism(Xgrid,dato,i);
%     end
% else
%     cd(path_forder_DR);
%     cx=cell(dato.n_gene,1);
%     for i=1:dato.n_gene
%         cx{i}=IF_FM_user(Xgrid,i,IPTG);
%     end
%     cd(PathCurrent_SL);
% end
beta=2.5;
gamma=1;
eta=2.0015;
K=2.9618e-5;
%IPTG=40e-6;
phi=(1+IPTG/K).^eta;
eps1=0.1;
eps2=0.075;
cx_cell{1}=1./(1+(Xgrid{2}).^beta)+ (1-1./(1+(Xgrid{2}).^beta))*eps1 ;
cx_cell{2}=1./(1+((Xgrid{1})./phi).^gamma)+ (1-1./(1+((Xgrid{1})./phi).^gamma))*eps2;

for i=1:dato.n_gene
    cx{i}=cx_cell{i};
end


% Other time independent functions
sumkmcx=0;
for i=1:dato.n_gene
    sumkmcx = sumkmcx + dato.R_constants(i,1)*cx{i};
end
sumprotdeg=0;
for i=1:dato.n_gene
    sumprotdeg = sumprotdeg + dato.R_constants(i,4);
end
expl_den = 1+(sumkmcx - sumprotdeg)*deltat;


% Saving only dato.Time_mesh(4) simulated times 
nt_sol = SLdato.Time_mesh(4) + 1;
PX_sol      = cell(nt_sol,1);
TT          = zeros(nt_sol,1);
PX_sol{1}   = PX;
TT(1)       = tl(1);    
kk_sol = 1;
jjkk=(nt-1)/(nt_sol-1)+1:(nt-1)/(nt_sol-1):nt;


for j = 2:nt
    %tic
    % PX_bar construction using interpolation 
    if dato.n_gene==1
        PX_bar = interp1(Xgrid{1},PX,Xbargrid{1});
        PX_bar(isnan(PX_bar)==1)=0;
    elseif dato.n_gene>1
        PX_bar = interpn(Xgrid{1:dato.n_gene},PX,Xbargrid{1:dato.n_gene});
        PX_bar(isnan(PX_bar)==1)=0;
    end
       
     % Integral term computation by numerical integration
    Lix=0;
    for i=1:dato.n_gene
        Lix = Lix + dato.R_constants(i,1)/b{i}*e_lx{i}.*cumtrapz(x{i},e_x{i}.*cx{i}.*PX,i);
    end
       
    % Explicit method    
    PX = (PX_bar+deltat*Lix)./expl_den; 
    
    % Zero boundary condition
    CFaux=cell(dato.n_gene,1);
    for i=1:dato.n_gene
        CFaux{i}=':';
    end
    for i=1:dato.n_gene
        CF=CFaux;
        CF{i}=iN{i};
        PX(CF{1:dato.n_gene})=zeros(size(PX(CF{1:dato.n_gene})));
    end
    
    % Normalization: int_xmin^xmax(PX)dx=1
    auxnorpx=PX;
    for i=1:dato.n_gene
        auxnorpx = trapz(x{i},auxnorpx);
    end
    PX=PX/auxnorpx;
            
    % Saving the solution fo the current time step
    if j==jjkk(kk_sol)
%         fprintf('Time = %f \n',tl(j))
        PX_sol{kk_sol+1} = PX;
        TT(kk_sol+1)       = tl(j);
        kk_sol = kk_sol+1;
    end
    %toc
end


end



