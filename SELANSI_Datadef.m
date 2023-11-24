% Reactions (i=1,...,n, with n=number of genes involved)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %       0 - c_i(x)km_i -> mRNA_i        %
        %  mRNA_i -    kx_i    -> mRNA_i + X_i  %
        %  mRNA_i - gamma_m_i  -> 0             %
        %     X_i -gamma_x_i(x)-> 0             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_i(x)       := input function which collects the feedback mechanism
% km_i         := mRNA_i production rate
% kx_i         := X_i production rate
% gamma_m_i    := mRNA_i degradation rate
% gamma_x_i(x) := X_i degradation rate (can be a variable function)
close all
clear all

% Number of genes considered

n_gene=2;

% Degradation functions are considered contants
DF_Type=0;

% Obtain the actual path
PathCurrent = pwd;

% Create a forder to save data and results

Folder_name='Simulacion3';
path_forder_DR=fullfile(PathCurrent,'DATA',Folder_name);
mkdir(path_forder_DR);
mkdir(fullfile(path_forder_DR,'Reaction_data'));
mkdir(fullfile(path_forder_DR,'Mesh_data'));
mkdir(fullfile(path_forder_DR,'Results'));

%Crear carpeta auxiliar para almacenar el resultado final


% Feedback mechanism constants cx=IF_FeedbackMechanism(IF_Type,n,x,H,K,epsilon)
% IFT_define=input('Choose your feedback mechanism input function \n Insert 1 if you want to use a predefined Hill function \n Insert 0 if you define your input function \n');
% if IFT_define==1
%     IF_Type = 'Hill'; % this option can be variable in the future
% else
%     IF_Type = 'user';
% end
IF_Type='user';
IFT_define=0;
copyfile(fullfile(pwd,'SELANSI_Files','InputFunction','IF_FM_user.m'),path_forder_DR);



% Matrix (nx4) of input variables (reaction contants)
R_constants=ones(n_gene,4);
R_constants(1,4)=1;
R_constants(1,3)=50;
R_constants(1,2)=450;
R_constants(1,1)=7.812500;
R_constants(2,4)=1;
R_constants(2,3)=10;
R_constants(2,2)=50;
R_constants(2,1)=3.120000;
FID=fopen(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'),'w+');
for i=1:n_gene
   fprintf(FID,'%f %f %f %f \n',R_constants(i,:));
end
fclose(FID);

% Mesh for the semilagrangian method
% Number of proteins mesh, [Xi_min, Xi_max, deltaxi]
Prot_mesh=zeros(n_gene,3);
Prot_mesh(1,2)=150;
Prot_mesh(1,3)=750;
Prot_mesh(2,2)=50;
Prot_mesh(2,3)=250;

FID=fopen(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'),'w+');
for i=1:n_gene
    fprintf(FID,'%f %f %f \n',Prot_mesh(i,:));
end
fclose(FID);

% Mesh of the semilagrangian method
Time_mesh=[0 50 100 50];
FID=fopen(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'),'w+');
fprintf(FID,'%f %f %f %f \n',Time_mesh);
fclose(FID);

% Initial condition for the semiLagrangian method
copyfile(fullfile(pwd,'SELANSI_Files','InitialCondition','IC_Function.m'),path_forder_DR)



% Loading the data to be saved in the variables
dato.R_constants=load(fullfile(path_forder_DR,'Reaction_data','R_constants.txt'));
dato.n_gene=n_gene;
dato.DF_Type=DF_Type;
dato.IF_Type=IF_Type;
SLdato.Prot_mesh=load(fullfile(path_forder_DR,'Mesh_data','Prot_mesh.txt'));
SLdato.Time_mesh=load(fullfile(path_forder_DR,'Mesh_data','Time_mesh.txt'));
save(fullfile(path_forder_DR,'Reaction_data','parameters.mat'),'dato');
save(fullfile(path_forder_DR,'Mesh_data','SL_parameters.mat'),'SLdato');
% end

%fprintf('\n Your parameters are saved in %s \n',path_forder_DR)

