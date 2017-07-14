function [ output_args ] = createStochDecomp(name,scalevar,scalepar)
% function symbolically generates equations describing variance contributions
% resulting from individual reactions

%% reading the model definition files 
disp('reading the model definition files') 
DefsDir = [pwd, '/models/', name];  %where model definition files are found
[parn, parv, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); %reading information about parameters
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']); %reading stoichiometry matrix

%% creating folder to write
addpath([pwd, '/models','/',name]); % adding paths to the model's folder
mkdir([pwd,'/models/','/',name,'/symbolic/']); %creating a folder to write results of symbolic computations 
addpath([pwd,'/models/','/',name,'/symbolic/']); % adding paths to results of symbolic computations

size_sto=size(stoichiometry); %dimension of stoichiomety matrix
varn=size_sto(1); %number of model variebles
rean=size_sto(2); %number of model reactions

%%symbolic definition of stoichiometry matrix

syms S; % creating symbolic stoichiometry matrix
for i=1:size_sto(1),
    for j=1:size_sto(2),
        S(i,j)=stoichiometry(i,j); % filling symbolic stoichiometry matrix
    end
end


disp('creating definition of the model');
y=sym('y','positive'); % symbolic definition of the state vector


% symbolic definition of the mean and variances variables
for i=1:varn,   
    y(i) = sym((['y',num2str(i)]),'positive'); %symbolic variable 
    var_st{i}=['y(',num2str(i),')']; %name using brackets
    var_st_sym{i}=['y',num2str(i)];  %name without brackets
end


% symbolic definition of the stimulus
stim=sym('stim','positive' );
stm_st{1}=[name,'_stimulus(t)'];  %name using brackets
stm_st_sym{1}='stim'; %name without brackets


%symbolic definition of the parameter vector
param=sym('param','positive' );

% symbolic definition of time
t=sym('t','positive' );

pnum = length(parn);
for i = 1:pnum,
    param(i) = sym(['p',num2str(i)],'positive');  %symbolic variable 
    par_st{i}=['p(',num2str(i),')']; %name using brackets
    par_st_sym{i}=['p',num2str(i)];  %name without brackets
end

if nargin > 1 %if scaled variables are present
    for i = 1:length(scalepar),
      S(scalevar(i),:)=param(scalepar(i))*S(scalevar(i),:); %correcting stoichiometry matrix to include scaled variables    
    end
end
rates = str2func([name,'_rates']);  %name of the function contiaing reaction rates

%% creating  symbolic macroscopic rate equation (MRE)
disp('creating deterministic equations')
symrates=feval(rates,y,param,t,stim);
MRE=S*symrates;

%% creating symbolic Jacobian of MRE
disp('creating A matrix')
 
J=jacobian(MRE,y(1:varn));
J_pom=substitution(J,var_st_sym,var_st);
J_pom=substitution(J_pom,stm_st_sym,stm_st);
savefunction(substitution(J_pom,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,'_MRE_jacobian.m']);


%% creating symbolic variance equations
disp('creating variance equations')
E=sym('E','positive');
for j=1:rean,
    E(j,j)=(sqrt(symrates(j)));        
end
EE=sym('EE','positive');

for j=1:rean,
    EE(j,j)=((symrates(j)));        
end

D=sym('D','real'); % Creating symbolic matrix D
D=(S*EE)*(S');
Sigma=sym('Sigma','real');

totdim=2*varn+varn*(varn-1)/2; % calculating dimension of concatenation of deterministic state variables and variances

y(totdim)=0;

% giving names to symbolic variables
k=1;
for i=1:varn,
    for j=(i+1):varn,
        y(varn+varn+k) = sym(['y',num2str(varn+varn+k)],'real');
        var_st{varn+varn+k}=['y(',num2str(varn+varn+k),')'];
        var_st_sym{varn+varn+k}=['y',num2str(varn+varn+k)];
        k=k+1;   
    end
end


for i=1:varn,
    y(varn+i) = sym(['y',num2str(varn+i)],'real');
    var_st{varn+i}=['y(',num2str(varn+i),')'];
    var_st_sym{varn+i}=['y',num2str(varn+i)];
end

% assigning names to symbolic matrix Sigma
for i=1:varn,%diagonal elements
    Sigma(i,i)=y(varn+i);
end


%%% lower diagonal
k=1;
for i=1:varn,
    for j=(i+1):varn,
        Sigma(i,j)=y(varn+varn+k);
        k=k+1;   
    end
end
%using symmetry
Sigma=Sigma+Sigma';

%%%diagonal
for i=1:varn,
   Sigma(i,i)=y(varn+i);
end


Sigma_dot=J*Sigma+ Sigma*(J')+D; % Sigma equation


% rewriting variance equations from a matrix format to a vector format
sym('variances_dot','real');
% extracting variance equations into a cncatenated form
k=1;
for i=1:varn,
    variances_dot(i)=Sigma_dot(i,i);
    k=k+1;
end


for i=1:varn,
    for j=(i+1):varn,
        variances_dot(k)=Sigma_dot(i,j);
        k=k+1;   
    end
end

% create joint set of equations (deterministic + variances)
all_eq=[MRE; variances_dot'];
all_eq_pom=substitution(all_eq,var_st_sym,var_st);
all_eq_pom=substitution(all_eq_pom,stm_st_sym,stm_st);
savefunction(substitution(all_eq_pom,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,'_all_equations.m']);

% creating equations describing individual contributions
for zz=1:size_sto(2), % for each reaction
    syms Si;
    for i=1:size_sto(1),
        for j=1:size_sto(2),        
            if(j==zz)    
            Si(i,j)=S(i,j); %rewriting original stoichiometry matrix in jth column
            else   
                Si(i,j)=0; %setting remaining elemens to zero
            end
        end
    end
        % Symbolic Di matrix
        Di=sym('D1','real');
        Di=(Si*E)*((Si*E)');


        sym('variances_doti','real');
        Sigma_doti=J*Sigma+ Sigma*(J')+Di; %equation describing ith each contributions


        % extracting contribution equations into a cncatenated form
        k=1;
        for i=1:varn,
            variances_doti(i)=Sigma_doti(i,i);
            k=k+1;
        end


        for i=1:varn,
            for j=(i+1):varn,
                variances_doti(k)=Sigma_doti(i,j);
                k=k+1;   
            end
        end


        eq_name=['_all',int2str(zz),'_equations.m'];
        % create joint set of equations (deterministic + ith contributions)
        all_eqi=[MRE; variances_doti'];
        all_eq_pomi=substitution(all_eqi,var_st_sym,var_st);
        all_eq_pomi=substitution(all_eq_pomi,stm_st_sym,stm_st);
                
        savefunction(substitution(all_eq_pomi,par_st_sym,par_st),[pwd, '/models/','/',name,'/symbolic/',name,eq_name]); 
    end

end