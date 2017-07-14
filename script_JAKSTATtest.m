close all;
clear all;

% model name
name='JAKSTAT';
max(cell2mat(cellfun(@(x) str2num(x{1}),regexp(regexp(ratestxt,'stimulus{[0-9]+}','match'),'[0-9]+','match'), 'UniformOutput', false)))

% Creating StochDecomp readable files from a SBML file
createSBML(name); % creates StochDecomp-readable model definition from a

%SBML format does not allow to encode a time dependent stimulus
%If modelling involves a time dependent stimulation, as is the case in the JAK-STAT example, we suggest to repleace an auxiliary parameter with a time dependent stimulus.
% Here parmater 5 is replaced. One must execute lines below after running createSBML(name);
% 
% par_sub='par(5)';
% 
% fin = fopen([pwd, '/models','/',name,'/',name,'_rates.m']);
% fout = fopen([pwd, '/models','/',name,'/',name,'_ratesnew.m'],'w');
% 
% while ~feof(fin)
%    s = fgetl(fin);
%    s = strrep(s, par_sub, 'stimulus');
%    fprintf(fout,'%s',s);
%    fprintf(fout,'\n');
% end
% 
% fclose(fin)
% fclose(fout)
% delete([pwd, '/models','/',name,'/',name,'_rates.m']);
% movefile([pwd, '/models','/',name,'/',name,'_ratesnew.m'],[pwd, '/models','/',name,'/',name,'_rates.m'])



%createStochDecomp(name,[14,15], [6,7]); % symbolic computations to create appropiate equations stored in the
% SYMBOLIC  folder, run only once - can take several minutes


addpath([pwd, '/models','/',name]); % adding the path to the model's folder
addpath([pwd,'/models/','/',name,'/symbolic/']); % adding the path to the folders containing results of symbolic computations
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); % reading in parameters


% converting model paramars into molecular numbers units
W_C=6.02*(10^23)*1400*10^(-18)*10^(-6); % conversion constant (from mM to molecular numbers)
par(2)=par(2)/W_C;
par(7)=par(7)/W_C;

% defining  times to calculate contributions
times=50+[0,2,4,6,8,10,15,20, 30,40,50, 60]; % initial constant is set to 50 in order to allow the system to converge to a quasi-stationary state 


% defining initial conditions
% initial condition
% first 17 elements of the vector y0 are mean initial molecular numnbers,
% elements 18-34 are initial variances in molecular numnbers, 
% remaining elemens are initial covariances of molecular numnbers stored as
% above-diagonal concatenation of the initial covariances matrix

y0=zeros(1,170);
y0(1)=1.99989483593865340000e+002*W_C; % litrature value from Swameye et al. 2003
y0(15)=1.99989483593865340000e+002*W_C*par(7); % litrature value from Swameye et al. 2003

par(8)=0; % setting extrinsic noise variance to zero 
R=CalcContrib(name,times,y0,par);  % calculating contributions of all reactions into variability of all variables

kk=16
figure;
bar(([R{1}(kk,:);R{2}(kk,:);R{3}(kk,:);R{4}(kk,:)+R{5}(kk,:)+R{6}(kk,:)+R{7}(kk,:)+R{8}(kk,:)+R{9}(kk,:)+R{10}(kk,:)+R{11}(kk,:)+R{12}(kk,:);R{13}(kk,:);R{14}(kk,:)+R{15}(kk,:)])','stack')
% Alternatively
%PlotTimeContrib(R,16); % plotting contributions of ALL reactions to the 16th variable at all time points 

print -painters -dpdf -r600 NoExt.pdf %saving the plot 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Extrinsic noise analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extrinsic noise half life 0.01 min

R={}; %initialisation of a list to store computations

%initialising extrinsing noise parameters
ext_variance=(10/100)^2; % extrinsic noise variance 
T_12=0.01; % extrinsic noise half-life
par(9)=log(2)/T_12; % extrinsic noise death
par(8)=ext_variance*par(9); % extrinsic noise birth

%intial values for variables used to model extrinsic noise 
y0(17)=par(8)/par(9); %initial mean
y0(34)=par(8)/par(9); %initial variance

R=CalcContrib(name,times,y0,par);  % calculating contributions of all reactions into variability of all variables
kk=16; % variable index, x_16 denotes the variable describingtotal stat in the nucleus (sum of x_4 to x_13)

% plotting contributions of 1st,2nd and 3rd reaction; joint contribution of
% 4th, 5th, 6th, 7th,  8th, 9th, 10th, 11th, 12th, 13th; and joint contribution
% of 14th and 15th (extrinsic noise)
figure;
bar(([R{1}(kk,:);R{2}(kk,:);R{3}(kk,:);R{4}(kk,:)+R{5}(kk,:)+R{6}(kk,:)+R{7}(kk,:)+R{8}(kk,:)+R{9}(kk,:)+R{10}(kk,:)+R{11}(kk,:)+R{12}(kk,:);R{13}(kk,:);R{14}(kk,:)+R{15}(kk,:)])','stack')
print -painters -dpdf -r600 T001.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extrinsic noise half life 0.1 min T_12=0.1

R={};
% re-initialising extrinsing noise parameters to set noise half life 0.1 minutes
ext_variance=(10/100)^2;
T_12=0.1;
par(9)=log(2)/T_12;
par(8)=ext_variance*par(9);

%intial values for veriables used to model extrinsic noise 
y0(17)=par(8)/par(9); %initial mean
y0(34)=par(8)/par(9); %initial variance

R=CalcContrib(name,times,y0,par);  % calculating contributions of all reactions into variability of all variables


kk=16; % variable index, x_16 denotes the variable describingtotal stat in the nucleus (sum of x_4 to x_13)

% plotting contributions of 1st,2nd and 3rd reaction; joint contribution of
% 4th, 5th, 6th, 7th,  8th, 9th, 10th, 11th, 12th, 13th; and joint contribution
% of 14th and 15th (extrinsic noise)
figure
bar(([R{1}(kk,:);R{2}(kk,:);R{3}(kk,:);R{4}(kk,:)+R{5}(kk,:)+R{6}(kk,:)+R{7}(kk,:)+R{8}(kk,:)+R{9}(kk,:)+R{10}(kk,:)+R{11}(kk,:)+R{12}(kk,:);R{13}(kk,:);R{14}(kk,:)+R{15}(kk,:)])','stack')
print -painters -dpdf -r600 T01.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extrinsic noise half life 1 minute T_12=1

R={};
% re-initialising extrinsing noise parameters to set noise half life 1 minute
ext_variance=(10/100)^2;
T_12=1;
par(9)=log(2)/T_12;
par(8)=ext_variance*par(9);

%intial values for variables used to model extrinsic noise 
y0(17)=par(8)/par(9); %initial mean
y0(34)=par(8)/par(9); %initial variance

R=CalcContrib(name,times,y0,par);  % calculating contributions of all reactions into variability of all variables

kk=16; % variable index, x_16 denotes the variable describingtotal stat in the nucleus (sum of x_4 to x_13)
% plotting contributions of 1st,2nd and 3rd reaction; joint contribution of
% 4th, 5th, 6th, 7th,  8th, 9th, 10th, 11th, 12th, 13th; and joint contribution
% of 14th and 15th (extrinsic noise)
figure;
bar(([R{1}(kk,:);R{2}(kk,:);R{3}(kk,:);R{4}(kk,:)+R{5}(kk,:)+R{6}(kk,:)+R{7}(kk,:)+R{8}(kk,:)+R{9}(kk,:)+R{10}(kk,:)+R{11}(kk,:)+R{12}(kk,:);R{13}(kk,:);R{14}(kk,:)+R{15}(kk,:)])','stack')
print -painters -dpdf -r600 T1.pdf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Generating the posterior distributions of individual contributions
%%%%%% from experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('posterior.txt'); %reading in a sample from the posterior ditribution ( generated using PUA package : http://bmi.bmt.tue.nl/sysbio/software/pua.html )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lines 164 - 192 take several hours to complete 
%
% The result of the computations has been precalculated and saved as 'DENS.txt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% numbparticles=size(posterior,1); % number of particles to be sampled from the posterior 
% calculating contributions for each parameter values in the posterior
% R={}
% for(i=1:999)  
% y0loop=y0; %initial condition
% parchain=posterior(i,:);
% y0loop(15)=1.99989483593865340000e+002*Conversion*parchain(7); % scalled observable
% parchain(2)=parchain(2)/Conversion;
% parchain(6)=parchain(6)/Conversion;
% parchain(7)=parchain(7)/Conversion;
% R{i}=CalcContrib(name,times,y0loop,[parchain,0,1]) % converting parameters
% end
% numbreactions=15; %number of reactions in the model
% variable=16; % variable of interest (nuclear STAT molecules)
% 
% 
% timepoint=5; % timepoint to plot analysis for
% % extracting values of each of the contributions for each value of the
% % paramters in the posterior for the variable 16 and timepoint 5
% 
% DENS=zeros(numbreactions,numbparticles);
% for(k=1:numbreactions) %for every reaction     
%         hmatrix=zeros(timepoint,numbparticles);
%         for(i=1:numbparticles) % every particle
%           hvar=R{i};
%           hmatrix(i)=hvar{k}(variable,timepoint);       
%         end        
%         DENS(k,:)=hmatrix; %everage over particles       
% end

load('DENS.txt'); % reading in the precaluated results coded above

% Plotting posterior distributions of the total variance as well as of the individual contributions as in the Figure 2C of SI

figure;
hist(DENS(1,:)+DENS(2,:)+DENS(3,:)+DENS(4,:)+DENS(5,:)+DENS(6,:)+DENS(7,:)+DENS(8,:)+DENS(9,:)+DENS(10,:)+DENS(11,:)+DENS(12,:)+DENS(13,:))
figure;
ksdensity(DENS(1,:)+DENS(2,:)+DENS(3,:)+DENS(4,:)+DENS(5,:)+DENS(6,:)+DENS(7,:)+DENS(8,:)+DENS(9,:)+DENS(10,:)+DENS(11,:)+DENS(12,:)+DENS(13,:))
legend('total variance')

figure;
hist(DENS(1,:))
legend('Reaction 1')
figure;
ksdensity(DENS(1,:))
legend('Reaction 1')

figure;
hist(DENS(2,:))
legend('Reaction 2')
figure;
ksdensity(DENS(2,:))
legend('Reaction 2')

figure;
hist(DENS(3,:))
legend('Reaction 3')
figure;
ksdensity(DENS(3,:))
legend('Reaction 3')

figure;
hist(DENS(4,:)+DENS(5,:)+DENS(6,:)+DENS(7,:)+DENS(8,:)+DENS(9,:)+DENS(10,:)+DENS(11,:)+DENS(12,:))
legend('Reaction 4-12')
figure;
ksdensity(DENS(4,:)+DENS(5,:)+DENS(6,:)+DENS(7,:)+DENS(8,:)+DENS(9,:)+DENS(10,:)+DENS(11,:)+DENS(12,:))
legend('Reaction 4-12')


figure;
hist(DENS(13,:))
legend('Reaction 13')
figure;
ksdensity(DENS(13,:))
legend('Reaction 13')