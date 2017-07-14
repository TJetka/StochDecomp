function traj = CalcContrib_i(name,i,times,y0,par)
% function numerically calculates variance contributions resulting from
% i-th reaction based on symbolic calculations performed by the function
% createStochDecomp()
addpath([pwd, '/models','/',name,'/symbolic/' ]); %adding path 
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']); %reading stoichiometry matrix

npar=length(par); %number of paramters 
nvar=size(stoichiometry); %number of variables
nvar=nvar(1);
tspan=[0,max(times)]; %time interval of interest

all_equations=str2func([name,'_all',int2str(i),'_equations']); % name of equation to be solved
mysolution=ode15s(all_equations,tspan,y0,[],par);  % solving equation
traj=deval(mysolution,times); % solution
traj=traj( (nvar+1):(2*nvar),:); % extracting variances of interest from the solution
end
