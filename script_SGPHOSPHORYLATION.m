clear all;
close all;

% name of the model
name='SGPHOSPHORYLATION';

% symbolic computations to create appropiate equations stored in the
% SYMBOLIC  folder
createStochDecomp(name); %run only once



addpath([pwd, '/models','/',name]); % adding the path to the model's folder
addpath([pwd,'/models/','/',name,'/symbolic/']); % adding the path to the folders containing results of symbolic computations

[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); % reading in parameters

times=20; %times vector 

% initial condition
% first 3 elements of the vector y0 are mean initial molecular numnbers,
% elements 4-6 are initial variances in molecular numnbers, 
% remaining elemens are initial covariances of molecular numnbers stored as
% above-diagonal concatenation of the initial covariance matrix

y0=zeros(1,9);


y0(1)=0;
y0(2)=0;
y0(3)=0;

y0(4)=0;
y0(5)=0;
y0(6)=0;

y0(7)=0;
y0(8)=0;
y0(9)=0;


par=[40, 2,2,1,0.5,1,1]; % redefining parameter values

R=CalcContrib(name,times,y0,par); %Calculation of the contributions of all reaction to all species at times give by the vector times

PlotStationaryContrib(R,3); %Plotting contribution of all reactions to the third variable.


% Calculating and plotting time dependent contributions. 
times=[0.0001, 1, 2, 5, 10, 100, 200];
R=CalcContrib(name,times,y0,par);
PlotTimeContrib(R,3);






% parameters for fast phosphorylation regime 
par(1)=40;
par(2)=2;
par(3)=2;
par(4)=1;
par(5)=0.5;
par(6)=1.0;
par(7)=1;

R=CalcContrib(name,times,y0,par); % calculating contributions of all reactions into variability of all variables
whichvar=3; 
figure;
PlotTimeContrib(R,whichvar)
figure;
PlotStationaryContrib(R,whichvar)  %plotting contributions to the 3rd variable

% parameters for slow phosphorylation regime 
par(1)=40;
par(2)=2;
par(3)=2;
par(4)=1;
par(5)=0.5;
par(6)=0.1;
par(7)=1;

R=CalcContrib(name,times,y0,par); % calculating contributions of all reactions into variability of all variables

whichvar=3;
figure;
PlotTimeContrib(R,whichvar)
figure;
PlotStationaryContrib(R,whichvar)  %plotting contributions to the 3rd variable


% parameters for "no phosphorylation" regime 
par(1)=40;
par(2)=2;
par(3)=2;
par(4)=1;
par(5)=0.5;
par(6)=0.0;
par(7)=1;

whichvar=3; 
figure;
PlotStationaryContrib(R,whichvar)  %plotting contributions to the 3rd variable

