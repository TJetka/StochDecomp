close all;
clear all;

% name of the model
name='MICHAELISMENTEN';

% symbolic computations to create appropiate equations stored in the
% SYMBOLIC  folder
%createStochDecomp(name);  %run only once

addpath([pwd, '/models','/',name]); % adding the path to the model's folder
addpath([pwd,'/models/','/',name,'/symbolic/']); % adding the path to the folders containing results of symbolic computations
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); % reading in parameters


% initial condition
% first 4 elements of the vector y0 are mean initial molecular numnbers,
% elements 5-8 are initial variances in molecular numnbers, 
% remaining elemens are initial covariances of molecular numnbers stored as
% above-diagonal concatenation of the initial covariance matrix

y0(1)=3;
y0(2)=50;
y0(3)=0;
y0(4)=3;

y0(5)=0;
y0(6)=0;
y0(7)=0;
y0(8)=0;

y0(9)=0;
y0(10)=0;
y0(11)=0;

y0(12)=0;
y0(13)=0;
y0(14)=0;

par(1)=0.1;
par(2)=50;
par(3)=1;
par(4)=20;
par(5)=0.1;

times=[0.1,2,3, 130];  %time at which contributions are calculated (large enough to ensure stationarity)


% legend descibing individual reactions
leg{1}='subtrate +';
leg{2}='complex forward';
leg{3}='complex backward';
leg{4}='product +';
leg{5}='prod. degradation';

R=CalcContrib(name,times,y0,par);  % calculating contributions of all reactions into variability of all variables

whichvar=1;
PlotStationaryContrib(R,whichvar); %plotting contributions to the 1st variable
legend(leg);
title('Substrate')

figure;
whichvar=2;
PlotStationaryContrib(R,whichvar); %plotting contributions to the 2nd variable
legend(leg);
title('Enzyme')

figure;
whichvar=3;
PlotStationaryContrib(R,whichvar); %plotting contributions to the 3rd variable
legend(leg);
title('Compex')

figure;
whichvar=4;
PlotStationaryContrib(R,whichvar); %plotting contributions to the 4th variable
legend(leg);
title('Product')