close all;
clear all;

% name of the model
name='CONVERSION';

% symbolic computations to create appropiate equations stored in the
% SYMBOLIC  folder
% createStochDecomp(name);  %run only once

addpath([pwd, '/models','/',name]); % adding the path to the model's folder
addpath([pwd,'/models/','/',name,'/symbolic/']); % adding the path to the folders containing results of symbolic computations
[parn, par, parnames] = textread([pwd, '/models/','/',name,'/'  ,name, '.par'], '%s %f %q'); % reading in parameters


% initial condition
% first 3 elements of the vector y0 are mean initial molecular numnbers,
% elements 4-6 are initial variances in molecular numnbers, 
% remaining elemens are initial covariances of molecular numnbers stored as
% above-diagonal concatenation of the initial covariance matrix

y0(1)=3;
y0(2)=3;
y0(3)=3;

y0(4)=0;
y0(5)=0;
y0(6)=0;

y0(7)=0;
y0(8)=0;
y0(9)=0;

times=[1,2,30]; %time at which contributions are calculated (large enough to ensure stationarity)


parslow=par; % parameters for slow reactions regime


R=CalcContrib(name,times,y0,parslow);  % calculating contributions of all reactions into variability of all variables

whichvar=3;
PlotStationaryContrib(R,whichvar);  %plotting contributions to the 3rd variable

% parameters for fast reactions regime
parfast=parslow;
parfast(2)=10;
parfast(3)=10;

R=CalcContrib(name,times,y0,parfast);  % calculating contributions of all reactions into variability of all variables
whichvar=3;
figure;
PlotStationaryContrib(R,whichvar);
