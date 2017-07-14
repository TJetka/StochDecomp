function [Stoichiometry_matrix species_names] = createSBML(varargin)
% reads vector of initial values from SBML model in XML
% and generates symbolic equations and derivatives to files (optional)

% usage: fromSBMLToCluster(name, [optional, default=1] generateODE)

% the SBML model needs to be in /models/name/name.xml
% if generateODE=1 we need also functions 'createODE' and 'WriteODEFunctionWithParams'
% and installed SBMLToolbox-4.1.0 or higher
% if generateODE=0 the ODE files need to be in /models/name

switch(nargin)
    case 0
        error('fromSBMLToCluster(name, (optional) generateODE)\n%s', 'must have at least one argument');
    case 1
        name = varargin{1};
    otherwise
        error('fromSBMLToCluster(name, (optional) isOdeGenerated)\n%s', 'does not take more than one argument');
end;

disp('reading the model from XML...')
model = TranslateSBML(['models/', name, '/', name, '.xml']);
disp('analysing species and parameters...')
if (isfield(model, 'time_symbol'))
        if (~isempty(model.time_symbol))
            timeVariable = model.time_symbol;
        else
            timeVariable = 't';
        end;
    else
        timeVariable = 't';
    end;

[param_Names, param_Values] = GetAllParametersUnique(model);
param_num=length(param_Values);

species = AnalyseSpecies(model);
species_names=GetSpecies(model);
species_num=length(species);

%% Rates file
disp('generating rates file...')
fileName = fullfile(pwd,'models',name,[name,'_rates.m']);
reaction=Model_getListOfReactions(model);
reaction_Num=length(reaction);

fileID = fopen(fileName, 'w');
fprintf(fileID,  'function R = %s_rates(x, par, %s, stimulus)\n', name, timeVariable);
fprintf(fileID,  'R =[\n');
 for i=1:reaction_Num
    kineticLaw=Reaction_getKineticLaw(reaction(i));
    formula=KineticLaw_getMath(kineticLaw);
    formula=ChangeNames(formula,species_names,'species');
    formula=ChangeNames(formula,param_Names,'par');
    %Array{i} = sprintf('%s', );
    fprintf(fileID, '%s;\n',char(formula));
 end
 fprintf(fileID,  '];\n end');
fclose(fileID);

%% stoich file
disp('generating file containing stoichiometry matrix...')
[Stoichiometry_matrix species_names]=GetStoichiometryMatrix(model);
dlmwrite(fullfile(pwd,'models',name,[name,'_stoich.txt']),Stoichiometry_matrix,' ');

%% par file
disp('generating file containing parameter set...')
param_Values=param_Values';
fName=fullfile(pwd,'models',name,[name,'.par']);
fid = fopen(fName,'w');            %# Open the file
  for i=1:param_num
  fprintf(fid,'%s %f \n',param_Names{i}, param_Values(i));       %# Print the string
  end
fclose(fid);                     %# Close the file

%% species file
disp('generating file containing species order...')
fName=fullfile(pwd,'models',name,[name,'_species.txt']);
fid = fopen(fName,'w');            %# Open the file
  for i=1:species_num
  fprintf(fid,'%s x(%u) \n',species_names{i}, i);       %# Print the string
  end
fclose(fid);                     %# Close the file

end

