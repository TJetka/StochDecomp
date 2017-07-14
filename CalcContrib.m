function R = CalcContrib(name,times,y0,parv)
% function numerically calculates variance contributions resulting from
% all of the model reactions based on symbolic calculations performed by the function
% createStochDecomp()

addpath([pwd, '/models','/',name,'/symbolic/' ]); %adding path 
stoichiometry=load([pwd, '/models/','/',name,'/' ,name,'_stoich.txt']);  %reading stoichiometry matrix
nreac=size(stoichiometry); % number of reactions
nreac=nreac(2);


parfor i=1:nreac, % calculate contribution of each reaction
    R{i}=CalcContrib_i(name,i,times,y0,parv);
end

end
