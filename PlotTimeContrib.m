function PlotTimeContrib( R,whichvar )
%presenting variance decomposition in the out-of-steady-state (time
%dependent) decomposition as a barplot

% extracting elements correspond to th e whichvar variable
M=R{1}(whichvar,:);

for i =2:length(R)    
   M=[M;R{i}(whichvar,:)];
end

bar(M','stack')

%% legend cereating
%defining elements of legend
for i=1:length(R)
    leg{i}=['reaction ',int2str(i)];
end
%creating legend
legend(leg, 'Location','Best');

end


