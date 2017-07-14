function PlotStationaryContrib(R,whichvar )
%presenting variance decomposition in the stationary state as a piechart

% extracting elements correspond to the whichvar variable
M=R{1}(whichvar,:); 
for i =2:length(R)    
   M=[M;R{i}(whichvar,:)];
end

pie((M(:,size(M,2))/sum(M(:,size(M,2)))))

%% legend cereating
%defining elements of legend
disp(length(R));

for i=1:length(R)
    leg{i}=['reaction ',int2str(i)];
end
%creating legend
legend(leg, 'Location','Best');

end


