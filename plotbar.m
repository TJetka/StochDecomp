function plotbar( R,whichvar )

M=R{1}(whichvar,:);

for i =2:length(R)    
   M=[M;R{i}(whichvar,:)];
end

bar(M','stack')

end


