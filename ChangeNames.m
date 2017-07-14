function out = ChangeNames(in,model_Names,mode)
%CHANGESPECIENAMES Summary of this function goes here
%Detailed explanation goes here
n=length(model_Names);
out=in;

if strcmp(mode,'species')
    
for i=n:-1:1
  out = strrep(out,model_Names{i},['x(',num2str(i),')']) ;
end

elseif strcmp(mode,'par')

for i=n:-1:1
  out = strrep(out,model_Names{i},['par(',num2str(i),')']) ;
end

else error('incorrect input')
end

