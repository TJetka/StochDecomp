function savefunction(rhs,name)
%
% creates the file 'name'_jac.m with the model's jacobian
%
% rhs is the right hand side of the model



%namer = [mdir1,name,'_jac.m'];
% if ~CheckFileExists(namer)
%     return;
% end

sizef = size(rhs);

dim1=sizef(1);

dim2=sizef(2);

R=rhs


file = fopen(name,'w');

% writing the jacobian
fprintf(file,'%s\n\n','function R = f(t,y,p)');
%add the force options
%fprintf(file, '%%Model jacobian dy/dy\n\n');
fprintf(file,'%s\n','R = [');
for i=1:dim1  
    fprintf(file,'     [');
    for j=1:dim2
        fprintf(file,' (%s) ',char(R(i,j)));
    end
    fprintf(file,'];\n');
end   
fprintf(file,'];\n');
fclose(file);

