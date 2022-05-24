function eps=getdata(x,materials,name,lab)

index = find(strcmp(materials, name));
ind=find(x{1,index}(:,1)==lab);

eps=(x{1,index}(ind,2)-1i*x{1,index}(ind,3))^2;

end