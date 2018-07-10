function pts =densify(p1,spacing)
np=length(p1);
p2=zeros(np,2); padd=zeros(floor(np/3),2); 
k=0;
for i = 1:np-1
    if isnan(p1(i,:)) 
        k=k+1;
        p2(k,1:2)=NaN; 
        continue;
    end
    if isnan(p1(i+1,:))
        k=k+1; 
        p2(k,1:2)=NaN;
        continue
    end; 
    d = sqrt((p1(i,1)-p1(i+1,1)).^2 + (p1(i,2)-p1(i+1,2)).^2);
    nsplit = ceil(d/spacing);
    k=k+1;
    p2(k,:)=p1(i,:); 
    padd(1:nsplit,1)=(linspace(p1(i,1),p1(i+1,1),nsplit))';
    padd(1:nsplit,2)=(linspace(p1(i,2),p1(i+1,2),nsplit))';
    k=k+nsplit; 
    p2(k-nsplit+1:k,:)=padd(1:nsplit,:);
    k=k+1; 
    p2(k,:)=p1(i+1,:);
    padd=padd*0; 
end
x=p2(1:k,1); y=p2(1:k,2);
pts = [x,y];
end