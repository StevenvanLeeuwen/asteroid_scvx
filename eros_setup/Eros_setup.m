function [vert,facet]=Eros_setup(vf,x,y,z) %codegen

l=1;
m=1;
vert=zeros(856,3);
facet=zeros(1708,3);
for i=1:length(vf)
    if vf(i)=='v'
        vert(l,:)=[x(i) y(i) z(i)];
        l=l+1;
    end
    if vf(i)=='f'
        facet(m,:)=[x(i)+1 y(i)+1 z(i)+1];
        m=m+1;
    end
end
return