n=10;
p=0.5;
nsort=10000;
for s=1:nsort
    vector(s,1)=sum(rand(n,1)>p);
end
edges=0:10;
bar(edges,histc(vector,edges)/nsort)