<<<<<<< HEAD
n=10;
p=0.5;
nsort=10000;
for s=1:nsort
    vector(s,1)=sum(rand(n,1)>p);
end
edges=0:10;
=======
n=10;
p=0.5;
nsort=10000;
for s=1:nsort
    vector(s,1)=sum(rand(n,1)>p);
end
edges=0:10;
>>>>>>> f6aa4e5e0f096b95898ede76ce749f7364d07626
bar(edges,histc(vector,edges)/nsort)