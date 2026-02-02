function res = ballsum(a)
%this helps a bit with precision, without costing too much in computing
%res=sum(a);
res=sort(sum(a(a>0)))+sort(sum(a(a<0)),'descend');
end
