%%DistanceMatrix: r^2
function DM = DistanceMatrix(dsites, ctrs)
[M,s] = size(dsites); N = length(ctrs); 
DM = zeros(M,N);
for i=1:s
    [dr,cc]=ndgrid(dsites(:,i),ctrs(:,i));
    DM = DM +(dr-cc).^2;
end
DM=sqrt(DM);
return

