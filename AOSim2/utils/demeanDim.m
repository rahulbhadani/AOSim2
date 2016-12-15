function V = demeanDim(V,dim)

% V = demeanDim(V,dim)
%  remove the mean in the specified direction. 
%  JLCodona: 20160909

if(nargin<2)
    dim = 1;
end

if(dim==1)
    V = V - repmat(mean(V,1),[size(V,1) 1]);
else
    V = V - repmat(mean(V,2),[1 size(V,2)]);
end

return;
