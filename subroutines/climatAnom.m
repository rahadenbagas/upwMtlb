function [varClimat,varAnom] = climatAnom (var, time)

[~,month,~] = datevec(time); 

for k = 1:12
   ind = month==k; % Indices of month k
   varClimat(:,:,k) = nanmean(var(:,:,ind),3);
   varAnom(:,:,ind) = bsxfun(@minus,var(:,:,ind),varClimat(:,:,k));
end
