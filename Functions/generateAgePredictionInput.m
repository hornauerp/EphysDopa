function [mat,age,cl] = generateAgePredictionInput(input,tp,age,cl)
sample_idx = repmat(1:tp,1,ceil(size(input,1)/tp^2));
sample_idx = sample_idx(randperm(size(input,1)/tp))+(0:5:height(input)-1);
mat = input(sample_idx,:);
age = age(sample_idx);
cl = cl(sample_idx);
end