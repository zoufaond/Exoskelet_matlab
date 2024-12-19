clc
aaa = [1,2,3,4,5,6, NaN,8,9;
       1,2,3,4,5,6,7,8,9];

for k = 1:2
    res = funn(aaa(k,:));
    % aaa(k,:)
end

function res = funn(aaa)
for i = 1:length(aaa)
    aaa(i);
    if anynan(aaa(i))
        error('Nan')
    end
    res = aaa(i)
end
end