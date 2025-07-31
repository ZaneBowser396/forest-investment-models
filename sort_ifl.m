function sorted = sort_ifl(filtered)

IFL_IDs = {filtered.IFL_ID};
primaryNumbers = zeros(size(IFL_IDs));
secondaryNumbers = zeros(size(IFL_IDs));
for i = 1:length(IFL_IDs)
    parts = sscanf(IFL_IDs{i}, 'SEA_%d_%d');
    primaryNumbers(i) = parts(1);
    if length(parts) > 1
        secondaryNumbers(i) = parts(2);
    else
        secondaryNumbers(i) = 0;
    end
end

numericMatrix = [primaryNumbers', secondaryNumbers'];
[~, sortIdx] = sortrows(numericMatrix);
sorted = filtered(sortIdx);

end