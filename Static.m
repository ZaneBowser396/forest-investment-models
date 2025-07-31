%%% Static Solution
%% Ranking
clear, clc

% Read data from Excel file using readtable
filename = 'tabulated_ifl_data.xlsx';   
dataTable = readtable(filename);

% Sort dataTable by 'OverallRankingValue' in descending order
sortedDataTable = sortrows(dataTable, 'OverallRankingValue', 'descend');

% Assign new ranking positions
sortedDataTable.RankingOrder = (1:height(sortedDataTable))';

% Create a table with IDs and their new ranking order
rankingTable = sortedDataTable(:, {'IFL_ID', 'RankingOrder'});

% Overwrite the original Excel file with the sorted data in Sheet 1
writetable(sortedDataTable, filename, 'Sheet', 1);

% Display the ranking table in MATLAB
disp('IDs and their Ranking Order:');
disp(rankingTable);

% Confirmation message
disp('Original Excel file has been overwritten with the sorted data and new ranking order.');

%% Partial Allocation

% Set number of top IFLs to be selected for further inspections
topN = 10;

% Read in newly sorted data from excel sheet
rankedData = readtable(filename, 'Sheet', 'Sheet1');

selectedIFLs = rankedData(1:topN, :);

% Extract variables from the table using column names
IDs = selectedIFLs.IFL_ID;      
B_i = selectedIFLs.Species_Richness;      
p_i = selectedIFLs.LocalSupport;
q_i = selectedIFLs.ChanceOfSuccess;
c_i = selectedIFLs.Travel_Time;
T_i = selectedIFLs.A_2020;
I = 1;
z = 0.3;

% Compute K_i
K_i = B_i .* ((p_i .* q_i * I) ./ (c_i .* T_i)) .^ z;

% Compute alpha
alpha = 1 / (1 - z);

% Compute numerator and denominator
numerator = (K_i * z) .^ alpha;
denominator = sum(numerator);

% Compute u_i
u_i = numerator / denominator;

% Verify the budget constraint
sum_u = sum(u_i);
fprintf('Sum of u_i: %.6f\n', sum_u);

% Compute total biodiversity benefit
beta_total = sum(K_i .* u_i .^ z);
fprintf('Total biodiversity benefit (beta): %.6f\n', beta_total);

%% Sorted Table to Track Changes from Ranking
% Assign original positions
original_positions = (1:length(IDs))';

% Create the resultsTable with original positions
resultsTable = table(IDs, u_i, original_positions, 'VariableNames', ...
    {'ID', 'OptimalAllocation', 'OriginalPosition'});

% Sort the resultsTable by 'OptimalAllocation' in descending order
sortedResults = sortrows(resultsTable, 'OptimalAllocation', 'descend');

% Assign new positions after sorting
sortedResults.NewPosition = (1:height(sortedResults))';

% Compute the change in position
sortedResults.ChangeInPosition = sortedResults.OriginalPosition - sortedResults.NewPosition;

% Rearrange columns for better readability
sortedResults = sortedResults(:, {'ID', 'OptimalAllocation', ...
    'OriginalPosition', 'NewPosition', 'ChangeInPosition'});

% Select only the top N IFLs for saving
topSortedResults = sortedResults(1:topN, :);

% Display the top sorted results in MATLAB
disp(['Top ', num2str(topN), ' Optimal Allocations with Change in Position:']);
disp(topSortedResults);

%% Compute Biodiversity Benefit β_i for Each IFL
% Compute β_i (contribution of each IFL to the total benefit)
beta_i = K_i .* u_i .^ z;

% Compute total biodiversity benefit (β_total)
beta_total = sum(beta_i);
fprintf('Total biodiversity benefit (β_total): %.6f\n', beta_total);

%% Create a Summary Table with β_i
partialAllocationTable = table(IDs, K_i, u_i, beta_i, ...
    'VariableNames', {'ID', 'K_i', 'OptimalAllocation', 'Beta_i'});

% Save the partial allocation results to a new sheet in the Excel file
writetable(partialAllocationTable, filename, 'Sheet', 'Top_IFL_Allocation');

% Display the summary table
disp('Partial Allocation Results with Biodiversity Benefit (β_i):');
disp(partialAllocationTable);

%% Clear Previous Content from 'Optimal_Allocation' Sheet

% Create a blank table to clear old data from the sheet
writetable(table(), filename, 'Sheet', 'Optimal_Allocation', 'WriteMode', 'overwrite');

%% Save the Top N Results to the 'Optimal_Allocation' Sheet

% Save only the top N sorted results to the 'Optimal_Allocation' sheet
writetable(topSortedResults, filename, 'Sheet', 'Optimal_Allocation');

% Confirmation message
disp(['Top ', num2str(topN), ' optimal allocation results saved to the "Optimal_Allocation" sheet.']);

%% Scatter Plots to Determine Effects on Optimal Allocation

% B_i vs u_i
figure;
scatter(B_i, u_i, 'filled');
xlabel('Biodiversity Benefit (B_i)');
ylabel('Optimal Allocation (u_i)');
title('Relationship Between Biodiversity and Optimal Allocation');
hold on;
p = polyfit(B_i, u_i, 1); 
y_fit = polyval(p, B_i);
plot(B_i, y_fit, 'r--', 'LineWidth', 1.5);
legend('Data Points', 'Trend Line');
hold off;

% c_i vs u_i
invc_i = 1 ./ c_i;
figure;
scatter(invc_i, u_i, 'filled');
xlabel('Cost (1 / c_i)');
ylabel('Optimal Allocation (u_i)');
title('Relationship Between Cost and Optimal Allocation');
hold on;
p = polyfit(invc_i, u_i, 1);  
y_fit = polyval(p, invc_i);
plot(invc_i, y_fit, 'r--', 'LineWidth', 1.5);
legend('Data Points', 'Trend Line');
hold off;

% Plot Optimal Allocation (u_i) vs. Biodiversity Benefit (β_i)
figure;
scatter(u_i, beta_i, 'filled');

% Add labels and title to the plot
xlabel('Optimal Allocation (u_i)');
ylabel('Biodiversity Benefit (β_i)');
title('Optimal Allocation vs Biodiversity Benefit');

% Add a best-fit line
hold on;
p = polyfit(u_i, beta_i, 1);  % Linear fit (degree 1 polynomial)
y_fit = polyval(p, u_i);  % Evaluate the polynomial fit
plot(u_i, y_fit, 'r--', 'LineWidth', 1.5);  % Plot the best-fit line

% Add legend and grid
legend('Data Points', 'Best-Fit Line', 'Location', 'best');
grid on;
hold off;

%% Pie Chart of Optimal Allocations (u_i)
figure;

% Generate the pie chart using optimal allocation (u_i)
h = pie(u_i);

% Convert IDs to strings and escape underscores for proper display
labels = replace(string(IDs), "_", "\_");  % Replace underscores with escaped versions

% Use a preinstalled colormap (e.g., 'parula')
colormapName = 'parula';  % Choose any MATLAB colormap: 'parula', 'jet', 'turbo', etc.
numColors = length(u_i);  % Number of slices

% Extract colors from the chosen colormap
colors = colormap(colormapName);  % Get full colormap
customColors = colors(round(linspace(1, size(colors, 1), numColors)), :);  % Sample required colors

% Apply the colormap to the pie slices
for k = 1:2:length(h)  % Pie slices are every other element in 'h'
    idx = (k + 1) / 2;  % Calculate the index for the color array
    h(k).FaceColor = customColors(idx, :);
end

% Adjust the text objects in the pie chart with the corresponding labels
for k = 1:length(u_i)
    hText = h(k * 2);  % Every second element is a text object
    hText.String = labels{k};  % Assign the escaped label
    hText.FontSize = 12;  % Optional: Increase font size for better readability
end

% Set the title for the pie chart
title('Pie Chart of Optimal Allocations', 'FontSize', 14);

%% Uncertainty Analysis

% Read the optimal allocation data from the 'Optimal_Allocation' sheet
filename = 'tabulated_ifl_data.xlsx';  
optimalData = readtable(filename, 'Sheet', 'Optimal_Allocation');

% Extract IDs, original optimal allocations, and their original positions
IDs = optimalData.ID;
originalAllocations = optimalData.OptimalAllocation;

% Use the row order in the original data as the "original positions"
originalPositions = (1:height(optimalData))';

% Number of trials for uncertainty analysis
numTrials = 100;

% Initialize storage for rank changes and shifted allocations across trials
rankingChanges = zeros(height(optimalData), numTrials);
allShiftedAllocations = zeros(height(optimalData), numTrials);

for trial = 1:numTrials
    % Randomly shift the optimal allocations by 1-5%
    randomShift = 1 + (0.01 + (0.05 - 0.01) * rand(height(optimalData), 1));
    shiftedAllocations = originalAllocations .* randomShift;

    % Store the shifted allocations for this trial
    allShiftedAllocations(:, trial) = shiftedAllocations;

    % Create a table with the shifted allocations (same order as original)
    shiftedTable = table(IDs, shiftedAllocations, ...
        'VariableNames', {'ID', 'ShiftedAllocation'});

    % Find the new rank positions based on shifted allocations
    [~, newOrder] = sort(shiftedAllocations, 'descend');
    newPositions = zeros(height(shiftedTable), 1);
    newPositions(newOrder) = 1:height(shiftedTable);  % Assign new ranks

    % Store the rank changes relative to the original positions
    rankChange = originalPositions - newPositions;
    rankingChanges(:, trial) = rankChange;
end

%% Analyze the Results

% Compute the mean and standard deviation of rank changes across trials
meanRankChange = mean(rankingChanges, 2);
stdRankChange = std(rankingChanges, 0, 2);

% Compute the mean and standard deviation of shifted allocations
meanShiftedAllocations = mean(allShiftedAllocations, 2);

% Create a summary table with results
summaryTable = table(IDs, originalAllocations, meanShiftedAllocations, ...
    meanRankChange, stdRankChange, ...
    'VariableNames', {'ID', 'OriginalAllocation', 'MeanShiftedAllocation', ...
                       'MeanRankChange', 'StdRankChange'});

%% Clear Previous Content from 'Uncertainty_Analysis' Sheet

% Create a blank table to clear old data from the sheet
writetable(table(), filename, 'Sheet', 'Uncertainty_Analysis', 'WriteMode', 'overwrite');

%% Save the Results to the 'Uncertainty_Analysis' Sheet

% Save the summary table to the 'Uncertainty_Analysis' sheet
writetable(summaryTable, filename, 'Sheet', 'Uncertainty_Analysis');

% Display the summary table
disp('Uncertainty Analysis Results:');
disp(summaryTable);
