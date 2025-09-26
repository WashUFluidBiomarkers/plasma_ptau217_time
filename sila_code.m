function sila_ADNI
% Originally created by Kellen Petersen
% Updated by Kellen Petersen, September 12, 2025
% SILA analysis for ADNI dataset - estimates age at plasma p-tau217 positivity

% Add SILA toolbox to path
addpath(pwd + "/SILA-AD-Biomarker-main")
addpath(pwd + "/SILA-AD-Biomarker-main/demo")

% Set dataset name and load corresponding data
dataset_name = 'ADNI';  % Options: 'KADRC' or 'ADNI'

% Load dataset based on selection
if strcmp(dataset_name, 'KADRC')
    df = readtable("");
elseif strcmp(dataset_name, 'ADNI')
    df = readtable("");
else
    error('Invalid dataset name. Choose either ''KADRC'' or ''ADNI''');
end

% Select relevant columns and sort by participant and age
df = df(:, {'ID', 'EXAMDATE', 'C2N_plasma_ptau217_ratio', 'AGE'});
df = sortrows(df, {'ID','AGE'});

% Initialize analysis variables
df.ref_EAOA(:) = false;  % Reference point for onset age estimation
udf = unique(df.ID);
ptau_threshold = 4.06;   % Plasma p-tau217 positivity threshold
df.one_pos = zeros(height(df), 1);

% Identify participants who ever exceed the biomarker threshold
for i = 1:length(udf)
    current_id = udf(i);
    id_rows = df.ID == current_id;
    if any(df.C2N_plasma_ptau217_ratio(id_rows) >= ptau_threshold)
         df.one_pos(id_rows) = 1;
    end
end

% Calculate visit numbers and mark reference scans
for i = 1:numel(udf)
    ids = df.ID == udf(i);
    df.nvis(ids) = nnz(ids);
    df.visnum(ids) = 1:nnz(ids);
    % Use last scan as reference point for onset age estimation
    df.ref_EAOA(find(ids,1,'last')) = true;
end

% Train SILA model using optimized parameters
[tsila,tdrs] = SILA(df.AGE, ...
    df.C2N_plasma_ptau217_ratio, ...
    df.ID, ...
    0.1, ...           % Learning rate
    ptau_threshold, ...
    500);             % Number of iterations

% Estimate onset ages using trained SILA model
test_df = SILA_estimate(tsila, ...
    df.AGE(df.ref_EAOA), ...
    df.C2N_plasma_ptau217_ratio(df.ref_EAOA), ...
    df.ID(df.ref_EAOA), ...
    'align_event', 'last');

% Flag extrapolated estimates beyond training data range
test_df.extrapolated = test_df.val > max(tsila.val);

% Merge onset age estimates with original data
df_out = outerjoin(df, test_df, ...
    'Type','left', ...
    'leftKeys','ID', ...
    'RightKeys','subid', ...
    'RightVariables',{'estaget0','ageref','minage','maxage','truncated','extrapolated'});

% Create visualization plots
create_sila_plots(df, df_out, tsila, tdrs, ptau_threshold);

% Prepare output data structures
clock = tsila(:,{'val','time','adtime'});
clock = renamevars(clock, "val", "plasma");

df_output = df_out(:,{'ID','EXAMDATE','AGE','C2N_plasma_ptau217_ratio', 'estaget0', 'one_pos'});
df_output = renamevars(df_output, "C2N_plasma_ptau217_ratio", "plasma");

% Save results to files
save_sila_results(clock, df_output, tdrs, dataset_name);

end

function create_sila_plots(df, df_out, tsila, tdrs, ptau_threshold)
% Generate standard SILA analysis plots

figure(1000);
clf()
set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 24, 16]);

% Plot 1: Biomarker trajectories by chronological age
subplot(2,2,1)
spaghetti_plot2(df.AGE, df.C2N_plasma_ptau217_ratio, df.ID, df.one_pos)
hold on, plot(xlim, ptau_threshold*[1,1], '--k')
title('Plasma p-tau217 by Age')
xlabel('Age (years)')
ylabel('Plasma p-tau217 (pg/mL)')

% Plot 2: Biomarker trajectories by SILA time
subplot(2,2,2)
spaghetti_plot2(df_out.AGE-df_out.estaget0, df.C2N_plasma_ptau217_ratio, df.ID, df.one_pos)
hold on
plot(tsila.adtime, tsila.val, '-k', 'LineWidth', 4)
hold on, plot(xlim, ptau_threshold*[1,1], '--k')
title('Plasma p-tau217 by SILA Time')
xlabel('SILA time (years)')
ylabel('Plasma p-tau217 (pg/mL)')

% Plot 3: Discrete rate sampling curve
subplot(2,2,3)
plot(tdrs.val, tdrs.rate, '-'), hold on
plot(tdrs.val, tdrs.rate + tdrs.ci, '--r')
plot(tdrs.val, tdrs.rate - tdrs.ci, '--r')
title('Rate of Change vs Biomarker Level')
xlabel('Plasma p-tau217 (pg/mL)')
ylabel('Rate of Change per Year')

% Plot 4: SILA modeled trajectory
subplot(2,2,4)
plot(tsila.adtime, tsila.val, '-'), hold on
plot(xlim, ptau_threshold*[1,1], '--k')
title('SILA Model Trajectory')
xlabel('Time from Threshold (years)')
ylabel('Plasma p-tau217 (pg/mL)')
legend({'Modeled curve','Threshold'}, 'Location', 'northwest')

end

function save_sila_results(clock, df_output, tdrs, dataset_name)
% Save SILA analysis results to CSV files

sila_folder = '~/WASHU/WASHU_Projects/039_FNIH_FINAL/results_SILA';
timestamp = '2025-04-27';

% Save clock lookup table
writetable(clock, sila_folder + "/SILA_clock_" + dataset_name + "_" + timestamp + ".csv")

% Save individual participant data with onset age estimates
writetable(df_output, sila_folder + "/SILA_df_" + dataset_name + "_" + timestamp + ".csv")

% Save discrete rate sampling results
writetable(tdrs, sila_folder + "/SILA_tdrs_" + dataset_name + "_" + timestamp + ".csv")

end

function spaghetti_plot2(age, val, subid, ref_EAOA, varargin)
% Create spaghetti plot showing individual biomarker trajectories
% Colors trajectories based on whether participant ever exceeded threshold

tin = table(age, val, subid, ref_EAOA, 'VariableNames', {'age', 'val', 'subid', 'ref_EAOA'});
tin = sortrows(tin, {'subid', 'age'});

% Apply optional filtering
if nargin == 5
    tin = tin(varargin{1}, :);
end

subs = unique(tin.subid);

% Plot individual trajectories
for i = 1:numel(subs)
    ts = tin(tin.subid == subs(i), :);
    
    % Color coding: grey for never positive, red for ever positive
    if all(~ts.ref_EAOA)
        plot(ts.age, ts.val, '.-', 'Color', [0.5 0.5 0.5]), hold on
    else
        plot(ts.age, ts.val, '.-', 'Color', [0.85 0.1 0.1]), hold on
    end
end

end
