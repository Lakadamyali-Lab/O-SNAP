% -------------------------------------------------------------------------
% preprocess_OSNAP_feature_data.m
% -------------------------------------------------------------------------
% Performs different preprocessing on O-SNAP feature data based on flags
% passed to the function, including normalization and filtering.
%
% Example on how to use it:
%   [feature_data_norm,group_values] = preprocess_OSNAP_feature_data(feature_data)
% -------------------------------------------------------------------------
% Input:
%   feature_data: The feature data table, where each row represents 
%                 a sample(nucleus). The first three columns represent (1) 
%                 Group/Phenotype, (2) replicate, and (3) Sample 
%                 Identifier. Each subsequent column is  an O-SNAP feature.
% Output:
%   feature_data_norm: Output table with corresponding modifications based
%                      on flags
%   group_values: Returns the group/phenotype labels for every sample. Is
%                 helpful when numeric_only flag is set to true
% Options:
%   groups: Cell array containing char array of the identifiers of the 
%           phenotypes/cell states
%   replicates: Cell array containing char array of the identifiers of the 
%               replicates.
%   normalize: Flag to indicate whether data should be z-score normalized
%              across features
%   remove_NaN: Flag to indicate whether NaN values should be removed
%               (data is removed row-wise if a value is missing, 
%               i.e. the entire sample is omitted)
%   keep_rep_sample_info: Flag to indicate whether replicate information
%                         should be preserved in output
%   numeric_only: Flag to indicate whether only numeric values should  be
%                 preserved in output
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   hannah.kim3@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   H. H. Kim, J. A. Martinez-Sarmiento, F. R. Palma, A. Kant, E. Y. Zhang,
%   Z. Guo, R. L. Mauck, S. C. Heo, V. Shenoy, M. G. Bonini, M. Lakadamyali,
%   O-SNAP: A comprehensive pipeline for spatial profiling of chromatin
%   architecture. bioRxiv, doi: 10.1101/2025.07.18.665612 (2025).
% -------------------------------------------------------------------------
function [feature_data_norm,group_values] = preprocess_OSNAP_feature_data(feature_data,options)
arguments
    feature_data table
    options.groups cell = {};
    options.replicates cell = {};
    options.normalize logical = true;
    options.remove_NaN logical = true;
    options.keep_rep_sample_info logical = false;
    options.numeric_only logical = false; 
end
feature_data_norm = feature_data;
%% Select only desired groups and replicates for analysis
if ~isempty(options.replicates)
    feature_data_norm = feature_data_norm(ismember(feature_data_norm.biological_replicate,options.replicates),:);
end
if ~isempty(options.groups)
    feature_data_norm = feature_data_norm(ismember(feature_data_norm.group,options.groups),:);
end
if ismember('group', feature_data_norm.Properties.VariableNames)
    group_values = feature_data_norm.group;
else
    group_values = [];
end
%% Remove unnecessary columns (group, bio replicate)
if options.numeric_only
    feature_data_norm = feature_data_norm(:,vartype('numeric'));
elseif ~options.keep_rep_sample_info
    if ismember('biological_replicate',feature_data_norm.Properties.VariableNames)
        feature_data_norm(:,'biological_replicate') = [];
    end
    if ismember('name',feature_data_norm.Properties.VariableNames)
        feature_data_norm(:,'name') = [];
    end
end

%% Normalize T
if options.normalize
    % T_norm{:,vartype('numeric')} = normalize(abs(T_norm{:,vartype('numeric')}),1);
    feature_data_norm{:,vartype('numeric')} = normalize(feature_data_norm{:,vartype('numeric')},1);
    % remove NaNs
    if options.remove_NaN
        is_nan_all_row = all(ismissing(feature_data_norm{:,vartype('numeric')}),2);
        feature_data_norm = feature_data_norm(~is_nan_all_row,:);
        if ~isempty(group_values)
            group_values = group_values(~is_nan_all_row);
        end
        is_nan_any_col = any(ismissing(feature_data_norm),1);
        feature_data_norm = feature_data_norm(:,~is_nan_any_col);
        is_nan_any_row = any(ismissing(feature_data_norm{:,vartype('numeric')}),2);
        feature_data_norm = feature_data_norm(~is_nan_any_row,:);
        if ~isempty(group_values)
            group_values = group_values(~is_nan_any_row);
        end
    end
end
%% Remove NaNs
if options.remove_NaN
    is_nan_all_row = all(ismissing(feature_data_norm{:,vartype('numeric')}),2);
    feature_data_norm = feature_data_norm(~is_nan_all_row,:);
    if ~isempty(group_values)
        group_values = group_values(~is_nan_all_row);
    end
    is_nan_any_col = any(ismissing(feature_data_norm),1);
    feature_data_norm = feature_data_norm(:,~is_nan_any_col);
    is_nan_any_row = any(ismissing(feature_data_norm{:,vartype('numeric')}),2);
    feature_data_norm = feature_data_norm(~is_nan_any_row,:);
    if ~isempty(group_values)
        group_values = group_values(~is_nan_any_row);
    end
end
end

