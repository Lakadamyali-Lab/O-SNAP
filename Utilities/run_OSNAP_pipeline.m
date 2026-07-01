% -------------------------------------------------------------------------
% run_OSNAP_pipeline.m
% -------------------------------------------------------------------------
% Function that runs the suite of analysis available in O-SNAP, excluding
% the pseudotime analysis that implemented in R.
% 
% Most descriptors/features that are calculated are calculated with respect to
% the point cloud or an accurate representation of a pointcloud (based on
% the polyshape functions of MATLAB).
%
% Example on how to use it:
%   results = run_OSNAP_pipeline("C:\Users\JohnSmith\Documents\O-SNAP_Analysis\,...
%                      "SmithJohn_HeLa_ProteinKO_H2BHistone",...
%                      {'Control','KO'},...
%                      {'20250101','20250108','20250201'},...
%                      "run_extract_features",1,
%                      "run_generate_table",1,...
%                      "run_comparison",1,...
%                      "run_generate_batches",1,
%                      "run_feature_selection",1,
%                      "run_PCA",1, ...
%                      "run_classification_batch",1, ...
%                      "run_plot_violin",1,...
%                      "run_plot_radial",1,...
%                      "run_FSEA",1,...
%                      "save",1);
% -------------------------------------------------------------------------
% Input:
%   root_dir: The parent directory that contains the directory where the
%             analysis results are saved
%   analysis_name: Analysis identifier string
%   groups: Cell array containing char array of the identifiers of the
%           phenotypes/cell states
%           *ENSURE THAT EACH SAMPLE FILENAME CONTAINS EXACTLY ONE
%            IDENTIFIER FROM groups*
%   replicates: Cell array containing char array of the identifiers of the
%               replicates.
%               *ENSURE THAT EACH REPLCIATE IS STORED IN A SEPARATE FOLDER
%                PREFIXED BY "Rep"*
% Output:
%   OSNAP_data: A struct containing the following fields related to the
%              O-SNAP analysis:
%                 - analysis_name: Analysis identifier string
%                 - classification_summary_all: Table of classification
%                                               accuracies using "all"
%                                               approach
%                 - classification_summary_batch_all: Table of
%                                                     classification
%                                                     accuracies using
%                                                     "batch_all" approach
%                 - classification_summary_batch_each: Table of
%                                                      classification
%                                                      accuracies using
%                                                      "batch_each"
%                                                      approach
%                 - classifiers_all: Cell array with all instances of
%                                    trained models using "all" approach
%                 - classifiers_batch_all: Cell array with all instances
%                                          of trained models using
%                                          "batch_all" approach
%                 - classifiers_batch_each: Cell array with all instances
%                                           of trained models using
%                                           "batch_each" approach
%                 - date: The most recent date the analysis was run
%                 - feature_comparisons: A cell array where every cell is
%                                        a pair-wise combination of the
%                                        fold-change analysis used in the
%                                        volcano plots
%                 - feature_data: The table containing the O-SNAP feature
%                                 values where each row represents a
%                                 sample (nucleus) and each column is an
%                                 O-SNAP feature
%                 - feature_set_coverage: A cell array with a breakdown of
%                                         how many features in a set are
%                                         deemed significant by pairwise
%                                         comparison
%                 - groups: Cell array containing char array of the
%                           identifiers of the phenotypes/cell states
%                           *ENSURE THAT EACH SAMPLE FILENAME CONTAINS
%                           EXACTLY ONE IDENTIFIER FROM groups*
%                 - options: Struct array of parameter values for
%                            O-SNAP analysis
%                 - pca_result_all: Info on PCA transformation in the
%                                   classification pipeline using the
%                                   "all" approach
%                 - pca_result_batch_all: Info on the PCA transformation
%                                         in the classification pipeline
%                                         using the "batch_all" approach
%                 - pca_result_batch_each: Info on the PCA transformation
%                                         in the classification pipeline
%                                         using the "batch_each" approach"
%                 - replicates: Cell array containing char array of the
%                               identifiers of the replicates.
%                               *ENSURE THAT EACH REPLCIATE IS STORED IN
%                                A SEPARATE FOLDER PREFIXED BY "Rep"*
%                 - starttime: The starttime identifier for the analysis
%                              run
%                 - test_idxs: A cell array where each cell contains a
%                              logical array indicating samples for the
%                              test data of each fold
%                 - train_idxs: A cell array where each cell contains a
%                               logical array indicating samples for the
%                               training data of each fold
%                 - vars_select_result_all: MRMR scores of features for
%                                           ranking performed on entire
%                                           feature data
%                 - vars_select_result_batch: Table with MRMR scores of
%                                             each fold and the sum of the
%                                             scores across all folds
%                 - vars_selected_all: Features selected from MRMR
%                                      performed on entire dataset
%                 - vars_selected_batch: Table of boolean values
%                                        indicating whether a feature is
%                                        selected within a batch and from
%                                        the aggregated selection
%                 - venn_data: Data for Venn diagram to compare changes
%                              between 3-4 comparisons of phenotype pairs
% Options
%   run_extract_features: Calculate relevant features from each OSNAP
%                         sample MAT file
%   run_generate_table: Create the feature table where each row is an OSNAP
%                       sample and each column is a feature
%   run_comparison: Perform pair-wise comparison of OSNAP features between
%                   different phenotypes
%   run_generate_batches: When performing classification pipeline, create
%                         multiple folds to validate model performance
%   run_feature_selection: Perform feature selection for classification
%   run_PCA: PCA-transform data
%   run_classification_all: Create a classification model using all samples
%                           for training (including feature selection and
%                           PCA transformation), using 5-fold cross
%                           validation to evaluate model. Vulnerable to
%                           data leakage
%   run_classification_batch: Perform classification using multiple folds,
%                             both totally independently and using an 
%                             aggregate of feature selection results
%   run_venn: If there are 3-4 phenotypes being compared, create a venn
%             diagram to see what features change in common to each other
%   run_plot_violin: Create violin plots comparing feature distributions by
%                    phenotype
%   run_plot_radial: Create images for every sample depicting the radial
%                    density analysis
%   run_FSEA: Perform Feature Set Enrichment Analysis
%   split_method: The method to generate the folds in classification
%   proportion: If using the "bootstrap" method to create batches, defines
%               the ratio of test to total samples
%   test_train_k: If using the "k-folds" or "bootstrap" method to create
%                 batches, defines the number of folds generated
%   feature_select_max_idx: The maximum number of features that may be
%                           included for feature selection
%   num_components_explained: Defines the cutoff for the number of PCA
%                             components to keep. If < 1, the parameter is 
%                             the minimum amount of variance that the
%                             number of components must explain. If > 1,
%                             then it is the number of components that are
%                             kept.
%   n_models_per_type: The number of instances of a model of a given
%                      architecture to train.
%   alpha: The p-value cutoff for significance testing
%   fold_change_threshold: The fold change cut off when performing
%                          pair-wise analysis of groups
%   feature_universe_names: If multiple feature universes (ie. 
%                           categorization of features into sets) are 
%                           defined for FSEA, this parameter is the list of
%                           identifiers for those to include
%   FSEA_rank_type: The ranking metric approach
%   n_processes: The number of workers to use for parallel processing
%   filter: Flag as to whether to filter out samples that do not contain
%           either the specified phenotype or replicate labels
%   suffix: A suffix to the identifier for the analysis
%   save: Logical to save the analysis output to a file
%   save_if_error: Logical on whether to save the analysis if it encounters
%                  an error
%   check_overwrite: Check if a file is going to be overwritten and wait
%                    for user input
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
function OSNAP_data = run_OSNAP_pipeline(root_dir,analysis_name,groups,replicates,options)
arguments
    root_dir char
    analysis_name char
    groups cell
    replicates cell
    % run options (will overwrite)
    options.run_extract_features logical = 1
    options.run_generate_table logical = 1
    options.run_comparison logical = 1
    options.run_generate_batches logical = 1
    options.run_feature_selection logical = 1
    options.run_PCA logical = 1
    options.run_classification_all logical = 1
    options.run_classification_batch logical = 1
    options.run_venn logical = 0
    options.run_plot_violin logical = 1
    options.run_plot_radial logical = 1
    options.run_FSEA logical = 0
    % generate batch options
    options.split_method string = "k-fold"
    options.proportion double = 0.2;
    options.test_train_k double = 5;
    % feature select options
    options.feature_select_max_idx = 12;
    % PCA options
    options.num_components_explained double = 0.75
    % classification options
    options.n_models_per_type = 5;
    % comparison options
    options.alpha double = 0.05
    options.fold_change_threshold double = 2;
    % FSEA options
    options.feature_universe_names = ["universe_1"];
    options.FSEA_rank_type string = "S2N";
    % other options
    options.n_processes double = maxNumCompThreads;
    options.filter logical = true;
    options.suffix = ""
    options.save = 1;
    options.save_if_error = 1;
    options.check_overwrite = 0
end
%% Setup O-SNAP run
work_dir = fullfile(root_dir,analysis_name);
if options.suffix == ""
    save_analysis_path = fullfile(work_dir,analysis_name+".mat");
else
    save_analysis_path = fullfile(work_dir,analysis_name+"_"+options.suffix+".mat");
end
log_file = fullfile(work_dir,sprintf('%s_%s.txt',analysis_name,datetime('now','Format','y-MM-dd_HH-mm-ss')));
diary off
diary(log_file)
fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
fprintf("Running analysis for: " + analysis_name + "\n")
fprintf("RUN START: %s\n", string(datetime))
fprintf("   GROUPS: %s\n", string(join(groups(:),', ')))
fprintf("   REPS: %s\n", string(join(replicates(:),', ')))
warning('off','all');
close('all');
rng('default');
%% Check which analyses to run depending on the run options
try
    % Load existing data
    if exist(save_analysis_path,"file")
        starttime_step = tic;
        fprintf("  Loading data from %s...\n",save_analysis_path)
        try
            OSNAP_data = load(save_analysis_path);
        catch ME
            OSNAP_data = struct;
            handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"stop_run",0,"save",0);
            options.run_extract_features = 1;
        end
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        OSNAP_data.analysis_name = analysis_name;
        OSNAP_data.groups = unique(groups);
        OSNAP_data.replicates = unique(replicates);
        % Generate O-SNAP feature table data if not loaded
        if isfield(OSNAP_data,'feature_data')
            feature_data = OSNAP_data.feature_data;
        else
            options.run_generate_table = 1;
        end
        % Run feature comparisons if not loaded
        if ~isfield(OSNAP_data,'feature_comparisons')
            options.run_comparison = 1;
        end
        % Run steps prior to classification
        if options.run_classification_batch
            % Batch generation
            if ~all(isfield(OSNAP_data,{'train_idxs','test_idxs'}))
                options.run_generate_batches = 1;
            elseif isempty(OSNAP_data.train_idxs) || isempty(OSNAP_data.test_idxs)
                options.run_generate_batches = 1;
            end
            % Feature selection
            if ~all(isfield(OSNAP_data,{'vars_select_result_batch','vars_selected_batch'}))
                options.run_feature_selection = 1;
            elseif isempty(OSNAP_data.vars_select_result_batch) || isempty(OSNAP_data.vars_selected_batch)
                options.run_feature_selection = 1;
            end
            % PCA
            if ~isfield(OSNAP_data,'pca_result_batch_each')
                options.run_PCA = 1;
            elseif isempty(OSNAP_data.pca_result_batch_each)
                options.run_PCA = 1;
            end
            if ~isfield(OSNAP_data,'pca_result_batch_all')
                options.run_PCA = 1;
            elseif isempty(OSNAP_data.pca_result_batch_all)
                options.run_PCA = 1;
            end
        end
        OSNAP_data.options = options;
    end
catch ME
    options.run_generate_table = 1;
    options.run_comparison = 1;
    options.run_generate_batches = 1;
    options.run_feature_selection = 1;
    options.run_PCA = 1;
    options.run_classification_batch = 1;
    OSNAP_data.options = options;
    if exist('OSNAP_data','var')
        handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    end
end
%% Get metadata for the run
OSNAP_data.date = datetime;
OSNAP_data.log_file = log_file;
OSNAP_data.starttime = tic;
OSNAP_data.groups = unique(groups);
OSNAP_data.replicates = unique(replicates);
%% Extract O-SNAP features per sample (nucleus)
try
    if options.run_extract_features
        fprintf("  Generating features from scratch...\n")
        % extract features
        starttime_step = tic;
        OSNAP_sample_file_list = get_valid_OSNAP_samples(work_dir, groups, replicates,...
            {'x','y'},'filter',options.filter);
        generate_OSNAP_features_samplelist(OSNAP_sample_file_list,'overwrite',false,'n_processes',options.n_processes,'filter',true);
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    fprintf("%s\n",getReport(ME));
    return
end
%% Generate O-SNAP feature table data
try
    if options.run_generate_table
        fprintf("  Creating table from generated features...\n");
        % coallate features
        starttime_step = tic;
        feature_data = extract_OSNAP_features_batch(work_dir,groups,replicates);
        OSNAP_data.feature_data = feature_data;
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
%% Filter data
if exist("feature_data","var")
    try
        feature_data_filtered = filter_OSNAP_feature_data(feature_data,groups,replicates);
        writetable(feature_data_filtered,replace(save_analysis_path,".mat",".csv"));
    catch ME
        handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
        return
    end

else
    fprintf("No feature data, ending run...\n");
    return
end
%% Violin plots
try
    if options.run_plot_violin
        fprintf("  Plotting violin plots...\n")
        starttime_step = tic;
        plot_OSNAP_violin(feature_data_filtered,work_dir)
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
%% Compare all groups to each other
try
    if options.run_comparison || (options.run_venn && ~isfield(OSNAP_data,"feature_comparisons"))
        fprintf("  Comparing features...\n")
        starttime_step = tic;
        save_path = fullfile(work_dir, "features_volcano");
        if options.suffix ~= ""
            save_path = join(save_path,options.suffix,"_");
        end
        OSNAP_data.feature_comparisons = compare_OSNAP_group_pair_batch(feature_data_filtered,...
            "alpha",options.alpha,...
            "fold_change_threshold",options.fold_change_threshold,...
            "save_path",save_path);
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
if isfield(OSNAP_data,"feature_comparisons")
    for i=1:numel(OSNAP_data.feature_comparisons)
        disp(join(OSNAP_data.feature_comparisons{i}.groups," vs "))
        fprintf("UP: %3.0f\n",sum(all([OSNAP_data.feature_comparisons{i}.feature_table{:,5}>1, OSNAP_data.feature_comparisons{i}.feature_table{:,7} < 0.05],2)))
        fprintf("DOWN: %3.0f\n",sum(all([OSNAP_data.feature_comparisons{i}.feature_table{:,5}<-1, OSNAP_data.feature_comparisons{i}.feature_table{:,7} < 0.05],2)))
    end
end
%% Venn diagram
try
    if options.run_venn
        fprintf("  Creating venn diagrams...\n")
        starttime_step = tic;
        OSNAP_data.venn_data = plot_OSNAP_venn(OSNAP_data.feature_comparisons,"save_path",fullfile(work_dir, "features_venn"));
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
%% Run classification steps if needed
feature_data_filtered = filter_OSNAP_feature_data(feature_data_filtered,groups,replicates,...
    'remove_NaN',true);
%% Generate train and test batches
try
    if options.run_generate_batches
        fprintf("  Creating train/test batches...\n")
        starttime_step = tic;
        [OSNAP_data.train_idxs, OSNAP_data.test_idxs,by] = generate_OSNAP_train_test_sets(...
            feature_data_filtered,...
            "split_method",options.split_method,...
            "k",options.test_train_k,...
            "proportion",options.proportion);
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
        options.split_method = by;
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
%% Feature selection
try
    if options.run_feature_selection
        % Clear variable for new run
        if isfield(OSNAP_data,'vars_select_result_batch')
            OSNAP_data = rmfield(OSNAP_data,"vars_select_result_batch");
        end
        fprintf("  Selecting features (batch-wise)...\n")
        starttime_step = tic;
        save_path = fullfile(work_dir, join(['feature_selection_batch',groups],'_'));
        % Split batches for parallel
        train_idx_p = split_data_to_n_OSNAP(OSNAP_data.train_idxs,options.n_processes,"shuffle",false);
        n_processes = numel(train_idx_p);
        feature_data_p = cell(1,n_processes);
        var_scores_p = cell(1,n_processes);
        for p=1:n_processes
            feature_data_p{p} = feature_data_filtered;
        end
        parfor p=1:n_processes
            train_idx_b = train_idx_p{p};
            for b=1:numel(train_idx_b)
                train_data = feature_data_p{p}(train_idx_b{b},:);
                [~, var_scores_p{p}{b}] = fscmrmr(train_data(:,vartype('numeric')),train_data.group);
            end
        end
        % Bring together batches from each parallel process
        scores = [var_scores_p{:}];
        scores = vertcat(scores{:})';
        % Store values
        [OSNAP_data.vars_select_result_batch,OSNAP_data.vars_selected_batch] = select_OSNAP_features( ...
            scores, ...
            feature_data_filtered(:,vartype('numeric')).Properties.VariableNames, ...
            "max_idx",options.feature_select_max_idx, ...
            "save_path",save_path);
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
if isfield(OSNAP_data,'vars_selected_batch')
    vars_display = OSNAP_data.vars_selected_batch.Properties.RowNames(logical(OSNAP_data.vars_selected_batch.batch_sum));
    fprintf("      Selected features (n = %.0f):\n",numel(vars_display))
    for i=1:numel(vars_display)
        fprintf("        - %s\n", vars_display{i})
    end
end
%% PCA
try
    if options.run_PCA
        save_path = fullfile(work_dir, string(join(['PCA',groups],'_')));
        if options.suffix ~= ""
            save_path = string(join([save_path options.suffix],'_'));
        end
        fprintf("  Performing PCA (batch-wise)...\n")
        starttime_step = tic;
        % Split batches for parallel
        num_components_explained = options.num_components_explained;
        train_idx_p = split_data_to_n_OSNAP(OSNAP_data.train_idxs,options.n_processes,"shuffle",false);
        n_processes = numel(train_idx_p);
        n_batch = numel(OSNAP_data.train_idxs);
        if ~exist('feature_data_p','var')
            feature_data_p = cell(1,n_processes);
            for p=1:n_processes
                feature_data_p{p} = feature_data_filtered;
            end
        end
        % PCA based on feature selection by EACH batch
        if isfield(OSNAP_data,'vars_select_result_batch')
            vars_selected_b = split_data_to_n_OSNAP(arrayfun(@(x) OSNAP_data.vars_selected_batch.Properties.RowNames(logical(OSNAP_data.vars_selected_batch{:,x})), 1:n_batch,'uni',0),n_processes,"shuffle",false);
        else
            vars_selected_b = split_data_to_n_OSNAP(cell(1,n_batch),n_processes,"shuffle",false);
        end
        pca_result_b = cell(1,n_processes);
        parfor p=1:n_processes
            train_idx_b = train_idx_p{p};
            for b=1:numel(train_idx_b)
                train_data = feature_data_p{p}(train_idx_b{b},:);
                pca_result_b{p}{b} = run_OSNAP_PCA(train_data,...
                    "vars_sel",vars_selected_b{p}{b},...
                    "save_path",sprintf("%s_%02.0f_%02.0f",save_path,p,b),...
                    "num_components_explained",num_components_explained);
            end
        end
        OSNAP_data.pca_result_batch_each = [pca_result_b{:}];
        % PCA based on feature selection on ALL batches
        if isfield(OSNAP_data,'vars_select_result_batch')
            vars_selected_b = split_data_to_n_OSNAP(repmat({OSNAP_data.vars_selected_batch.Properties.RowNames(logical(OSNAP_data.vars_selected_batch{:,end}))},1,n_batch),n_processes,"shuffle",false);
        else
            vars_selected_b = split_data_to_n_OSNAP(cell(1,n_batch),n_processes,"shuffle",false);
        end
        pca_result_b = cell(1,n_processes);
        parfor p=1:n_processes
            train_idx_b = train_idx_p{p};
            for b=1:numel(train_idx_b)
                train_data = feature_data_p{p}(train_idx_b{b},:);
                pca_result_b{p}{b} = run_OSNAP_PCA(train_data,...
                    "vars_sel",vars_selected_b{p}{b},...
                    "save_path",save_path+"_all",...
                    "num_components_explained",num_components_explained);
            end
        end
        OSNAP_data.pca_result_batch_all = [pca_result_b{:}];
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
end
%% Run classification batch-wise
try
    if options.run_classification_batch
        n_models_per_type = options.n_models_per_type;
        fprintf("  Running classification (batch-wise)...\n")
        starttime_step = tic;
        % cell_dims = [numel(OSNAP_data.train_idxs),numel(OSNAP_data.vars_selected_batch),numel(OSNAP_data.pca_result)];
        % assert(all(cell_dims == cell_dims(1)),'OSNAP:classification_length_mismatch','number of batches not consistent with classification structures (feature selection, PCA)')
        % Split batches for parallel
        n_batch = numel(OSNAP_data.train_idxs);
        % Classification based on feature selection by EACH batch
        if isfield(OSNAP_data,'vars_select_result_batch')
            vars_selected_b = arrayfun(@(x) OSNAP_data.vars_selected_batch.Properties.RowNames(logical(OSNAP_data.vars_selected_batch{:,x}))',1:n_batch,'uni',0);
        else
            vars_selected_b = cell(1,n_batch);
        end
        if isfield(OSNAP_data,'pca_result_batch_each')
            if ~isempty(OSNAP_data.pca_result_batch_each)
                pca_result_b = OSNAP_data.pca_result_batch_each;
            else
                pca_result_b = cell(1,n_batch);
            end
        else
            pca_result_b = cell(1,n_batch);
        end
        classifiers = cell(1,numel(OSNAP_data.train_idxs));
        for b=1:n_batch
            train_data = feature_data_filtered(OSNAP_data.train_idxs{b},:);
            test_data = feature_data_filtered(OSNAP_data.test_idxs{b},:);
            classifiers{b} = run_OSNAP_classification_batch(...
                train_data,...
                'n_models_per_type',n_models_per_type,...
                'vars_selected',vars_selected_b{b},...
                'pca_result',pca_result_b{b},....
                'test_data',test_data,...
                'verbose',1);
        end
        OSNAP_data.classifiers_batch_each = classifiers;
        % Plot figures
        if options.save
            OSNAP_data.classification_summary_batch_each = evaluate_OSNAP_classifiers(OSNAP_data.classifiers_batch_each,"save_path",fullfile(work_dir,"batch_each"));
        else
            OSNAP_data.classification_summary_batch_each = evaluate_OSNAP_classifiers(OSNAP_data.classifiers_batch_each);
        end
        % Classification based on feature selection on ALL batches
        if isfield(OSNAP_data,'vars_select_result_batch')
            vars_selected_b = OSNAP_data.vars_selected_batch.Properties.RowNames(logical(OSNAP_data.vars_selected_batch{:,end}))';
        else
            vars_selected_b = [];
        end
        if isfield(OSNAP_data,'pca_result_batch_all')
            if ~isempty(OSNAP_data.pca_result_batch_all)
                pca_result_b = OSNAP_data.pca_result_batch_all;
            else
                pca_result_b = cell(1,n_batch);
            end
        else
            pca_result_b = cell(1,n_batch);
        end
        classifiers = cell(1,n_batch);
        for b=1:n_batch
            train_data = feature_data_filtered(OSNAP_data.train_idxs{b},:);
            test_data = feature_data_filtered(OSNAP_data.test_idxs{b},:);
            classifiers{b} =...
                run_OSNAP_classification_batch(...
                train_data,...
                'n_models_per_type',n_models_per_type,...
                'vars_selected',vars_selected_b,...
                'pca_result',pca_result_b{b},....
                'test_data',test_data,...
                'verbose',0);
        end
        OSNAP_data.classifiers_batch_all = classifiers;
        if options.save
            OSNAP_data.classification_summary_batch_all = evaluate_OSNAP_classifiers(OSNAP_data.classifiers_batch_all,"save_path",fullfile(work_dir,"batch_all"));
        else
            OSNAP_data.classification_summary_batch_all = evaluate_OSNAP_classifiers(OSNAP_data.classifiers_batch_allOSNAP_data.classifiers_batch_all);
        end
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
end
if isfield(OSNAP_data,"classification_summary_batch_each")
    fprintf('    Classification summary (batch - fully independent feature selection):\n');
    head(OSNAP_data.classification_summary_batch_each,3)
end
if isfield(OSNAP_data,"classification_summary_batch_all")
    fprintf('    Classification summary (batch - aggregated feature selection):\n');
    head(OSNAP_data.classification_summary_batch_all,3)
end
%% Run classification on whole data
try
    if options.run_classification_all
        n_models_per_type = options.n_models_per_type;
        fprintf("  Running classification (whole dataset)...\n")
        starttime_step = tic;
        % cell_dims = [numel(OSNAP_data.train_idxs),numel(OSNAP_data.vars_selected_batch),numel(OSNAP_data.pca_result)];
        % assert(all(cell_dims == cell_dims(1)),'OSNAP:classification_length_mismatch','number of batches not consistent with classification structures (feature selection, PCA)')
        % Split batches for parallel
        % Variable selection
        save_path = fullfile(work_dir, join(['feature_selection_all',groups],'_'));
        if options.suffix ~= ""
            save_path = string(join([save_path options.suffix],'_'));
        end
        [~, OSNAP_data.vars_select_result_all] = fscmrmr(feature_data_filtered(:,vartype('numeric')),feature_data_filtered.group);
        [~,OSNAP_data.vars_selected_all] = select_OSNAP_features(OSNAP_data.vars_select_result_all',feature_data_filtered(:,vartype('numeric')).Properties.VariableNames,"max_idx",options.feature_select_max_idx,"save_path",save_path);
        fprintf("      Selected features - ALL (n = %.0f):\n",numel(OSNAP_data.vars_selected_all))
        for i=1:numel(OSNAP_data.vars_selected_all)
            fprintf("        - %s\n", OSNAP_data.vars_selected_all{i})
        end
        % PCA
        save_path = fullfile(work_dir, string(join(['PCA',groups],'_')));
        if options.suffix ~= ""
            save_path = string(join([save_path options.suffix],'_'));
        end
        OSNAP_data.pca_result_all = run_OSNAP_PCA(feature_data_filtered,...
            "vars_sel",OSNAP_data.vars_selected_all,...
            "save_path",save_path,...
            "num_components_explained",options.num_components_explained);
        % Classification
        classifiers =...
            run_OSNAP_classification_batch(...
            feature_data_filtered,...
            'n_models_per_type',n_models_per_type,...
            'vars_selected',OSNAP_data.vars_selected_all,...
            'pca_result',OSNAP_data.pca_result_all,....
            'verbose',0);
        OSNAP_data.classifiers_all = {classifiers};
        if options.save
            OSNAP_data.classification_summary_all = evaluate_OSNAP_classifiers(OSNAP_data.classifiers_all,"save_path",fullfile(work_dir,"all"));
        else
            OSNAP_data.classification_summary_all = evaluate_OSNAP_classifiers(OSNAP_data.classifiers_allOSNAP_data.classifiers_all);
        end
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
end
if isfield(OSNAP_data,"classification_summary_all")
    fprintf('    Classification summary (batch - aggregated feature selection):\n');
    head(OSNAP_data.classification_summary_all,3)
end
%% Plot radial
try
    if options.run_plot_radial
        fprintf("  Plotting radial densities...\n")
        plot_OSNAP_ellipses_batch(feature_data_filtered,work_dir,groups,replicates);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
%% FSEA
try
    if options.run_FSEA
        fprintf("  Running feature set enrichment analysis...\n")
        starttime_step = tic;
        OSNAP_data.feature_set_coverage = calculate_OSNAP_feature_set_coverage(OSNAP_data.feature_comparisons,work_dir,"plot",false,"feature_universe_names",options.feature_universe_names);
        run_OSNAP_FSEA(work_dir, feature_data_filtered, options.feature_universe_names, "FSEA_rank_type",options.FSEA_rank_type);
        fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
    end
catch ME
    handle_OSNAP_error(ME,save_analysis_path,OSNAP_data,"save",options.save_if_error);
    return
end
%% Finish run
OSNAP_data.options = options;
conclude_OSNAP_run(save_analysis_path,OSNAP_data,"save",options.save)
end

%% Performs wrap up tasks to end the OSNAP run
function conclude_OSNAP_run(save_path,data,options)
arguments
    save_path string
    data struct
    options.save logical = true
end
if options.save
    fprintf("  Saving run to: %s...\n",save_path)
    starttime_step = tic;
    save_OSNAP_run(save_path,data,"append",false)
    fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
end
runtime = toc(data.starttime)/60;
fprintf("RUN END: %s\n", string(datetime))
fprintf("TIME ELAPSED: %.2f min\n", runtime)
fprintf("Done!\n")
fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
fclose('all');
warning('on','all');
set(findall(groot,'Type','Figure'),'visible','on');
diary off
return
end
%% Handles OSNAP errors
function handle_OSNAP_error(ME,filepath,data,options)
arguments
    ME
    filepath string
    data struct
    options.stop_run logical = true
    options.save logical = true
end
if nargin == 3
    fprintf("- - - - - - - - - - - OSNAP ERROR - - - - - - - - - - -\n")
    fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
    fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
elseif nargin == 4
    fprintf("- - - - - - - - - - - OSNAP ERROR - - - - - - - - - - -\n")
    fprintf("%s\n",getReport(ME,'extended','hyperlinks','off'))
    fprintf("- - - - - - - - - - - - - - - - - - - - - - - - - - - - \n")
end
if options.save
    fprintf("  Saving run to: %s...\n",filepath)
    starttime_step = tic;
    save_OSNAP_run(filepath,data,"append",false)
    fprintf("      Completed %s (%.2f min)...\n",string(datetime),toc(starttime_step)/60);
end
if options.stop_run
    conclude_OSNAP_run(filepath,data,"save",false);
end
end
%% Save OSNAP run
function save_OSNAP_run(save_path,data,options)
arguments
    save_path string
    data struct
    options.check_overwrite logical = 0
    options.append logical = 1
end
% If structure is empty, return
if isempty(fieldnames(data))
    return
end
% Append option
if options.append
    if exist(save_path,'file')
        save(save_path,'-struct','data','-append');
    else
        fprintf("Warning: File does not exist, creating file: "+save_path)
        save(save_path,'-struct','data');
    end
    % Overwrite file
elseif options.check_overwrite
    if exist(save_path,'file')
        overwrite_dialog();
    else
        save(save_path,'-struct','data');
    end
    % Overwrite file with no check
else
    save(save_path,'-struct','data');
end
    function overwrite_dialog()
        response = questdlg(sprintf('Warning: %s\nalready exists. Overwrite?',save_path),...
            'Warning: Overwrite',...
            'Overwrite','New file...','Cancel','Cancel');
        switch response
            case 'Overwrite'
                fprintf('      Overwriting...')
                save(save_path,'-struct','data');
            case 'New file...'
                new_path = uigetfile();
                save(new_path,'-struct','data');
            case 'Cancel'
                fprintf('Save canceled')
        end
    end
end
