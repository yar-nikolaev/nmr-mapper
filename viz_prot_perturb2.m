function viz_prot_perturb2()

% TODO:
% - Automatically? find out which expnos used for the mapper?!? (Is this possible?)

% Versions:
% viz_prot_perturb:   initial version, based on 
%   Projects/2015_NMR_peak_tracking/site_mapper/190411_TDP_FUS_SMN2_LLPS/analysis/runall_TDP_190411.m

% viz_prot_perturb2:
%   - includes 31P data import and plotting against it
%   - partially based on CUG / MBNL analysis /Volumes/Data/yar/Dropbox/_eth2/project_CUG_MBNL/code/180926_mapper_with_31P/analysis/runall_var_MBNL_31P_v180930_f_report.m
%   - option to plot both CSP and Intensity
%   - option to plot chemical shifts of peaks
%   - fixed opt_to_duplicate

%% Settings
%=================================
% Sections to run:
flag_run_assignments = 0;
flag_load_assignments = 1;

% Switching Y-data
flag_plot_both_intensity_and_CSP_on_one_plot = 1; % This forces plotting of the first dataset!

flag_plot_intensity_instead_of_traj = 0;
if flag_plot_intensity_instead_of_traj
    optns.yaxis_label = 'peak amplitude';
    optns.normalize_intensity = 1; % is only set if plotting intensity
else
    optns.yaxis_label = 'CSP';
end;
optns.same_yscale = 1;

flag_show_peak_CS = 1;

% Switching X-Axis data
xaxis = 3; % 1 - spectrum id; 2 - time; 3 - 31P(RNA)

%%% Selecting subsets (if have a list of sets defined, but want to plot
%%% only a few):
% select_dsets = [1:2 4]; % 1:2 4 - compares 3rd nucleotide
% select_dsets = [1 3 5]; % 1:2 4 - compares 4th nucleotide

% select_residues = [39 5 19 27 30 35 65 5 19 27 30 35 65 5 19 27 30 35 65 5 19 27 30 35];
% select_residues = [1:48 50:61];
% select_residues = [1:5];

% "select_spectra" - in reverse - QnD solution to fix nc_proc
% drop_timepoints = [2 3 6];
% drop_timepoints = [2 3 6 22:35]; % tmp hack - to match MAPPER vector to P50N525 time

optns.colors = [...
        [0, 0, 0]; % black
        [0.7, 0, 0]; % dark red
        [0, 0.6, 0]; % dark green
        [1, 0, 0]; % red
        [0, 0.9, 0]; % green        
        ];
    
%============================================================
%%% Inputs: mapper data, peaklists, RNA integration, protein sequences
%============================================================
current_dir = fileparts(mfilename('fullpath'));
nmr_dataset_dir = fullfile(current_dir,'..');

dset_names = {...
%     'XXX'
    'P50N525 A20'
    };

%%% Mapper data
%============
mapper_data_dirs = {
%     'XXX'
%     fullfile(fileparts(mfilename('fullpath')), '..', 'mapper_data')
    fullfile(fileparts(mfilename('fullpath')), '..', 'IN115a_mapper')
    };
% s3_tracking dir - also contains intensities
mapper_data_paths = cellfun(@(x) fullfile(x,'s3_tracking','trajectories_filtered.csv'), mapper_data_dirs, 'un', 0);

% Need to provide these manually at the moment to tell where to read the
% times of experiments from.
mapper_expnos = [4003 4008 4013 4018 4023 4028 4033 4038 4043 4078 4083 4088];

n_sets = numel(mapper_data_paths);

%%% Peaklists (can provide multiple, if have multiple datasets)
%============
peaklist_dirs = {...
%     fullfile( current_dir, '..', 'peaklists' )
    fullfile( current_dir, '..', 'peaklists' )
    };
peaklist_names = {...
%     'TDPNTD_peaklist.peaks'
    'peaks_start_mapper.peaks'
};
peaklist_paths = cellfun(@(x,y) fullfile(x,y), peaklist_dirs, peaklist_names, 'un', 0);


%%% RNA data (optional, if not trying to plot time/RNA on xaxis)
%============
RNA_data_dirs = {...
    'datasave'
    };
RNA_data_filenames = {...
    'RNA_A20_P50N525.csv'
    };
RNA_data_paths = cellfun(@(x,y) fullfile(x,y), RNA_data_dirs, RNA_data_filenames, 'un', 0);


%%% Protein sequences (optional)
%============
% prot_sequence_dirs = repmat({'data_test'}, n_sets, 1);
% prot_seq_filenames = repmat({'UP1_1-196_1-letter.aa'}, n_sets, 1);
% prot_seq_paths = cellfun(@(x,y) fullfile(x,y), prot_sequence_dirs, prot_seq_filenames, 'un', 0);
prot_seq_paths = repmat({''}, n_sets, 1);


%================================================================
%%%% EDITING BELOW THIS LINE IS NORMALLY NOT REQUIRED
%================================================================

%% Select only some datasets
%================
if exist('select_dsets','var') && ~isempty(select_dsets)
    dset_names = dset_names(select_dsets);
    mapper_data_paths = mapper_data_paths(select_dsets);
    peaklist_paths = peaklist_paths(select_dsets);
    RNA_data_paths = RNA_data_paths(select_dsets);
    n_sets = numel(mapper_data_paths);
end

%% Import RNA conc and time. Interpolate RNA conc to mapper 2DHN times.
%======================================================================
if xaxis > 1 
    f_rna_data_exists = (exist('RNA_data_paths','var') && ~isempty(RNA_data_paths));
    assert(f_rna_data_exists, 'Abort: file with RNA & time data seems missing.');    
    
    tmp = cellfun(@(x) mapper.import_31P_n_time_from_csv(x), RNA_data_paths, 'un', 0);
    RNA_time = cellfun(@(x) x(:,1), tmp, 'un', 0);
    RNA = cellfun(@(x) x(:,2), tmp, 'un', 0);
    clear tmp;
    
    if exist('mapper_expnos','var') && ~isempty(mapper_expnos)
        time0 = getTime0(nmr_dataset_dir, '.');
        HN_time = {getNMRTime(nmr_dataset_dir,mapper_expnos,time0)};
    else
        fprintf('Using time or RNA in xaxis requires mapper_expnos for interpolation.');
    end;
    
    RNA_in_HN_time = cellfun(@(x,y,z) interp1(x,y,z, 'linear', 'extrap'), RNA_time, RNA, HN_time, 'un', 0);
    
end;

%% Data to use for x-axis
%================
switch xaxis
    case 1
        optns.xaxis_label = 'spectrum id';
    case 2
        optns.xaxis_label = 'time [min]';        
        optns.xaxis = HN_time; % Needs to be a CELL ARRAY (if have many datasets!)
    case 3
        optns.xaxis_label = 'RNA conc [mM]';
        optns.xaxis = RNA_in_HN_time;
end;

%% Assign trajectories to peaks
%=================================
if flag_run_assignments
    assign = cell(n_sets,1);
    for i=1:n_sets                        
        assign{i} = mapper.match_assign_with_traj_v04_intensity(dset_names{i},...
            peaklist_paths{i},...
            mapper_data_paths{i},...
            prot_seq_paths{i}...
            );
    end    
end

if flag_load_assignments
    assign = cell(n_sets,1);
    for i=1:n_sets
        tmp = load( fullfile('datasave', sprintf('assign_%s.mat', dset_names{i})) );
        assign{i} = tmp.assign;
    end
    clear tmp;
end
        

%% Visualization
%=================================
optns.plotLW = 2;
% SMN1 vs SMN2 for poster
% + need to use plot_traj3_rb - for color!
% optns.plotSym = 'o-';
% optns.markerSize = 3;

optns.plotSym = '-';
optns.markerSize = 4;

optns.same_min_xscale = 1;
optns.legend = dset_names;

% arrayfun(@(x) find([8 9]==x,1), assign{1}.ids, 'unif', 0)
% [tf, idx] = ismember([8,12], assign{1}.ids)    

% seq_end = numel(assign{1}.prot_seq);
% peak_names = assign{1}.prot_seq(assign{1}.ids(assign{1}.ids < seq_end));
% optns.titles = peak_names;
optns.titles = assign{1}.names;

if ~flag_plot_intensity_instead_of_traj
    traj_select = cellfun(@(x) x.traj_HNcsp, assign, 'unif', 0);
else
    traj_select = cellfun(@(x) x.intens, assign, 'unif', 0);
end

if flag_plot_both_intensity_and_CSP_on_one_plot
    traj_select{1} = assign{1}.traj_HNcsp;
%     traj_select{2} = assign{1}.traj_HNcsp;
    traj_select{2} = assign{1}.intens;    

    n_sets = 2;    
    optns.legend = cellfun(@(x) strcat(optns.legend{1}, x), {': CSP', ': int'}, 'un', 0)';
    
    % Scale intensities such that they'd show normalized:
    traj_select{2} = bsxfun(@rdivide, traj_select{2}, max(traj_select{2},[],1));
    traj_select{2} = traj_select{2} .* max(max(traj_select{1},[],2)); % max value of all CSPs.    
end

%%% Replicate plotting options -- TODO: remove? (don't really need to
%%% define this here?    
opt_to_duplicate = {'xaxis', 'plotLW', 'plotSym'};
for iOpt=1:numel(opt_to_duplicate)
    if ~iscell(optns.(opt_to_duplicate{iOpt})) % encapsulate into cell if not a cell yet.
        optns.(opt_to_duplicate{iOpt}) = {optns.(opt_to_duplicate{iOpt})};
    end;            
    optns.(opt_to_duplicate{iOpt}) = repmat(optns.(opt_to_duplicate{iOpt}),n_sets,1);        
end;

% Select only subset of residues -- works only for TRAJECTORIES at the
% moment!
if exist('select_residues','var') && ~isempty(select_residues)
    %%% Use IDS if the peaklist had them defined
    %%% otherwise just peak numbers "as is"
    if isfield(assign{1},'ids')
        [~, idx] = ismember(select_residues, assign{1}.ids);
        optns.titles = peak_names( idx );
    else
        idx = select_residues;
        optns.titles = assign{1}.names( idx );    
    end;
    
    if flag_show_peak_CS
        optns.titles = arrayfun(@(name,H,N) {name{:}, sprintf('%.2f-%.1f', H, N)}, optns.titles, assign{1}.traj_H(idx,1), assign{1}.traj_N(idx,1), 'un', 0);
    end;
    
    if flag_plot_both_intensity_and_CSP_on_one_plot
        traj_select = cellfun(@(x) x(select_residues,:), traj_select, 'unif', 0);        
    else
        traj_select = cellfun(@(x) x.traj_HNcsp(idx,:), assign, 'unif', 0);        
    end    
    
else    
    if flag_show_peak_CS
        optns.titles = arrayfun(@(name,H,N) {name{:}, sprintf('%.2f-%.1f', H, N)}, optns.titles, assign{1}.traj_H(:,1), assign{1}.traj_N(:,1), 'un', 0);
    end;
end

if exist('drop_timepoints','var') && ~isempty(drop_timepoints)
    % TODO - vectorize this
    for iSet = 1:n_sets
        traj_select{iSet}(:,drop_timepoints) = [];
    end;
end

% fig_handle = plot_traj3(traj_select, optns);
% fig_handle = mapper.plot_traj3b(traj_select, optns);
fig_handle = mapper.plot_traj4b(traj_select, optns);

base_name = [sprintf('%s_',dset_names{1:end-1}),dset_names{end}];
if ~flag_plot_intensity_instead_of_traj
    mapper.save_figure(fig_handle, base_name);
else
    mapper.save_figure(fig_handle, sprintf('%s_intens', base_name));
end

end