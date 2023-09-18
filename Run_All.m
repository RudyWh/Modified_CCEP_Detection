addpath 'H:\Rudy\CHOP\CCEPS-main\CCEPS-main'
addpath 'H:\Rudy\CHOP\CCEPS-main\CCEPS-main\do_run\calculating_cceps_and_network'
addpath 'H:\Rudy\CHOP\CCEPS-main\CCEPS-main\do_run\eeg_processing'
addpath 'H:\Rudy\CHOP\CCEPS-main\RW_files'
% out = BB_filtered_sweep_analysis(out);
out = AA_alternative_filtering(out);
out = AA_Running_RejectOrKeep_RW(out);
out = AA_new_build_network(out,0);
out = AA_distance_vs_amp_and_lat(out);
out = AA_require_both_Ns(out);
out = AA_modify_analysisdata(out);
% AA_random_rejections_keeps(out);
%%
out = AA_grand_avg(out);
%%
out = AA_not_GA(out);
%%
out = RW_Area_Magnitude(out)
%% 

%plot(crp_parms.C)
out.elecs = rmfield(out.elecs,'detrend_flipped');
out.elecs = rmfield(out.elecs,'deriv');
out.elecs = rmfield(out.elecs,'avg_weight');
out.elecs = rmfield(out.elecs,'channel_avg');
out.elecs = rmfield(out.elecs,'detrend_filt_avgs');

out = rmfield(out,'AnalysisData_N1');
out = rmfield(out,'AnalysisData_N2');
%%
num_channels = size(out.chLabels,1);
% Finds first channel with keeps to reference for stim index, etc. 
kept_channel_found = 0;
channel_check = 1;
while kept_channel_found == 0 
    if ~isempty(out.elecs(channel_check).stim_idx)
        interpolation_period_start = out.elecs(channel_check).stim_idx - 3;
        kept_channel = channel_check; 
        kept_channel_found = 1;
    end
    channel_check = channel_check + 1;
end

out.CRP = struct([]);
out.CRP(1).data = cell(num_channels,num_channels);
out.CRP(1).time = cell(num_channels,num_channels);
out.CRP(1).C = cell(num_channels,num_channels);
out.CRP(1).p_val_tr = cell(num_channels,num_channels);

for ich = 1:num_channels
    for jch = 1:num_channels
        if ~isempty(out.elecs(ich).raw_sweeps)
            data_temp = out.elecs(ich).raw_sweeps(:,jch,:);
            data = reshape(data_temp, [size(data_temp,1) size(data_temp,3)]);
            time_points = size(out.elecs(kept_channel).avg,1);
            start_time = out.elecs(kept_channel).times(1);
            end_time = out.elecs(kept_channel).times(2);
            t = linspace(start_time,end_time,time_points);

            out.CRP.data{ich,jch} = data;
            out.CRP.time{ich,jch} = t;
        end
    end
end
%% 
for ich = 1:num_channels
    ich
    for jch = 1:num_channels
        V = out.CRP.data{ich,jch};
        t_win = out.CRP.time{ich,jch};
        if ~isempty(V) && ~isempty(t_win) 
            if ~isnan(out.CRP.data{ich,jch})
                [crp_parms, crp_projs]=CRP_method(V,t_win);
                out.CRP.parms{ich,jch}.C = crp_parms.C;
                out.CRP.projs{ich,jch}.p_val_tr = crp_projs.p_value_tR;
            end
        end
    end
end
