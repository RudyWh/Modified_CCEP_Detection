% addpath 'H:\Rudy\CHOP\CCEPS-main\CCEPS-main'
% addpath 'H:\Rudy\CHOP\CCEPS-main\CCEPS-main\do_run\calculating_cceps_and_network'
% addpath 'H:\Rudy\CHOP\CCEPS-main\CCEPS-main\do_run\eeg_processing'
% addpath 'H:\Rudy\CHOP\CCEPS-main\RW_files'
% out = BB_filtered_sweep_analysis(out);
% addpath '\\ressmb05.research.chop.edu\marsh_lab2\Caren\CCEP_files-CCEP_code\CCEPS-main\do_run\calculating_cceps_and_network'
out = AA_alternative_filtering(out);
out = AA_Running_RejectOrKeep_RW(out);
out = AA_new_build_network(out,0);
out = AA_distance_vs_amp_and_lat(out);
out = AA_require_both_Ns(out);
out = AA_modify_analysisdata(out);
% AA_random_rejections_keeps(out);
%%
addpath '\\ressmb05.research.chop.edu\marsh_lab2\Caren\CCEP_files-CCEP_code'
% BP_adjusted_results = cell(46,6);
%%
disp('Select the appropriate .mat files')
[files,path] = uigetfile('.mat','MultiSelect','on');
if ischar(files)
    nofiles = 1;
else
    nofiles = size(files,2);
end

SOZ_restriction_control = 0;
saved_info_cell = cell(46,1);

for file_num=1:nofiles
    file_num
    if nofiles == 1
        out_nested = load([path files]);
    else
        out_nested = load([path files{file_num}]); %name the out file
    end
    out = out_nested.out;
    % [real_file, out] = RP_get_filtered_file(file_num);
    % if real_file == 1 
        out = AA_alternative_filtering(out);
        out = AA_Running_RejectOrKeep_RW(out);
        out = AA_new_build_network(out,0);
        out = AA_distance_vs_amp_and_lat(out);
        out = AA_require_both_Ns(out);
        out = AA_modify_analysisdata(out);

        % % Reject White Matter
        % for ich = 1:size(out.chLabels,1)
        %     for jch = 1:size(out.chLabels,1)
        %         if out.chLabels{ich,6}~=1
        %             out.rejection_details(1).reject.keep(ich,jch) = 0;
        %             out.rejection_details(2).reject.keep(ich,jch) = 0;
        %             out.rejection_details(1).reject.not_grey(ich,jch) = 1;
        %             out.rejection_details(2).reject.not_grey(ich,jch) = 1;
        %         end
        %         if out.chLabels{jch,6}~=1
        %             out.rejection_details(1).reject.keep(ich,jch) = 0;
        %             out.rejection_details(2).reject.keep(ich,jch) = 0;
        %             out.rejection_details(1).reject.not_grey(ich,jch) = 1;
        %             out.rejection_details(2).reject.not_grey(ich,jch) = 1;
        %         end
        %     end
        % end

        %%% out = reject_SOZ_stim(out); %% COMMENT BACK OUT LATER (rejects all SOZ channels)
        
        if SOZ_restriction_control==1
            out = reject_NON_SOZ_stim(out);
        end
        if SOZ_restriction_control==2
            out = reject_SOZ_resp(out);
        end
        if SOZ_restriction_control==3
            out = reject_NON_SOZ_resp(out);
        end
        if SOZ_restriction_control==4
            out = reject_SOZ_both(out);
        end
        if SOZ_restriction_control==5
            out = reject_NON_SOZ_both(out);
        end
        if SOZ_restriction_control==6
            out = reject_SOZ_stim(out);
        end
        
        % NumRows = size(out.elecs,2);
        % for ich = 1:NumRows
        %     if isempty(out.elecs(ich).arts)
        %         ppp = 9;
        %         out.elecs(ich).N1 = [];
        %         out.elecs(ich).N2 = [];
        %     end
        % end
        
        num_sum = 0;
        dist_sum = 0;

        N1_amp_sum = 0;
        N2_amp_sum = 0;
        N1_lat_sum = 0;
        N2_lat_sum = 0;

        saved_N1_amp_counter = 1;
        saved_N2_amp_counter = 1;
        saved_N1_lat_counter = 1;
        saved_N2_lat_counter = 1;

        saved_info = zeros(2,4);

        for ich = 1:size(out.elecs,2)
            for jch = 1:size(out.elecs,2)
                if out.rejection_details(1).reject.keep(ich,jch)==1 && out.rejection_details(2).reject.keep(ich,jch)==1 && ~isempty(out.chLabels{ich,2}) && ~isempty(out.chLabels{ich,3}) && ~isempty(out.chLabels{ich,4}) && ~isempty(out.chLabels{jch,2}) && ~isempty(out.chLabels{jch,3}) && ~isempty(out.chLabels{jch,4})

                    stim_x = out.chLabels{ich,2};
                    stim_y = out.chLabels{ich,3};
                    stim_z = out.chLabels{ich,4};
                    resp_x = out.chLabels{jch,2};
                    resp_y = out.chLabels{jch,3};
                    resp_z = out.chLabels{jch,4};
                    val = sqrt(((stim_x - resp_x)^2)+((stim_y - resp_y)^2)+((stim_z - resp_z)^2));

                    dist_sum = dist_sum + val;
                    num_sum = num_sum + 1;

                    N1_indx = out.elecs(ich).N1(jch,2);
                    N2_indx = out.elecs(ich).N2(jch,2);

                    N1_amp_sum = N1_amp_sum + abs(out.elecs(ich).avg(N1_indx,jch));
                    N2_amp_sum = N2_amp_sum + abs(out.elecs(ich).avg(N2_indx,jch));
                    N1_lat_sum = N1_lat_sum + (N1_indx/out.other.stim.fs);
                    N2_lat_sum = N2_lat_sum + (N2_indx/out.other.stim.fs);

                    saved_info(saved_N1_amp_counter,1) = abs(out.elecs(ich).avg(N1_indx,jch));
                    saved_info(saved_N2_amp_counter,2) = abs(out.elecs(ich).avg(N2_indx,jch));
                    saved_info(saved_N1_lat_counter,3) = (N1_indx/out.other.stim.fs);
                    saved_info(saved_N2_lat_counter,4) = (N2_indx/out.other.stim.fs);

                    saved_N1_amp_counter = saved_N1_amp_counter + 1;
                    saved_N2_amp_counter = saved_N2_amp_counter + 1;
                    saved_N1_lat_counter = saved_N1_lat_counter + 1;
                    saved_N2_lat_counter = saved_N2_lat_counter + 1;
                end
            end
        end
        % % Not distance adjusted
        % BP_adjusted_results{file_num,1} = dist_sum/num_sum;
        % BP_adjusted_results{file_num,2} = N1_amp_sum/num_sum;
        % BP_adjusted_results{file_num,3} = N2_amp_sum/num_sum;
        % BP_adjusted_results{file_num,4} = N1_lat_sum/num_sum;
        % BP_adjusted_results{file_num,5} = N2_lat_sum/num_sum;
        % BP_adjusted_results{file_num,6} = num_sum;

        % Distance adjusted
        BP_adjusted_results{file_num,1} = dist_sum/(num_sum);
        BP_adjusted_results{file_num,2} = N1_amp_sum/(num_sum*BP_adjusted_results{file_num,1});
        BP_adjusted_results{file_num,3} = N2_amp_sum/(num_sum*BP_adjusted_results{file_num,1});
        BP_adjusted_results{file_num,4} = N1_lat_sum/(num_sum*BP_adjusted_results{file_num,1});
        BP_adjusted_results{file_num,5} = N2_lat_sum/(num_sum*BP_adjusted_results{file_num,1});
        BP_adjusted_results{file_num,6} = num_sum;

        % saved_info_cell{file_num,1} = saved_info;
        % name_str = string(file_num);
        % save(name_str,'out',"-v7.3")
    % end
end















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
