%% To run modified detection, open an existing out file & run file - add path directories may need to be filled in with where files from the original pipeline are in your structure/commented in 

% addpath 'YOUR_DIRECTORY_OF\CCEPS-main'
% addpath 'YOUR_DIRECTORY_OF\CCEPS-main\do_run\calculating_cceps_and_network'
% addpath 'YOUR_DIRECTORY_OF\CCEPS-main\do_run\eeg_processing'

%% RUN SINGLE OUT FILE 
out = AA_alternative_filtering(out);
out = AA_Running_RejectOrKeep_RW(out);
out = AA_new_build_network(out,0);
out = AA_distance_vs_amp_and_lat(out);
out = AA_require_both_Ns(out);
out = AA_modify_analysisdata(out);

% Can be commented in to display 25 random keeps/rejects as detected by the pipeline
% AA_random_rejections_keeps(out);

%% RUN MULTIPLE OUT FILES & COMPILE INFO ABOUT AMPLITUDE/LATENCY
% (requires additional information about electrode coordinate location & clinically identified seizure onset zone contacts to be added to out structure) 
% BP_adjusted_results = cell(46,6);
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
