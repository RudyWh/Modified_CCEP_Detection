function out = AA_grand_avg(out)

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

NumRows = size(out.elecs,2);
for ich=1:NumRows
    ich;
    for jch=1:NumRows
        c = ((ich-1)*(NumRows))+jch;
        c;
        if size(out.elecs(ich).avg,1)>=1 && ~isnan(out.AnalysisData_N1(c,2)) &&   (out.AnalysisData_N1(c,2))~=0
            N1_chosen_idx = round(out.AnalysisData_N1(c,2)*(out.other.stim.fs)) + stim_start;
            N1_signed_value = out.elecs(ich).avg(N1_chosen_idx,jch);
            sign_itself = sign(N1_signed_value);
            out.AnalysisData_N1(c,7) = sign_itself;
            out.AnalysisData_N2(c,7) = sign_itself;
        end
    end
end
%% 
target_size = size(out.elecs(chosen_row).avg);
NumRows = size(out.elecs,2);
for ich=1:NumRows
    ich;
    if size(out.elecs(ich).avg,1)>=1
        out.elecs(ich).detrend_flipped = NaN(target_size);
        for jch=1:NumRows
            c = ((ich-1)*(NumRows))+jch;
            stim_channel = out.AnalysisData_N1(c,5);
            response_channel = out.AnalysisData_N1(c,6);
            orientation_val = out.AnalysisData_N1(c,7);
            if orientation_val == 1
                % out.elecs(stim_channel).avg(:,response_channel)
                out.elecs(stim_channel).detrend_flipped(:,response_channel) = out.elecs(stim_channel).avg(:,response_channel);
            end
            if orientation_val == -1
                out.elecs(stim_channel).detrend_flipped(:,response_channel) = -1*(out.elecs(stim_channel).avg(:,response_channel));
            end
        end
    end
end
%%
% Take SOZ channels out 
target_length = size(out.elecs(chosen_row).avg,1);
SOZ_electrodes_pt1 = [];
for ich = 1:NumRows
    
    for jch = 1:NumRows
        if ismember(ich,SOZ_electrodes)
            out.elecs(ich).detrend_flipped(:,jch) = NaN(target_length);
        end
    end
end
%% 
for ich = 1:NumRows
    if size(out.elecs(ich).detrend_flipped,1)>1
        out.elecs(ich).avg_weight = 0;
    end
end
%% 


NumRows = size(out.elecs,2);
for ich = 1:NumRows
    for jch = 1:NumRows
        if out.rejection_details(1).reject.keep(ich,jch) ~=1
            out.elecs(ich).detrend_flipped(:,jch) = nan;
        end
        if size(out.elecs(ich).deriv,1) <=1
            out.elecs(ich).detrend_flipped = [];
        end
        if size(out.elecs(ich).deriv,1) > 1
            if ~isnan(out.elecs(ich).detrend_flipped(1,jch))
                out.elecs(ich).avg_weight = out.elecs(ich).avg_weight + 1;
            end
        end
    end
end
%% 
% NumRows = size(out.elecs,2);
% for ich = 1:NumRows
%     for jch = 1:NumRows
%         if strcmp('temporal',out.chLabels{jch,9})
%         %if strcmp('m_temporal',out.chLabels{jch,9})
%         %if strcmp('insula',out.chLabels{jch,9})
%         %if strcmp('parietal',out.chLabels{jch,9})
%         %if strcmp('frontal',out.chLabels{jch,9})
%         %if strcmp('BG',out.chLabels{jch,9})
%             ok = 1;
%         else
%             out.elecs(ich).detrend_flipped(:, jch) = NaN((size(out.elecs(ich).detrend_flipped,1)),1);
%         end
%     end
% end

%% GETTING AVGERAGE FOR ALL RESPONSES WHEN A CHANNEL STIMMED 
NumRows = size(out.elecs,2);
for ich=1:NumRows
    ich;
    num_points = size(out.elecs(ich).detrend_flipped,1);
    for x = 1:num_points
        x_mean = nanmean(out.elecs(ich).detrend_flipped(x,:));
        out.elecs(ich).channel_avg(x,1) = x_mean;
    end
end
%%
chosen_row = 5;
NumRows = size(out.elecs,2);
row_idx = size(out.elecs(chosen_row).detrend_flipped,1)
column_idx = size(out.elecs(chosen_row).detrend_flipped,2)

out.all_ch_avgs = NaN(row_idx, column_idx);
%%
sum_total_weight = 0;
for ich=1:NumRows
    if size(out.elecs(ich).channel_avg,1)>=1
        col_vect = out.elecs(ich).channel_avg(:,1);
        out.all_ch_avgs(:,ich) = col_vect*(out.elecs(ich).avg_weight);
        sum_total_weight = sum_total_weight + out.elecs(ich).avg_weight
    end
end
%% 
NumRows = size(out.all_ch_avgs,1);
% for c=1:NumRows
   x_sum = nansum(out.all_ch_avgs,2);
   out.all_ch_avgs(:,(size(out.elecs,2))+1) = x_sum/sum_total_weight;
% end

%% THIS PART MAY BREAK FOR ANYONE BUT 29 CURRENTLY
% avg = out.elecs(10).avg(:,1);
% times = out.elecs(10).times;
% eeg_times = convert_indices_to_times(1:length(avg),out.other.stim.fs,times(1));
%% 
NumRows = size(out.elecs,2);
for ich=1:NumRows
    for jch = 1:NumRows
        if size(out.elecs(ich).detrend_flipped)>=1 
            if ~isnan(out.elecs(ich).detrend_flipped(1,jch))
                % plot(eeg_times,out.elecs(ich).detrend_flipped(:,jch), 'black')
                plot(out.elecs(ich).detrend_flipped(:,jch),'black')
                hold on
            end
        end
    end
end
%% 

NumRows = size(out.elecs,2);
out.dist_sum = 0;
out.kept_counter = 0;
for ich = 1:NumRows
    ich
    for jch = 1:NumRows
        if size(out.elecs(ich).detrend_flipped,1)>=2 
            if ~isnan(out.elecs(ich).detrend_flipped(1,jch)) && ~isnan(out.elecs(ich).N1(jch,3))
                out.dist_sum = out.dist_sum + out.elecs(ich).N1(jch,3);
                out.kept_counter = out.kept_counter + 1;
            end
        end
    end
end
%% 
% xline((out.other.stim.fs)/2, 'red')

avg_distance = out.dist_sum / out.kept_counter;

plot(out.all_ch_avgs(:,NumRows+1),'red','LineWidth',2)
% plot(eeg_times,out.all_ch_avgs(:,(NumRows+1)),'red', 'LineWidth',2)
smallFont = {'fontsize',8};
% title('CHOPCCEP\_029 Responses')
xlabel('Time (s)',smallFont{:}) 
ylabel('Amplitude (mA)',smallFont{:})
aaa = max(out.all_ch_avgs(:,(NumRows+1)));
% xlim([-0.5 0.8])
%ylim([-25 45])
N1_dat = out.all_ch_avgs([((1020/2048)*(out.other.stim.fs)):((1139/2048)*(out.other.stim.fs))],(NumRows+1));
[N1_peaks,N1_indicies] = findpeaks(N1_dat);
[N1_MaxPeak,N1_Index] = max(N1_peaks);
N1_MaxPeak
N1_result = find((out.all_ch_avgs(:,(NumRows+1)))==N1_MaxPeak);
final_N1 = (N1_result - stim_start)/(out.other.stim.fs)
% final_N1_adjusted = final_N1 / avg_distance
% N1_MaxPeak_adjusted = N1_MaxPeak / avg_distance

N2_dat = out.all_ch_avgs([((1140/2048)*(out.other.stim.fs)):((1840/2048)*(out.other.stim.fs))],(NumRows+1));
[N2_peaks,N2_indicies] = findpeaks(N2_dat);
[N2_MaxPeak,N2_Index] = max(N2_peaks);
N2_MaxPeak
N2_result = find((out.all_ch_avgs(:,(NumRows+1)))==N2_MaxPeak);
final_N2 = (N2_result - stim_start)/(out.other.stim.fs)
% final_N2_adjusted = final_N2 / avg_distance
% N2_MaxPeak_adjusted = N2_MaxPeak / avg_distance

avg_distance
keep_counter = out.kept_counter
end
