function out = AA_not_GA(out)

NumRows = size(out.elecs,2);
N1_number = 0;
N1_amp_sum = 0;
N1_lat_sum = 0;
N1_dist_sum = 0;
baseline_sum = 0;
for ich = 1:NumRows
    for jch = 1:NumRows
        if size(out.elecs(ich).N1,1)>=2
            if out.rejection_details(1).reject.keep(ich,jch) ==1 && out.rejection_details(2).reject.keep(ich,jch) ==1
                current_N1_amp = out.elecs(ich).N1(jch,1);
                current_N1_lat = out.elecs(ich).N1(jch,2);
                current_dist = out.elecs(ich).N1(jch,3);
                current_baseline_sd = out.elecs(ich).N1(jch,4);
    
                N1_number = N1_number + 1;
                N1_amp_sum = N1_amp_sum + current_N1_amp;
                N1_lat_sum = N1_lat_sum + current_N1_lat;
                N1_dist_sum = N1_dist_sum + current_dist;
                baseline_sum = baseline_sum + current_baseline_sd;
            end
        end
    end
end
dist = N1_dist_sum/N1_number;
N1_amp = N1_amp_sum/N1_number;
N1_lat = (N1_lat_sum/N1_number)/(out.other.stim.fs);
baseline_sd = baseline_sum / N1_number

%% 

NumRows = size(out.elecs,2);
N2_number = 0;
N2_amp_sum = 0;
N2_lat_sum = 0;
N2_dist_sum = 0;
for ich = 1:NumRows
    for jch = 1:NumRows
        if size(out.elecs(ich).N2,1)>=2
            if out.rejection_details(1).reject.keep(ich,jch) ==1 && out.rejection_details(2).reject.keep(ich,jch) ==1
                current_N2_amp = out.elecs(ich).N2(jch,1);
                current_N2_lat = out.elecs(ich).N2(jch,2);
                current_dist_N2 = out.elecs(ich).N2(jch,3);
    
                N2_number = N2_number + 1;
                N2_amp_sum = N2_amp_sum + current_N2_amp;
                N2_lat_sum = N2_lat_sum + current_N2_lat;
                N2_dist_sum = N2_dist_sum + current_dist_N2;
            end
        end
    end
end

N2_amp = N2_amp_sum/N2_number;
N2_lat = (N2_lat_sum/N2_number)/(out.other.stim.fs);
N2_dist = N2_dist_sum/N2_number;
if N1_number == N2_number
    num_keeps = N2_number;
end
%%
% chosen_row = 5;
% if size(out.elecs(5).stim_idx,1)~=0
%     stim_start = out.elecs(5).stim_idx;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)~=0
%     stim_start = out.elecs(7).stim_idx;
%     chosen_row = 7;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)~=0
%     stim_start = out.elecs(11).stim_idx;
%     chosen_row = 11;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)~=0
%     stim_start = out.elecs(85).stim_idx;
%     chosen_row = 85;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)~=0
%     stim_start = out.elecs(18).stim_idx;
%     chosen_row = 18;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)~=0
%     stim_start = out.elecs(9).stim_idx;
%     chosen_row = 9;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)~=0
%     stim_start = out.elecs(10).stim_idx;
%     chosen_row = 10;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)~=0
%     stim_start = out.elecs(23).stim_idx;
%     chosen_row = 23;
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)==0 && size(out.elecs(3).stim_idx,1)~=0
%     stim_start = out.elecs(3).stim_idx;
%     chosen_row = 3;
% end
% if size(out.chLabels,1)>=185 
%     if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)==0 && size(out.elecs(3).stim_idx,1)==0 && size(out.elecs(185),1)~=0
%         stim_start = out.elecs(185).stim_idx;
%         chosen_row = 185;
%     end
% end
% if size(out.elecs(5).stim_idx,1)==0 && size(out.elecs(7).stim_idx,1)==0 && size(out.elecs(11).stim_idx,1)==0 && size(out.elecs(85).stim_idx,1)==0 && size(out.elecs(18).stim_idx,1)==0 && size(out.elecs(9).stim_idx,1)==0 && size(out.elecs(10).stim_idx,1)==0 && size(out.elecs(23).stim_idx,1)==0 && size(out.elecs(3).stim_idx,1)==0 && size(out.elecs(185).stim_idx,1)==0 && size(out.elecs(13),1)~=0
%     stim_start = out.elecs(13).stim_idx;
%     chosen_row = 13;
% end
% 
% idx_before_stim = 30;
% n1_time = [11e-3 50e-3];
% loose_n1_time = [11e-3 50e-3];
% n2_time = [50e-3 300e-3];
% stim_time = [-5e-3 10e-3];
% tight_stim_time = [-5e-3 10e-3];
% stim_val_thresh = 1e3;
% rel_thresh = 3;
% fs = out.other.stim.fs;
% max_crossings = 3;
% 
% n1_idx = floor(n1_time*fs);
% loose_n1_idx = floor(loose_n1_time*fs);
% n2_idx = floor(n2_time*fs);
% stim_indices = floor(stim_time*fs);
% tight_stim_indices = floor(tight_stim_time*fs);
% stim_idx = stim_start;
% 
% temp_n1_idx = n1_idx + stim_idx - 1;
% temp_loose_n1_idx = loose_n1_idx + stim_idx - 1;
% temp_n2_idx = n2_idx + stim_idx - 1;
% temp_stim_idx = stim_indices + stim_idx - 1;
% temp_tight_stim = tight_stim_indices + stim_idx-1;
% 
% row_idx = size(out.elecs(chosen_row).detrend_flipped,2);
% eeg = out.all_ch_avgs(:,row_idx+1);
% 
% baseline = mean(eeg(1:stim_idx-idx_before_stim));
% n1_eeg = eeg(temp_n1_idx(1):temp_n1_idx(2));
% n2_eeg = eeg(temp_n2_idx(1):temp_n2_idx(2));
% both_eeg = eeg(temp_n1_idx(1):temp_n2_idx(2));
% 
% n1_eeg_abs = abs(n1_eeg-baseline);
% n2_eeg_abs = abs(n2_eeg-baseline);
% both_eeg_abs = abs(both_eeg);
% 
% baseline_sd = std(eeg(1:stim_idx-idx_before_stim))
% n1_z_score = n1_eeg_abs/baseline_sd;
% n2_z_score = n2_eeg_abs/baseline_sd;
% all_z_score = both_eeg_abs/baseline_sd;
% n1_x_vals = n1_idx+(out.other.stim.fs/2):(out.other.stim.fs/2)+length(n1_z_score)-1+n1_idx;
% plot (n1_x_vals, n1_z_score, 'green')
% n2_x_vals = n2_idx+(out.other.stim.fs/2):(out.other.stim.fs/2)+length(n2_z_score)-1+n2_idx;
% plot (n2_x_vals, n2_z_score, 'green')
% 
% % OPTION B: MAX PEAK 
% % find the identity of the peaks
% choose_qualifier_MaxAmp = 0;
% n1_type_MaxAmp = 3; % DNE to start
% [pks1_MaxAmp,locs1_MaxAmp] = findpeaks(n1_z_score,'MinPeakDistance',5e-3*fs);
% [n1_peak_MaxAmp,I1_MaxAmp] = max(pks1_MaxAmp); % find the biggest
% n1_peak_idx_MaxAmp = round(locs1_MaxAmp(I1_MaxAmp));
% if I1_MaxAmp>=0.01
%     n1_peak_MaxAmp = n1_peak_MaxAmp;
%     n1_peak_idx_MaxAmp = n1_peak_idx_MaxAmp;
% else
%     n1_peak_MaxAmp = nan;
%     n1_peak_idx_MaxAmp = nan;
% end
% if isempty(n1_peak_MaxAmp)
%     n1_peak_MaxAmp = nan;
%     n1_peak_idx_MaxAmp = nan;
% end
% 
% if ~isnan(n1_peak_idx_MaxAmp)
%     n1_peak_idx_MaxAmp = n1_peak_idx_MaxAmp + temp_n1_idx(1) - 1 - stim_idx - 1;
%     n1_objective_idx_MaxAmp = n1_peak_idx_MaxAmp + stim_start;
%     n1_y_MaxAmp = eeg(n1_objective_idx_MaxAmp,1);
% 
%     if n1_y_MaxAmp ==abs(n1_y_MaxAmp)
%         n1_type_MaxAmp = 1; %maximum
%     end
%     if n1_y_MaxAmp ~= abs(n1_y_MaxAmp)
%         n1_type_MaxAmp = 0; %minimum
%     end
% end
% 
% choose_qualifier_MaxAmp = choose_qualifier_MaxAmp + n1_peak_MaxAmp;
% 
% n2_control_var_MaxAmp = 0;
% [pks2_MaxAmp,locs2_MaxAmp] = findpeaks(n2_z_score,'MinPeakDistance',5e-3*fs);
% n2_type_MaxAmp = 4; % DNE to start
% 
% while n2_control_var_MaxAmp ==0
%     [n2_peak_MaxAmp,I2a_MaxAmp] = max(pks2_MaxAmp); % find the biggest
%     n2_peak_idx_MaxAmp = round(locs2_MaxAmp(I2a_MaxAmp));
%     if isempty(n2_peak_MaxAmp)
%         n2_peak_MaxAmp = nan;
%         n2_peak_idx_MaxAmp = nan;
%         n2_control_var_MaxAmp = 1;
%     end
%     if n2_peak_MaxAmp <5
%         n2_peak_MaxAmp = nan;
%         n2_peak_idx_MaxAmp = nan;
%         n2_control_var_MaxAmp = 1;
%     end
%     if isnan(n2_peak_MaxAmp)
%         aaaah = 99
%         n2_control_var_MaxAmp = 1;
%         n2_peak_MaxAmp = nan;
%         n2_peak_idx_MaxAmp = nan;
%     end
%     if ~isnan(n2_peak_idx_MaxAmp)
%         n2_peak_idx_MaxAmp = n2_peak_idx_MaxAmp + temp_n2_idx(1) - 1 - stim_idx - 1;
%         n2_objective_idx_MaxAmp = n2_peak_idx_MaxAmp + stim_start;
%         n2_y_MaxAmp = eeg(n2_objective_idx_MaxAmp,1);
% 
%         n2_control_var_MaxAmp = 1;
% 
%     end
% end
% choose_qualifier_MaxAmp = choose_qualifier_MaxAmp + n2_peak_MaxAmp;
% 
% 
% n1_peak = n1_peak_MaxAmp
% n1_peak_time = (n1_peak_idx_MaxAmp) / out.other.stim.fs
% n2_peak = n2_peak_MaxAmp
% n2_peak_time = (n2_peak_idx_MaxAmp) / out.other.stim.fs


end