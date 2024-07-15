function BP_random_rejections_keeps(out)

%% Parameters
pretty = 0;
n_to_plot = 25; % how many total to show
n_per_line = 5;
n_lines = 5;
n1_time = [15e-3 50e-3];
n2_time = [50e-3 300e-3];
zoom_times = [-300e-3 300e-3];
zoom_factor = 2;
which_n = 1;

%% Get various path locations
locations = cceps_files; % Need to make a file pointing to you own path
pwfile = locations.pwfile;
loginname = locations.loginname;
script_folder = locations.script_folder;
results_folder = locations.results_folder;

% add paths
addpath(genpath(script_folder));
if isempty(locations.ieeg_folder) == 0
    addpath(genpath(locations.ieeg_folder));
end
name = out.name;
out_folder = [results_folder,'validation/',name,'/'];
if ~exist(out_folder,'dir')
    mkdir(out_folder)
end

%% Pick intracranial chs with bipolar signal

% keep_chs = get_chs_to_ignore(out.bipolar_labels)
% save('keep_chs')
all_ones = ones(size(out.elecs,2), 1);
% keep_vals = zeros(size(out.elecs,2),1);
for ich = 1:size(out.elecs,2)
    if isempty(out.elecs(ich).arts)
        keep_vals(ich,1) = 0;
    else
        keep_vals(ich,1) = 1;
    end
end

keep_chs = (keep_vals == all_ones);
% save('keep_chs')


%% Get rejection details arrays
thresh = out.rejection_details(which_n).thresh;
which = out.rejection_details(which_n).which;

sig_avg = out.rejection_details(which_n).reject.sig_avg;
pre_thresh = out.rejection_details(which_n).reject.pre_thresh;
at_thresh = out.rejection_details(which_n).reject.at_thresh;
no_both = out.rejection_details(which_n).reject.no_both;
keep = out.rejection_details(which_n).reject.keep;
deriv = out.rejection_details(which_n).reject.deriv;
num1 = out.rejection_details(which_n).reject.num1;
num2 = out.rejection_details(which_n).reject.num2;
num3 = out.rejection_details(which_n).reject.num3;
num4 = out.rejection_details(which_n).reject.num4;
num5 = out.rejection_details(which_n).reject.num5;
not_grey = out.rejection_details(which_n).reject.not_grey;
SOZ_restriction = out.rejection_details(which_n).reject.SOZ_restrictions;
empty = out.rejection_details(which_n).reject.empty;

any_reject = sig_avg == 1| pre_thresh == 1 | at_thresh == 1 | no_both == 1 | num1 ==1 | num2 ==1 | num3 ==1 | num4 ==1| num5 ==1 | not_grey ==1| empty ==1 | SOZ_restriction ==1 | deriv ==1; 

% Calculate total numbers
nkeep = sum(keep(:) == 1);
nreject = sum(any_reject(:) == 1);
nunstim = sum(isnan(keep(:)));

if nunstim+nreject+nkeep ~= size(keep,1)*size(keep,1)
    a = nkeep+nreject+nunstim;
    b = size(keep,1)*size(keep,1);
    error('numbers do not add up');
end

% Loop through rejection types
for j = 1:2
    if j == 1
        thing = keep;
        cat = 'New Keep';
    else
        thing = any_reject;
        cat = 'Reject Any';
    end
    
    meet_criteria = find(thing==1);
    
    % Restrict to those on keep chs
    [row,col] = ind2sub(size(keep),meet_criteria);
    meet_criteria(keep_chs(row) == false) = [];
    col(keep_chs(row) == false) = [];
    meet_criteria(keep_chs(col) == false) = [];
    
    % Initialize figure
    figure
    set(gcf,'position',[100 100 1200 1000])
    t = tiledlayout(n_lines,n_per_line,'padding','compact','tilespacing','compact');
    
 
    % Pick a random N
    to_plot = randsample(meet_criteria,25);
    %to_plot=meet_criteria; %(use if <25 keeps, else comment out and use above)
    % to_plot = [15527, 5719, 42355, 33043, 54727, 55243, 55759, 15785, 15241, 5679, 15269, 3671, 16817, 7759, 15757, 7525, 16559, 62209, 54985, 15011, 8017, 8045, 62725, 61951, 33301];
    
    % Loop through these
    for i=1:length(to_plot)
        
        ind = to_plot(i);
        
        % convert this to row and column
        [row,col] = ind2sub(size(keep),ind);
        
        row;
        col;
        
        % get why it was rejected
        why = nan;
        if j == 2
            if sig_avg(row,col) == 1
                why = 'averaging';
            end
            if pre_thresh(row,col) == 1
                % if ~isnan(why)
                %     error('what');
                % end
                why = 'artifact';
            end
            if at_thresh(row,col) == 1
                % if ~isnan(why)
                %     error('what');
                % end
                why = 'threshold';
            end
            if no_both(row,col) == 1
                if isnan(why)
                    why = 'no both';
                end
            end
            % if no_n1(row,col) == 1
            %     why = 'no n1 peak';
            % end
            if no_both(row,col) == 1
                why = 'lacks N1 or N2';
            end
            if num1(row,col) == 1
                why = 'stim amp above thresh';
            end
            if num2(row,col) == 1
                why = 'signal amp above thresh';
            end
            if num3(row,col) == 1
                why = 'excess crossings in N1 period';
            end
            if num4(row,col) == 1
                why = 'lack baseline return bt. stim & N1';
            end
            if num5(row,col) == 1
                why = 'lack return to baseline after N1';
            end
            if not_grey(row,col) == 1
                why = 'not grey';
            end
            if SOZ_restriction(row,col) == 1
                why = 'SOZ restriction';
            end
            % if all_bad(row,col) == 1
            %     why = 'ch labeled bad';
            % end
            if empty(row,col) == 1
                why = 'empty';
            end
            if deriv(row,col) == 1
                if isnan(why)
                    why = 'deriv';
                end
            end
        end
        
        % Get the waveform
        row
        col
        avg = out.elecs(row).avg(:,col);
        times = out.elecs(row).times;
        eeg_times = convert_indices_to_times(1:length(avg),out.other.stim.fs,times(1));
        wav =  out.elecs(row).(which)(col,:);
        
        stim_idx = out.elecs(row).stim_idx;
        wav_idx = wav(2)+stim_idx+1;
        wav_time = convert_indices_to_times(wav_idx,out.other.stim.fs,times(1));
        n1_idx = floor(n1_time*out.other.stim.fs);
        n2_idx = floor(n2_time*out.other.stim.fs);
        temp_n1_idx = n1_idx + stim_idx - 1;
        temp_n2_idx = n2_idx + stim_idx - 1;
        
        if j==1
            if ~isnan(out.elecs(row).N1(col,1))
                % Plot
                nexttile
                plot(eeg_times,avg,'k','linewidth',2);
                hold on
                
                if (out.elecs(row).N1(col,1))~=0
                    x = (out.elecs(row).N1(col,2)/out.other.stim.fs);
                    x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1);
                    y = out.elecs(row).avg(x_indx,col);
                    plot(x,y,'bX','markersize',15,'linewidth',4);
                    if ~pretty
                    text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f',...
                        which,wav(1)), 'fontsize',9)
                    end
                end
                hold on
                if ~isnan(out.elecs(row).N2(col,2))
                    x = (out.elecs(row).N2(col,2)/out.other.stim.fs);
                    x_indx = out.elecs(row).N2(col,2)+stim_idx+1;
                    y = out.elecs(row).avg(x_indx,col);
                    plot(x,y,'rX','markersize',15,'linewidth',4);
                end
        
                %xlim([eeg_times(1) eeg_times(end)])
                xlim([zoom_times(1) zoom_times(2)]);
                
                % Zoom in (in the y-dimension) around the maximal point in the N1
                % time period
                height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
                if ~any(isnan(avg))
                    ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
                end
                
                
                labels = out.bipolar_labels;
                stim_label = labels{row};
                resp_label = labels{col};
                pause(0.1)
                xl = xlim;
                yl = ylim;
                if ~pretty
                    text(xl(1),yl(2),sprintf('Stim: %s\nResponse: %s',stim_label,resp_label),...
                        'horizontalalignment','left',...
                        'verticalalignment','top','fontsize',9);
                end
                plot([0 0],ylim,'k--');
                set(gca,'fontsize',9)
                if pretty
                    yticklabels([])
                    %xticklabels([])
                    xtl = xticklabels;
                    xtlc = cellfun(@(x) sprintf('%s s',x),xtl,'uniformoutput',false);
                    %xlabel('Time (s)')
                    xticklabels(xtlc)
                end
            end
        end
        if j == 2
            %if ~isnan(out.elecs(row).N1(col,1))
                % Plot
                nexttile
                plot(eeg_times,avg,'k','linewidth',2);
                hold on
                
                if (out.elecs(row).N1(col,1))~=0 && ~isnan(out.elecs(row).N1(col,1))
                    x = (out.elecs(row).N1(col,2)/out.other.stim.fs);
                    x_indx = round(out.elecs(row).N1(col,2)+stim_idx+1);
                    y = out.elecs(row).avg(x_indx,col);
                    plot(x,y,'bX','markersize',15,'linewidth',4);
                    if ~pretty
                    text(wav_time+0.01,avg(round(wav_idx)),sprintf('%s z-score: %1.1f',...
                        which,wav(1)), 'fontsize',9)
                    end
                end
                hold on
                if ~isnan(out.elecs(row).N2(col,2))
                    x = (out.elecs(row).N2(col,2)/out.other.stim.fs);
                    x_indx = out.elecs(row).N2(col,2)+stim_idx+1
                    y = out.elecs(row).avg(x_indx,col);
                    plot(x,y,'rX','markersize',15,'linewidth',4);
                end
        
                %xlim([eeg_times(1) eeg_times(end)])
                xlim([zoom_times(1) zoom_times(2)]);
                
                % Zoom in (in the y-dimension) around the maximal point in the N1
                % time period
                height = max(abs(avg(temp_n1_idx(1):temp_n1_idx(2))-median(avg)));
                if ~any(isnan(avg))
                    ylim([median(avg)-zoom_factor*height,median(avg)+zoom_factor*height]);
                end
                
                
                labels = out.bipolar_labels;
                stim_label = labels{row};
                resp_label = labels{col};
                pause(0.1)
                xl = xlim;
                yl = ylim;
                if ~pretty
                    text(xl(1),yl(2),sprintf('Stim: %s\nResponse: %s',stim_label,resp_label),...
                        'horizontalalignment','left',...
                        'verticalalignment','top','fontsize',9);
                end
                plot([0 0],ylim,'k--');
                set(gca,'fontsize',9)
                if pretty
                    yticklabels([])
                    %xticklabels([])
                    xtl = xticklabels;
                    xtlc = cellfun(@(x) sprintf('%s s',x),xtl,'uniformoutput',false);
                    %xlabel('Time (s)')
                    xticklabels(xtlc)
                end
                title(why)
            %end
        end
    end
    
    if pretty == 0
        title(t,sprintf('%s %s z-score threshold %1.1f',cat,which,thresh));
    end
    
    % Save the figure
    if pretty
        fname = sprintf('%s_%sthresh_%d_pretty',cat,which,thresh);
    else
        fname = sprintf('%s_%sthresh_%d',cat,which,thresh);
    end
    print(gcf,[out_folder,fname],'-dpng');
    

    
end
