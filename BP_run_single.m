out = AA_alternative_filtering(out);
out = AA_Running_RejectOrKeep_RW(out);
out = AA_new_build_network(out,0);
out = AA_distance_vs_amp_and_lat(out);
out = AA_require_both_Ns(out);
out = AA_modify_analysisdata(out);

 AA_random_rejections_keeps(out);
%% Get various path locations and save if desired
% locations = cceps_files; % Need to make a file pointing to you own path
% results_folder = locations.results_folder;
% outdir = [results_folder,'out_files/'];
% pt=sprintf('pt%c',(string(inputdlg('if there is a file with pt_1 or pt_2 enter the pt number here, else enter 0'))));
% dataName=out.name;
% if pt(3)~='0'
%     dataName=sprintf('%s_%s',dataName,pt);
% else
% end
% %%
% save([outdir,sprintf('out_bipolar_final_%s',dataName)],'out','-v7.3');