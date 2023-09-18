function out = AA_distance_vs_amp_and_lat(out)
NumRows= size(out.elecs,2);
for n=1:NumRows;
    if ~isempty(out.elecs(n).N1)
        n;
        stim_x = out.chLabels{n,2};
        stim_y = out.chLabels{n,3};
        stim_z = out.chLabels{n,4};
       
        for m=1:NumRows;
            resp_x = out.chLabels{m,2};
            resp_y = out.chLabels{m,3};
            resp_z = out.chLabels{m,4};
            val = sqrt(((stim_x - resp_x)^2)+((stim_y - resp_y)^2)+((stim_z - resp_z)^2));
        
            test_val = 2;
            test_check = size(test_val);
        
            class(val);
            checker = size(val);
        
            if (checker==test_check);
                out.elecs(n).N1(m,3) = val;
                out.elecs(n).N2(m,3) = val;
            end
        
            if (checker~=test_check);
                out.elecs(n).N1(m,3) = nan;
                out.elecs(n).N2(m,3) = nan;
            end
        end
    end
end
%% 
NumRows = size(out.elecs,2);
out.AnalysisData_N1 = NaN((NumRows*NumRows),6);
out.AnalysisData_N2 = NaN((NumRows*NumRows),6);
%%
for n=1:NumRows
    if ~isempty(out.elecs(n).N1)
        ich = n;
        for m=1:NumRows
            c = ((n-1)*(NumRows))+m;
            
            out.elecs(n);
    
            % if out.rejection_details(1).reject.keep(n,m)~=1
            %     out.elecs(n).N1(m,1) = nan;
            %     out.elecs(n).N1(m,2) = nan;
            %     out.elecs(n).N1(m,3) = nan;
            %     out.elecs(n).N1(m,4) = nan;
            %     out.AnalysisData_N1(c,1) = nan;
            %     out.AnalysisData_N1(c,2) = nan;
            %     out.AnalysisData_N1(c,3) = nan;
            %     out.AnalysisData_N1(c,4) = nan;
            %     out.AnalysisData_N1(c,5) = n;
            %     out.AnalysisData_N1(c,6) = m;
            % end
            % if out.rejection_details(2).reject.keep(n,m)~=1
            %     out.elecs(n).N1(m,1) = nan;
            %     out.elecs(n).N1(m,2) = nan;
            %     out.elecs(n).N1(m,3) = nan;
            %     out.elecs(n).N1(m,4) = nan;
            %     out.AnalysisData_N1(c,1) = nan;
            %     out.AnalysisData_N1(c,2) = nan;
            %     out.AnalysisData_N1(c,3) = nan;
            %     out.AnalysisData_N1(c,4) = nan;
            %     out.AnalysisData_N1(c,5) = n;
            %     out.AnalysisData_N1(c,6) = m;
            % end
            % if out.rejection_details(1).reject.keep(n,m)==1
                out.AnalysisData_N1(c,1) = out.elecs(n).N1(m,1);
                out.AnalysisData_N1(c,2) = (out.elecs(n).N1(m,2))/out.other.stim.fs;
                out.AnalysisData_N1(c,3) = out.elecs(n).N1(m,3);
                if size(out.elecs(ich).N1,2)==4
                    out.AnalysisData_N1(c,4) = out.elecs(n).N1(m,4);
                end
                if size(out.elecs(ich).N1,2)==3
                    out.AnalysisData_N1(c,4) = nan;
                end
                out.AnalysisData_N1(c,5) = n;
                out.AnalysisData_N1(c,6) = m;
            % end
            % if out.rejection_details(2).reject.keep(n,m)==1
                out.AnalysisData_N2(c,1) = out.elecs(n).N2(m,1);
                out.AnalysisData_N2(c,2) = (out.elecs(n).N2(m,2))/out.other.stim.fs;
                out.AnalysisData_N2(c,3) = out.elecs(n).N2(m,3);
                if size(out.elecs(ich).N2,2)==4
                    out.AnalysisData_N2(c,4) = out.elecs(n).N2(m,4);
                end
                if size(out.elecs(ich).N2,2)==3
                    out.AnalysisData_N2(c,4) = nan;
                end
                out.AnalysisData_N2(c,5) = n;
                out.AnalysisData_N2(c,6) = m;
            % end
    
            % stim_matter_type = out.chLabels{n,6};
            % stim_in_soz = out.chLabels{n,10};
            % if size (stim_matter_type,1) ~= 0
            %     for m=1:NumRows
            %         out.elecs(n).N1(m,4) = stim_matter_type;
            %         out.elecs(n).N2(m,4) = stim_matter_type;
            %         out.elecs(n).N1(m,6) = stim_in_soz;
            %         out.elecs(n).N2(m,6) = stim_in_soz;
            %     end
            % end
            % if size(stim_matter_type,1) == 0
            %     out.elecs(n).N1(m,4) = nan;
            %     out.elecs(n).N2(m,4) = nan;
            %     out.elecs(n).N1(m,6) = nan;
            %     out.elecs(n).N2(m,6) = nan;
            % end
        end
    
    
        % resp_matter_type = out.chLabels{n,6};
        % resp_in_soz = out.chLabels{n,10};
        % if size (resp_matter_type,1) ~= 0
        %     for m=1:NumRows
        %         out.elecs(m).N1(n,5) = resp_matter_type;
        %         out.elecs(m).N2(n,5) = resp_matter_type;
        %         out.elecs(m).N1(n,7) = resp_in_soz;
        %         out.elecs(m).N2(n,7) = resp_in_soz;
        %     end
        % end
        % if size(resp_matter_type,1) == 0
        %     for m=1:NumRows
        %         out.elecs(m).N1(n,5) = nan;
        %         out.elecs(m).N2(n,5) = nan;
        %         out.elecs(m).N1(n,7) = nan;
        %         out.elecs(m).N2(n,7) = nan;
        %     end
        % end
    end
end
end