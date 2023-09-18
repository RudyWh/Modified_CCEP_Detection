function out = AA_modify_analysisdata(out)
NumRows = size(out.elecs,2);
for n=1:NumRows
    for m=1:NumRows
        c = ((n-1)*(NumRows))+m;
        c;
        if out.AnalysisData_N1(c,2)>= 0.017
            out.AnalysisData_N1(c,4) = 2;
            out.AnalysisData_N2(c,4) = 2;
        end
        if out.AnalysisData_N1(c,2)<0.017
            out.AnalysisData_N1(c,4) = 1;
            out.AnalysisData_N2(c,4) = 1;
        end
        if isnan(out.AnalysisData_N1(c,2))
            out.AnalysisData_N2(c,1) = nan;
            out.AnalysisData_N2(c,2) = nan;
            out.AnalysisData_N2(c,3) = nan;
            out.AnalysisData_N2(c,4) = nan;
        end
    end
end
end