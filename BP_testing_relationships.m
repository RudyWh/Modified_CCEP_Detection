tiledlayout(2,2)

% adjusted N1 amplitude
nexttile
scatter(cell2mat(BP_adjusted_results(:,8)),cell2mat(BP_adjusted_results(:,2)))
hold on 
title("Adjusted N1 amplitude")
X = [ones(length(cell2mat(BP_adjusted_results(:,8))),1) cell2mat(BP_adjusted_results(:,8))];
b = X\(cell2mat(BP_adjusted_results(:,2)));
yCalc2 = X*b;
plot(cell2mat(BP_adjusted_results(:,8)),yCalc2)
Rsq = 1 - sum(((cell2mat(BP_adjusted_results(:,2))) - yCalc2).^2)/sum(((cell2mat(BP_adjusted_results(:,2))) - mean(cell2mat(BP_adjusted_results(:,2)))).^2);
[rho,pval] = corr((cell2mat(BP_adjusted_results(:,8))),(cell2mat(BP_adjusted_results(:,2))));
h(1) = plot(NaN, NaN, '-r');
h(2) = plot(NaN, NaN, 'ow');
h(3) = plot(NaN, NaN, 'ow');
pval_string = "pval: "+string(pval);
rsq_string = "r^2: "+ string(Rsq);
legend(h,"trendline",pval_string, rsq_string)


% adjusted N2 amplitude
nexttile
scatter(cell2mat(BP_adjusted_results(:,8)),cell2mat(BP_adjusted_results(:,3)))
hold on
title("Adjusted N2 amplitude")
X = [ones(length(cell2mat(BP_adjusted_results(:,8))),1) cell2mat(BP_adjusted_results(:,8))];
b = X\(cell2mat(BP_adjusted_results(:,3)));
yCalc2 = X*b;
plot(cell2mat(BP_adjusted_results(:,8)),yCalc2)
Rsq = 1 - sum(((cell2mat(BP_adjusted_results(:,3))) - yCalc2).^2)/sum(((cell2mat(BP_adjusted_results(:,3))) - mean(cell2mat(BP_adjusted_results(:,3)))).^2);
[rho,pval] = corr((cell2mat(BP_adjusted_results(:,8))),(cell2mat(BP_adjusted_results(:,3))));
h(1) = plot(NaN, NaN, '-r');
h(2) = plot(NaN, NaN, 'ow');
h(3) = plot(NaN, NaN, 'ow');
pval_string = "pval: "+string(pval);
rsq_string = "r^2: "+ string(Rsq);
legend(h,"trendline",pval_string, rsq_string)


% adjusted N1 latency
nexttile
scatter(cell2mat(BP_adjusted_results(:,8)),cell2mat(BP_adjusted_results(:,4)))
hold on
title("Adjusted N1 latency")
X = [ones(length(cell2mat(BP_adjusted_results(:,8))),1) cell2mat(BP_adjusted_results(:,8))];
b = X\(cell2mat(BP_adjusted_results(:,4)));
yCalc2 = X*b;
plot(cell2mat(BP_adjusted_results(:,8)),yCalc2)
Rsq = 1 - sum(((cell2mat(BP_adjusted_results(:,4))) - yCalc2).^2)/sum(((cell2mat(BP_adjusted_results(:,4))) - mean(cell2mat(BP_adjusted_results(:,4)))).^2);
[rho,pval] = corr((cell2mat(BP_adjusted_results(:,8))),(cell2mat(BP_adjusted_results(:,4))));
h(1) = plot(NaN, NaN, '-r');
h(2) = plot(NaN, NaN, 'ow');
h(3) = plot(NaN, NaN, 'ow');
pval_string = "pval: "+string(pval);
rsq_string = "r^2: "+ string(Rsq);
legend(h,"trendline",pval_string, rsq_string)


% adjusted N2 latency
nexttile
scatter(cell2mat(BP_adjusted_results(:,8)),cell2mat(BP_adjusted_results(:,5)))
hold on
title("Adjusted N2 latency")
X = [ones(length(cell2mat(BP_adjusted_results(:,8))),1) cell2mat(BP_adjusted_results(:,8))];
b = X\(cell2mat(BP_adjusted_results(:,5)));
yCalc2 = X*b;
plot(cell2mat(BP_adjusted_results(:,8)),yCalc2)
Rsq = 1 - sum(((cell2mat(BP_adjusted_results(:,5))) - yCalc2).^2)/sum(((cell2mat(BP_adjusted_results(:,5))) - mean(cell2mat(BP_adjusted_results(:,5)))).^2);
[rho,pval] = corr((cell2mat(BP_adjusted_results(:,8))),(cell2mat(BP_adjusted_results(:,5))));
h(1) = plot(NaN, NaN, '-r');
h(2) = plot(NaN, NaN, 'ow');
h(3) = plot(NaN, NaN, 'ow');
pval_string = "pval: "+string(pval);
rsq_string = "r^2: "+ string(Rsq);
legend(h,"trendline",pval_string, rsq_string)

% %%
% figure
% tiledlayout(2,2)
% 
% % unadjusted N1 amplitude
% nexttile
% scatter(cell2mat(Excel_dat_unadjusted(:,8)),cell2mat(Excel_dat_unadjusted(:,2)))
% hold on 
% title("Unadjusted N1 amplitude")
% X = [ones(length(cell2mat(Excel_dat_unadjusted(:,8))),1) cell2mat(Excel_dat_unadjusted(:,8))];
% b = X\(cell2mat(Excel_dat_unadjusted(:,2)));
% yCalc2 = X*b;
% plot(cell2mat(Excel_dat_unadjusted(:,8)),yCalc2)
% Rsq = 1 - sum(((cell2mat(Excel_dat_unadjusted(:,2))) - yCalc2).^2)/sum(((cell2mat(Excel_dat_unadjusted(:,2))) - mean(cell2mat(Excel_dat_unadjusted(:,2)))).^2);
% [rho,pval] = corr((cell2mat(Excel_dat_unadjusted(:,8))),(cell2mat(Excel_dat_unadjusted(:,2))));
% h(1) = plot(NaN, NaN, '-r');
% h(2) = plot(NaN, NaN, 'ow');
% pval_string = "pval: "+string(pval);
% legend(h,"trendline",pval_string)
% 
% 
% % unadjusted N2 amplitude
% nexttile
% scatter(cell2mat(Excel_dat_unadjusted(:,8)),cell2mat(Excel_dat_unadjusted(:,3)))
% hold on
% title("Unadjusted N2 amplitude")
% X = [ones(length(cell2mat(Excel_dat_unadjusted(:,8))),1) cell2mat(Excel_dat_unadjusted(:,8))];
% b = X\(cell2mat(Excel_dat_unadjusted(:,3)));
% yCalc2 = X*b;
% plot(cell2mat(Excel_dat_unadjusted(:,8)),yCalc2)
% Rsq = 1 - sum(((cell2mat(Excel_dat_unadjusted(:,3))) - yCalc2).^2)/sum(((cell2mat(Excel_dat_unadjusted(:,3))) - mean(cell2mat(Excel_dat_unadjusted(:,3)))).^2);
% [rho,pval] = corr((cell2mat(Excel_dat_unadjusted(:,8))),(cell2mat(Excel_dat_unadjusted(:,3))));
% h(1) = plot(NaN, NaN, '-r');
% h(2) = plot(NaN, NaN, 'ow');
% pval_string = "pval: "+string(pval);
% legend(h,"trendline",pval_string)
% 
% 
% % unadjusted N1 latency
% nexttile
% scatter(cell2mat(Excel_dat_unadjusted(:,8)),cell2mat(Excel_dat_unadjusted(:,4)))
% hold on
% title("Unadjusted N1 latency")
% X = [ones(length(cell2mat(Excel_dat_unadjusted(:,8))),1) cell2mat(Excel_dat_unadjusted(:,8))];
% b = X\(cell2mat(Excel_dat_unadjusted(:,4)));
% yCalc2 = X*b;
% plot(cell2mat(Excel_dat_unadjusted(:,8)),yCalc2)
% Rsq = 1 - sum(((cell2mat(Excel_dat_unadjusted(:,4))) - yCalc2).^2)/sum(((cell2mat(Excel_dat_unadjusted(:,4))) - mean(cell2mat(Excel_dat_unadjusted(:,4)))).^2);
% [rho,pval] = corr((cell2mat(Excel_dat_unadjusted(:,8))),(cell2mat(Excel_dat_unadjusted(:,4))));
% h(1) = plot(NaN, NaN, '-r');
% h(2) = plot(NaN, NaN, 'ow');
% pval_string = "pval: "+string(pval);
% legend(h,"trendline",pval_string)
% 
% 
% % unadjusted N2 latency
% nexttile
% scatter(cell2mat(Excel_dat_unadjusted(:,8)),cell2mat(Excel_dat_unadjusted(:,5)))
% hold on
% title("Unadjusted N2 latency")
% X = [ones(length(cell2mat(Excel_dat_unadjusted(:,8))),1) cell2mat(Excel_dat_unadjusted(:,8))];
% b = X\(cell2mat(Excel_dat_unadjusted(:,5)));
% yCalc2 = X*b;
% plot(cell2mat(Excel_dat_unadjusted(:,8)),yCalc2)
% Rsq = 1 - sum(((cell2mat(Excel_dat_unadjusted(:,5))) - yCalc2).^2)/sum(((cell2mat(Excel_dat_unadjusted(:,5))) - mean(cell2mat(Excel_dat_unadjusted(:,5)))).^2);
% [rho,pval] = corr((cell2mat(Excel_dat_unadjusted(:,8))),(cell2mat(Excel_dat_unadjusted(:,5))));
% h(1) = plot(NaN, NaN, '-r');
% h(2) = plot(NaN, NaN, 'ow');
% pval_string = "pval: "+string(pval);
% legend(h,"trendline",pval_string)

%%
% RW todo:
% - try and get ref working through this*
% - take the unfiltered ref (if made) and rerun this as well**
% - or try filtering again (probably not)
%       TRY REJECT 0 & 3 (could also try with bipolar) 
%       Or amplitude stuff from stim sheet 
% - Caren fixing too big to open, is for CRP with ref (checking that out) 