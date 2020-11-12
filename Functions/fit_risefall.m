function table = fit_risefall(wfs_table,table_wfs_features)

waveform = wfs_table.single_template_spatial;
%waveform = resample(waveform,2,1); % upsample
[da,~] = max(-waveform); % to find the waveform with the maximum amplitude
[~,wo_max] = max(da);
ts = waveform(:,wo_max)';
ts = ts/abs(min(ts)); %Normalize amplitude
%[trough_val,trough_loc] = max(-ts);
%ts_cut = ts(trough_loc:end);
%[peak_val,peak_loc] = max(ts_cut);

%% Rising slope
%Rp = []; Rstat = []; R2_rising = []; Cp = []; Cstat = [];  R2_decay =[];
% lm_rising = [];lm_decay = []; 

try
    r = slewrate(ts(table_wfs_features.trough(2):end));
    table.lm_rising = r(1);
%     %[Rp,Rstat] = robustfit(1:15,ts_cut(1:15)');
%     %[R2_rising] = getrsq_robustfit(Rp,Rstat,1:15);
%     lm_rising = fitlm(1:7,ts(table_wfs_features.trough(2):table_wfs_features.trough(2)+6));
%     if lm_rising.Rsquared.Adjusted > 0.8
%         table.lm_rising = lm_rising.Coefficients.Estimate(2); %Changed
%         else
%             table.lm_rising = [];
%     end
catch
    table.lm_rising= [];
end

%% Falling slope
if ~isempty(table_wfs_features.peak2)
    peak_loc = table_wfs_features.peak2(2);
end

try
    %[Cp,Cstat] = robustfit(1:15,ts(peak_loc:peak_loc+14)');
    %[R2_decay] = getrsq_robustfit(Cp,Cstat,1:15);
    %     lm_decay = fitlm(1:10,ts(table_wfs_features.peak2(2):table_wfs_features.peak2(2)+9));
    %     if lm_decay.Rsquared.Adjusted > 0.8
    %         table.lm_decay = lm_decay.Coefficients.Estimate(2); %Changed
    %         else
    %             table.lm_decay = [];
    %     end
    d = slewrate(ts(table_wfs_features.peak2(2):end));
    table.lm_decay = d(1);
catch
    table.lm_decay= [];
end

% Gather results
%table.trough_loc = trough_loc;
%table.peak_loc = peak_loc;
table.ts = ts;


% table.Rp = Rp;
% table.R2_rising = R2_rising;
% table.Rstat = Rstat;
% table.Cp = Cp;
% table.Cstat =  Cstat;
% table.R2_decay = R2_decay;
