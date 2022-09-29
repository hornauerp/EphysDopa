function table_wfs_features = wave_form_features(thiswave)

%% Normalize amplitude
thiswave = thiswave/abs(min(thiswave));

%% Function to find zero-crossings
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
tx = zci(thiswave);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spike half-width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p,l,hw] = findpeaks(-thiswave,'Minpeakheight',max(- thiswave)-0.001,'Annotate','extents','WidthReference','halfheight');
if length(p)>1
    [da,wo] = max(p);
    p = p(wo);
    hw = hw(wo);
end
table_wfs_features.Amplitude = p;
table_wfs_features.HalfWidth = hw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wave form features (Ampl, Hw, asymmetry, trough-to-peak-ratio, trough-to-peak delay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_wfs_features.assy = 0;
table_wfs_features.duration = 0;
table_wfs_features.t2pdelay = 0;
table_wfs_features.t2pratio = 0;
table_wfs_features.trough = 0;
table_wfs_features.peak1 = [];
table_wfs_features.peak2 = [];
table_wfs_features.AUCpeak1 = 0;
table_wfs_features.WidthPeak1 = 0;
table_wfs_features.AUCpeak2 = 0;
table_wfs_features.WidthPeak2 = 0;
table_wfs_features.AUCtrough = 0;
table_wfs_features.WidthTrough = 0;
peak1_wo = [];
peak2_wo = [];

% find peaks 1 and 2
[peaks_da,peaks_wo] = findpeaks(thiswave);
[trough_da,trough_wo] = min(thiswave);
limit_before_trough = trough_wo-19; % range to find them
limit_after_trough = trough_wo+20; 

[~, idx] = min(abs(tx-trough_wo(1)));
if tx(idx) > trough_wo(1)
    trough_top = tx(idx);
    if idx == 1
        trough_bottom = 1;
    else
        trough_bottom = tx(idx-1);
    end
else
    trough_bottom = tx(idx);
        if idx==length(tx)
            trough_top = length(thiswave);
        else
            trough_top = tx(idx+1);
        end
end
table_wfs_features.WidthTrough = thiswave(trough_bottom:trough_top);
table_wfs_features.AUCtrough = abs(trapz(table_wfs_features.WidthTrough));

peak1_da = max(peaks_da(find(peaks_wo < trough_wo & peaks_wo > limit_before_trough)));

if ~isempty(peak1_da)
    peak1_wo = peaks_wo(find(peaks_da == peak1_da));
    
    [~, idx] = min(abs(tx-peak1_wo(1)));
    if tx(idx) > peak1_wo(1)
        peak1_top = tx(idx);
        if idx == 1
            peak1_bottom = 1;
        else
            peak1_bottom = tx(idx-1);
        end

    else
        peak1_bottom = tx(idx);
        if idx==length(tx)
            peak1_top = length(thiswave);
        else
            peak1_top = tx(idx+1);
        end
    end
    WidthPeak1 = thiswave(peak1_bottom:peak1_top);
    table_wfs_features.AUCpeak1 = trapz(WidthPeak1);
    table_wfs_features.WidthPeak1 = WidthPeak1;
end
peak2_da =  max(peaks_da(find(peaks_wo > trough_wo & peaks_wo < limit_after_trough))); %max(peaks_da(find(peaks_wo > trough_wo & peaks_wo < limit_after_trough))); % peak post trough

if ~isempty(peak2_da)
    peak2_wo = peaks_wo(find(peaks_da == peak2_da));
    [~, idx] = min(abs(tx-peak2_wo(1)));
    if tx(idx) > peak2_wo(1)
        peak2_top = tx(idx);
        if idx == 1
            peak2_bottom = 1;
        else
            peak2_bottom = tx(idx-1);
        end
    else
        peak2_bottom = tx(idx);
        if idx==length(tx)
            peak2_top = length(thiswave);
        else
            peak2_top = tx(idx+1);
        end
    end
    WidthPeak2 = thiswave(peak2_bottom:peak2_top);
    table_wfs_features.AUCpeak2 = trapz(WidthPeak2);
    table_wfs_features.WidthPeak2 = WidthPeak2;
end



% only calculate asymmetry, if you find both peaks reliably
if ~isempty(peak1_wo) && ~isempty(peak2_wo)
    table_wfs_features.trough = [trough_da(1) trough_wo(1)];
    table_wfs_features.peak1 =  [peak1_da(1) peak1_wo(1)];
    table_wfs_features.peak2 =  [peak2_da(1) peak2_wo(1)];
    table_wfs_features.assy = (abs(peak2_da(1))-abs(peak1_da(1)))/(abs(peak2_da(1))+abs(peak1_da(1)));
end

if ~isempty(peak2_wo) 
    table_wfs_features.t2pdelay = length(trough_wo:peak2_wo);
    table_wfs_features.t2pratio = abs(trough_da)/abs(peak2_da);
end

end