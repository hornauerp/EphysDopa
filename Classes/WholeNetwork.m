classdef WholeNetwork < handle

    properties (SetAccess=immutable)
        lookup_path = '/cluster/project/bsse_sdsc/BELUB/Scripts/Cell_lines_lookup.mat';
    end
    
    properties (SetObservable)
        InputPath
        Params
        Duration
        Templates
        XYElectrodes
        ActiveChannels
        Active  %(1,:) = spike times (2,:) channels of Active templates (3,:) within burst (1)
        SpikeTimes
        SamplingRate
        TemplateMatrix
        BaselineFrequency
        BaselineDistribution
        
        %Burstdata
        Bursts
        IBI = 0;
        IBIVar = 0;
        IBIMean = 0;
        BD = 0;
        BDMean = 0;
        BDVar = 0;
        IntraBF = 0;
        InterBF = 0;
        RiseTime = 0;
        RiseVelocity = 0;
        DecayTime = 0;
        DecayVelocity = 0;
        Synchronicity
        RegularityFrequency
        RegularityMagnitude
        RegularityFit
        
        %Metadata
        Species
        StemCell
        CellType
        Mutation
        Neurons
        Astrocytes
        PlatingDate
        Institute
        RecordingDate
        DIV
        ChipID
        CellLine
        Treatment
        TPT %Time past treatment
    end
    
    properties (Dependent)
        
    end
    
    %% Methods to initialize
    methods (Hidden)
        function nw = WholeNetwork(input_path,params,subsampling,method)
            if nargin > 0
                nw.InputPath = input_path;
                nw.Params = params;
                nw.getMetaData();
                nw.XYElectrodes = nw.getCoordinates();
                [nw.SpikeTimes, nw.SamplingRate] = nw.getSpikeTimes();
                nw.TemplateMatrix = nw.getTemplateMatrix();
                nw.Duration = double(max(nw.SpikeTimes(1,:))); %seconds
                nw.Templates = nw.generateTemplates();
                if isempty([nw.ActiveChannels])
                    return
                end
                if nargin > 2
                    nw.subsampleNW(subsampling,method);
                end
                if ~isempty([nw.ActiveChannels])
                    nw.getActive;
                    nw.getBaselineFrequency;
                    nw.getRegularity;
                    nw.getBurstTimes;
                    nw.getBurstStatistics;
                    nw.getSynchronicity;
                    nw.getSingleCellFeatures;
                end
            end
        end
        
        function subsampleNW(obj,subsampling,method)
            if subsampling < length([obj.ActiveChannels])
                min_nsize = 16;
                if method == "Act"
                    tmp_act = [obj.Templates(obj.ActiveChannels).Activity];
                    [~,ind] = sort(tmp_act,'descend');
                    if subsampling > 1
                        ind = ind(1:subsampling);
                    else
                        rel_subsampling = round(subsampling*length([obj.ActiveChannels]));
                        if rel_subsampling < min_nsize && min_nsize < length([obj.ActiveChannels])
                            rel_subsampling = min_nsize;
                        elseif rel_subsampling < min_nsize && min_nsize < length([obj.ActiveChannels])
                            rel_subsampling = length([obj.ActiveChannels]);
                        end
                        ind = ind(1:rel_subsampling);
                    end
                elseif method == "Rand"
                    if subsampling > 1
                        draw_n = subsampling;
                    else
                        draw_n = round(subsampling*length([obj.ActiveChannels]));
                        if draw_n < min_nsize && min_nsize < length([obj.ActiveChannels])
                            draw_n = min_nsize;
                        elseif draw_n < min_nsize && min_nsize > length([obj.ActiveChannels])
                            draw_n = length([obj.ActiveChannels]);
                        end
                    end
                    ind = datasample(1:length([obj.ActiveChannels]),draw_n,'Replace',false);
                end
                obj.ActiveChannels = obj.ActiveChannels(ind);
            end
        end
        
        function tmp = generateTemplates(obj)
            a = [];
            for i = 1:size(obj.TemplateMatrix,1)
                if obj.passCriteria(i)
                    w = permute(obj.TemplateMatrix(i,:,:),[3 2 1]);
                    w = sparse(double(w));
                    tmp(i) = Template(obj,i,obj.SpikeTimes,w);
                    a = [a i];
                end
            end
            if ~exist('tmp','var')
                tmp = [];
            end
            obj.ActiveChannels = a;
        end
        
        function check = passCriteria(obj,id)
            spk = obj.SpikeTimes(1,obj.SpikeTimes(2,:)==id);
            countsi = histcounts(diff(spk)*1000,'BinEdges',0:100);
            rv= sum(countsi(1:2))/length(spk);
            rate = numel(spk)/obj.Duration;
            amp = min(min(obj.TemplateMatrix(id,:,:)));
            wf = permute(obj.TemplateMatrix(id,:,:),[3 2 1]);
            ind = find(wf);
            [r,~] = ind2sub(size(wf),ind);
            active = unique(r)';
            check = numel(spk)> 2 && rv < 0.02 &&  rate > 0.1 && rate < 10 && length(active)>1; %amp > -15 && amp < -5 && 
        end
        
        function xy = getCoordinates(obj)
            xy = readNPY(fullfile(obj.InputPath, 'channel_positions.npy'));
        end
        
        function [spk,sr] = getSpikeTimes(obj)
            params_path = fullfile(obj.InputPath, 'params.py');
            fileID = fopen(params_path,'r');
            params_txt = fscanf(fileID,'%c');
            params_parts = strsplit(params_txt,'\n');
            sr_idx = startsWith(params_parts,'sample');
            sr = strsplit(params_parts{sr_idx},' = ');
            sr = str2double(sr{2});
            spk(1,:) = double(readNPY(fullfile(obj.InputPath, 'spike_times.npy')))/sr;
            spk(2,:) = uint16(readNPY(fullfile(obj.InputPath, 'spike_templates.npy')));
        end
        
        function tm = getTemplateMatrix(obj)
            tmpfile = fullfile(obj.InputPath, 'templates.npy');
            tm = readNPY(tmpfile);
        end
        
        function md = getMetaData(obj)
            try
                load(obj.lookup_path)
                ids = strsplit(obj.InputPath,'/');
                while ~strcmpi(ids{end},'Network')
                    ids(end) =[];
                end
                cell_line = strsplit(ids{end-3},'_');
                cell_line = str2double(cell_line{end});
                md = table2struct(Cell_lines(cell_line,1:end-1));
                for fn = fieldnames(md)'
                    obj.(fn{1}) = md.(fn{1});
                end
                obj.ChipID = ids{end-1};
                obj.RecordingDate = ids{end-2};
                pd = datetime(string(obj.PlatingDate),'InputFormat','yyMMdd');
                rd = datetime(obj.RecordingDate,'InputFormat','yyMMdd');
                obj.DIV = days(rd-pd);
                obj.CellLine = cell_line;
            catch
                warning('No metadata reference file found')
            end
        end
        
        function getBaselineFrequency(obj,p)
            if nargin == 1
                p = 0;
            end
            tmp = obj.Templates(obj.ActiveChannels);
            isi = arrayfun(@(x) diff(x.SpikeTimes),tmp,'UniformOutput',false);
            if ~isempty(isi)
                hc = histcounts([isi{:}],round(max([isi{:}])));
                sm = smoothdata(hc);
                [~,loc,w,~] = findpeaks(sm,'SortStr','descend','NPeaks',1);
                if ~isempty(loc)
                    obj.BaselineFrequency = 1/loc;
                    obj.BaselineDistribution = w;
                else
                    obj.BaselineFrequency = 0;
                    obj.BaselineDistribution = 0;
                end
            else
                obj.BaselineFrequency = 0;
                obj.BaselineDistribution = 0;
            end
            
            if p
                %%
                fontsz = 8;
                zero_idx = find(sm==0,1,'first');
                sm = sm(1:zero_idx);
                [~,min_idx] = min(sm(1:loc));
%                 figure('Color','w');
                burst_part = semilogy(sm(1:min_idx),'k');
                hold on
                burst_part.Color(4) = 0.3;
                peak_part = semilogy(min_idx:length(sm),sm(min_idx:end),'k');
                xl = xline(min_idx,'k--','Alpha',0.3);
                scatter([loc 33],sm([loc 33]),'ko')
                text(loc,sm(loc)+2500,['Baseline' newline 'frequency'],'FontSize',fontsz)
                text(32,sm(33)-26,'Harmonic','FontSize',fontsz,'HorizontalAlignment','right')
                xlabel('Interspike interval')
                ylabel('# Interspike intervals')
                axis tight
                box off
                set(gca,'FontSize',fontsz)
            end

        end
        
    end
    
    
    % Methods for further analysis
    methods
        %% Get methods for dependent properties
        function getActive(obj)
            ids = obj.ActiveChannels;
            ind = ismember(obj.SpikeTimes(2,:),ids);
            obj.Active = obj.SpikeTimes(:,ind);
        end

        
        %% Actively callable methods
        function Burst_merged = getBurstTimes(obj)
            params = obj.Params;
            if params.Bursts.N < 1 %Absolute value if N>1, otherwise relative to total activity
                N = ceil(length([obj.Active])*params.Bursts.N);
            else
                N = params.Bursts.N;
            end
            
            ISI_N = params.Bursts.ISI_N;
            
            if N<100
                N = 65;
                ISI_N = 1;
            end
            
            t = params.Bursts.merge_t*ISI_N;
            Spike.T = obj.Active(1,:);
            Spike.C = obj.Active(2,:);
            [Burst, ~] = burst_detect_isin(Spike,N,ISI_N);
            if length(Burst.T_start)>1
                t_diff = zeros(1,length(Burst.T_end)-1);
                for i = 1:(numel(Burst.T_end)-1)
                    t_diff(i) = Burst.T_start(i+1) - Burst.T_end(i);
                end
                Burst_merged.T_start = Burst.T_start(1);
                Burst_merged.T_end = [];
                Burst_merged.S = [];
                n_spk = 0;
                for i = 1:length(t_diff)
                    n_spk = n_spk + Burst.S(i);
                    if t_diff(i) > t %Paramter for burst merge
                        Burst_merged.T_start = [Burst_merged.T_start Burst.T_start(i+1)];
                        Burst_merged.T_end = [Burst_merged.T_end Burst.T_end(i)];
                        Burst_merged.S = [Burst_merged.S n_spk];
                        n_spk = 0;
                    end
                end
                t_diff_merged = zeros(1,length(Burst_merged.T_end)-1);
                for i = 1:(numel(Burst_merged.T_end)-1)
                    t_diff_merged(i) = Burst_merged.T_start(i+1) - Burst_merged.T_end(i);
                end
                
                Burst_merged.T_end = [Burst_merged.T_end Burst.T_end(end)];
                last_burst = Burst.S(end) + n_spk;
                Burst_merged.S = [Burst_merged.S last_burst];
                obj.Bursts = Burst_merged;
            else
                Burst = rmfield(Burst,'C');
                obj.Bursts = Burst;
            end
        end
        
        function obj = getBurstStatistics(obj)
            BurstInd = zeros(1,size(obj.Active,2));
            if length(obj.Bursts.S) > 1

                for i = 1:(numel(obj.Bursts.T_end)-1)
                    obj.IBI(i) = obj.Bursts.T_start(i+1) - obj.Bursts.T_end(i);
                end
                
                obj.IBI = obj.IBI';
                obj.BD = (obj.Bursts.T_end - obj.Bursts.T_start)';
                
                % Remove outliers
                [~,rm_ibi] = rmoutliers(obj.IBI);
                [~,rm_bd] = rmoutliers(obj.BD);
                rm_idx = find(reshape(rm_ibi,1,[])&reshape(rm_bd(1:end-1),1,[]));
                obj.BD(rm_idx) = [];
                ibi = obj.IBI;
                merge_idx = rm_idx(rm_idx<length(ibi)); %Remove last idx
                merge = ibi(merge_idx)+ibi(merge_idx+1);
                ibi(merge_idx+1) = merge;
                ibi(rm_idx) = [];
                obj.IBI = ibi;
                obj.Bursts.T_start(rm_idx) = [];
                obj.Bursts.T_end(rm_idx) = [];

                
                obj.IBIMean = mean(obj.IBI);
                obj.IBIVar = var(obj.IBI);
                obj.BDMean = mean(obj.BD);
                obj.BDVar = var(obj.BD);
                
                for i = 1:numel(obj.Bursts.T_end)
                    idx = find(obj.Active(1,:) >= obj.Bursts.T_start(i)...
                        & obj.Active(1,:) <= obj.Bursts.T_end(i));
                    BurstInd(idx) = 1;
                end
                
                obj.IntraBF = sum(BurstInd)/(sum(obj.BD)*length([obj.ActiveChannels]));
                obj.InterBF = (numel(BurstInd)-sum(BurstInd))/((obj.Duration-sum(obj.BD))*length([obj.ActiveChannels]));
                binning = obj.Params.Bursts.binning;
                Fs = 1/binning;
                for i = 1:length(obj.Bursts.T_end)
                    t = obj.Active(1,:);
                    t_start = obj.Bursts.T_start(i);
                    t_end = obj.Bursts.T_end(i);
                    counts = histcounts(t,t_start:binning:t_end);
                    counts_s = smoothdata(counts,'gaussian');
                    [~,pk] = max(counts_s);
                    try
                        rt(i) = risetime(counts_s(1:pk),Fs,'StateLevels',[min(counts_s(1:pk)) max(counts_s(1:pk))]);
                    catch
                        rt(i) = NaN;
                    end
                    try
                        dt(i) = falltime(counts_s(pk:end),Fs,'StateLevels',[min(counts_s(pk:end)) max(counts_s(pk:end))]);
                    catch
                        dt(i) = NaN;
                    end
                    try
                        rv(i) = slewrate(counts_s(1:pk),Fs,'StateLevels',[min(counts_s(1:pk)) max(counts_s(1:pk))]);
                    catch
                        rv(i) = NaN;
                    end
                    try
                        dv(i) = slewrate(counts_s(pk:end),Fs,'StateLevels',[min(counts_s(pk:end)) max(counts_s(pk:end))]);
                    catch
                        dv(i) = NaN;
                    end
                end
                obj.RiseTime = mean(rt,'omitnan');
                obj.DecayTime = mean(dt,'omitnan');
                obj.RiseVelocity = mean(rv,'omitnan');
                obj.DecayVelocity = mean(dv,'omitnan');
            end
            obj.Active(3,:) = BurstInd;
            for i = obj.ActiveChannels
                t = obj.Templates(i);
                t.getInburstRatio(BurstInd(obj.Active(2,:)==i));
            end
        end
        
        function obj = getRegularity(obj)
            % Calculate Regularity
            binning = obj.Params.Regularity.binning;
            Fs = 1/binning;
            spk_nw = histc(obj.Active(1,:),0:binning:obj.Duration);
            spk_nw = spk_nw./max(spk_nw);
            NFFT = length(spk_nw);
            F = (0 : 1/NFFT : 1/2-1/NFFT)*Fs;
            TEMP = fft(spk_nw,NFFT);
            TEMP(1) = 0;
            [mag,idx] = findpeaks(abs(TEMP(1:NFFT/2)),'NPeaks',1,'SortStr','descend');
            obj.RegularityFrequency = F(idx);
            obj.RegularityMagnitude = mag;
            
            binned = abs(TEMP(1:NFFT/2));
            binned = binned./max(binned);
            pks = NaN;
            l = [];
            p = [];
            w = [];
            while ~isempty(pks)
                try
                    [pks,locs,wd,~] = findpeaks(binned,'NPeaks',1,'SortStr','descend');
                catch
                    break
                end
                binned = binned(locs:end);
                l = [l locs];
                p = [p pks];
                w = [w wd];
            end
            x = cumsum(l)';
            cutoff = find(x>max(x)/3,1,'first');
           if cutoff<4
               cutoff = 4;
           end
           x = x(1:cutoff);
           p = p(1:cutoff);
           log_p = log10(p)-min(log10(p));
            try
                f = fit(x,log_p','exp1');
                obj.RegularityFit = f.b;
            catch
                obj.RegularityFit = NaN;
            end
        end

        
        function nc = getSynchronicity(obj)
            if isempty(which('CCGHeart'))
                cd('Functions')
                mex('CCGHeart.c');%Might need to specify the correct path here
                cd ..
            end
            if length(obj.ActiveChannels)>1
                spk = obj.Active(1,:);
                ch = obj.Active(2,:);
                [ccg,~] = CCG(spk,ch,'duration',1,'binSize',0.01);
                ccg_max = max(ccg,[],1);
                cc = permute(ccg_max,[3 2 1]);
                arrayfun(@(x) x.normCon(cc,obj),obj.Templates([obj.ActiveChannels]),'UniformOutput',false);
                nc = vertcat(obj.Templates.ConVector);
                obj.Synchronicity = mean(nc,'all');
            else
                obj.Synchronicity = NaN;
            end
        end
        
        function getSingleCellFeatures(obj)
            arrayfun(@(x) x.inferFeatures,obj.Templates(obj.ActiveChannels));
        end
        
        function saveNW(obj)
           nw = obj;
           nw.SpikeTimes = [];
           nw.Amplitudes = [];
           nw.TemplateMatrix = [];
           nw.Templates =  [];
           nw.temp_x = [];
           nw.temp_y = [];
           nw.Active = [];
           nw.ActiveChannels = [];
           save('nw.mat','nw');
        end
        
        function saveTMP(obj,save_path)
            nw = obj;
            nw.SpikeTimes = [];
            nw.TemplateMatrix = [];
            nw.temp_x = [];
            nw.temp_y = [];
            nw.Active = [];
            arrayfun(@(x) x.prepareSave, nw.Templates)
            fname = 'tmp.mat';
            ffile = fullfile(save_path,fname);
            save(ffile,'nw');
        end
        
        %% Plots
        function f = plotNetworkActivity(obj,t_start,length,tmp_ID)
            if nargin < 3
                t_start = 0;
                length = obj.Duration;
            end
            t_end = t_start+length;
            spk = obj.Active(1,:);
            ch = obj.Active(2,:);
            act = spk>=t_start & spk<=t_end;
            ch = ch(act);
            ids = unique(ch);
            Z=[];
            for i = 1:numel(ids)
               Z(i) = sum(ids(i)==ch); %Fix
               el_id = obj.Templates(ids(i)).Electrodes;
               X = obj.XYElectrodes(el_id,1);
               Y = obj.XYElectrodes(el_id,2);
               X = mean(X);
               Y = mean(Y);
               x(i) = X;
               y(i) = Y;
            end
            if nargin==4
            for i = 1:numel(tmp_ID)
                ref_el = obj.Templates(tmp_ID(i)).RefElectrode;
                X_tmp = obj.XYElectrodes(ref_el,1);
                Y_tmp = obj.XYElectrodes(ref_el,2);
                x_tmp(i) = X_tmp;
                y_tmp(i) = Y_tmp;
            end
            else
                x_tmp = [];
                y_tmp = [];
            end
            Z = Z/length;
            num_el = {obj.Templates.Electrodes};
            num_el = cellfun('length',num_el(ids));
            sz = 4*num_el;
            f = scatter(x,y,sz,Z,'filled');
            hold on
            scatter(x_tmp,y_tmp,100,'k*')
            f.MarkerEdgeAlpha = 0.8;
            f.MarkerFaceAlpha = 0.5;

            c = colorbar;
            ll = ceil(abs(log10(min(Z))));
            c.Ticks = [10.^(-ll:0) ceil(max(Z))];
            c.Label.String = 'Frequency per Unit [Hz]';
            set(gca,'ColorScale','log','CLim',[10^-ll ceil(max(Z))])
        end
        

        
        function p = plotCoactivity(obj,binning,bursts)
            
            if nargin == 3
                N_bursts = length(obj.Bursts.T_start);
                if max(bursts)>N_bursts
                    warning('Your specified number of bursts is larger than the actual burst count. Using all available bursts instead.')
                    bursts = 1:N_bursts;
                end
                t_start = obj.Bursts.T_start(bursts(1))-1;
                t_end = obj.Bursts.T_end(bursts(end))+1;
                for i = bursts(1):bursts(end)
                    xline(obj.Bursts.T_start(i),'--b',{string(round(obj.Bursts.T_start(i)))},'LabelHorizontalAlignment','center');
                    hold on
                    xline(obj.Bursts.T_end(i),'--r',{string(round(obj.Bursts.T_end(i)))},'LabelHorizontalAlignment','center');
                end
            else
                t_start = 0;
                t_end = 180;
            end
            y = histc(obj.Active(1,:),t_start:binning:t_end);
            x = t_start:binning:t_end;
            y(1) = 0;
            y(end) = 0;
            colormap('gray')
            p = patch(x,y,y,'EdgeColor','interp','FaceColor','interp');
            p.CData = abs(p.CData-max(p.CData));
            h=gca;
            h.CLim = [max(y)/3 max(y)];
            ylabel('Coactivity');
            xlabel('Time [s]');
            axis tight
        end
        
        function plotSpikePerUnit(obj)
            ch = obj.Active(2,:);
            [C,~,ic] = unique(ch);
            a_counts = accumarray(ic,1);
            %             subplot(1,2,1)
            bar(C,a_counts)
            axis tight
            xlabel('Units');
            ylabel('# Spikes')
        end
        
        function plotWaveforms(obj)
            wfs = vertcat(obj.Templates([obj.ActiveChannels]).MaxWf)*6.2;
            plot(wfs')
            set(gca, 'Box', 'off','TickDir', 'out','TickLength', [.01 .01] ,'XMinorTick','off','YMinorTick', 'off','LineWidth', 1);
            axis tight
            ylabel('Amplitude (\muV)');
            xlabel('Time (samples)')
        end
        
        function p = plotActivityRaster(obj,t_start,t_end)
            if nargin == 1
                t_start = 0;
                t_end = obj.Active(1,end);
            end
            
            spk = obj.Active(1,:);
            ind = spk>=t_start & spk<=t_end;
            p = plot(spk(ind),obj.Active(2,ind),'k.','MarkerSize',0.0001);
            p.Color(4) = 0.1;
            ylabel('Units');
            xlabel('Time [s]');
            axis tight
        end

        function plotOverview(obj)
            name = sprintf('Chip: %s Date: %s',obj.ChipID, obj.RecordingDate);
            figure('Name',name,'Units','pixels','Position',[0 0 1900 1080]);
            subplot(3,2,1); obj.plotWaveforms; title('Average Waveforms')
            ax = subplot(3,2,2); colormap(ax,'gray'); obj.plotCoactivity(0.1);  title('Coactivity')
            ax = subplot(3,2,3); obj.plotNetworkActivity; colormap(ax,'parula'); title('Unit activity')
            subplot(3,2,4); obj.plotSpikePerUnit; title('Spikes per Unit')
            subplot(3,2,[5,6]); obj.plotActivityRaster(0,180); title('Activity Raster')
        end
    end
end

