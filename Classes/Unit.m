classdef Unit < handle

    properties(Hidden)
       SparseMatrix
       SpikeFrames % (2,:) is Burst indicator with 1=Inburst 0=Outburst
    end
    
    properties(Dependent)
        Waveform
        SpikeTimes
    end
    
    properties
        Network
        Amplitude
        ID
        Electrodes
        RefElectrode
        Frequency
        ConVector
        
        MeanCon
        Halfwidth
        Asymmetry
        T2Pratio
        T2Pdelay
        AUCpeak1
        WidthPeak1
        AUCpeak2
        WidthPeak2
        AUCtrough
        WidthTrough
        Rise
        Decay
        
        ISIMean
        ISIVar
        Fano
        PACF
        ISIMeanInburst
        ISIVarInburst
        FanoInburst
        PACFInburst
        ISIMeanOutburst
        ISIVarOutburst
        FanoOutburst
        PACFOutburst
        InburstRatio
        Burstiness
        MaxWf
        BlDev
        Activity
        
        TemplateArea
        Sholl
        Spread
        CritV
        DendMax
        
        ResonanceFrequency
        ResonanceMagnitude
        ResonanceFit
    end
    

    
    methods
        function tmp = Unit(network,id,st,w)
            
            if nargin > 0
                tmp.Network = network;
                tmp.ID = id;
                tmp.SpikeFrames= tmp.getSpikes(st);
                tmp.Frequency = tmp.calcFrequency(network);
                tmp.SparseMatrix = w;
                tmp.Electrodes = tmp.findElectrodes();
                tmp.TemplateArea = tmp.calcTemplateArea();
            end
        end
        
        function st = get.SpikeTimes(obj)
           st = double(obj.SpikeFrames(1,:))/obj.Network.SamplingRate; 
        end
        
        function wf = get.Waveform(obj)
            wf = full(obj.SparseMatrix);
        end
        
        function getInburstRatio(obj,burst)
            obj.SpikeFrames(2,:) = burst;
            obj.InburstRatio = sum(burst)/length(burst);
        end
                
        function active = findElectrodes(obj)
            wf = obj.Waveform;
            ind = find(wf);
            [r,~] = ind2sub(size(wf),ind);
            active = unique(r)';
            if isempty(obj.RefElectrode)
                [r,~] = find(wf==min(min(wf)));
                obj.RefElectrode = r;
            end
            
        end
        
        function ta = calcTemplateArea(obj)
            pos = obj.Network.XYElectrodes;
            if length(obj.Electrodes)>2
                X_all= pos(:,1);
                Y_all = pos(:,2);
                X_active = X_all(obj.Electrodes);
                Y_active = Y_all(obj.Electrodes);
                b = boundary(X_active,Y_active,1);
                ta = polyarea(X_active(b),Y_active(b));
            else
                ta = 612.5; %Area of 2 electrodes
            end
        end
        
        function plotActiveElectrodes(obj)
            nw = obj.Network;
            pos = nw.XYElectrodes;
%             [X_all,Y_all] = pos2cor(pos(:,1),pos(:,2));
            X_all = pos(:,1);
            Y_all = pos(:,2);
            X_active = X_all(obj.Electrodes);
            Y_active = Y_all(obj.Electrodes);
            X_ref = X_all(obj.RefElectrode);
            Y_ref = Y_all(obj.RefElectrode);
            name = sprintf('Active Electrodes Unit %d',obj.ID);
%             figure('Name',name,'Position',[10 10 880 480]);
            plot(X_all,Y_all,'b.');
            hold on
            b = boundary(X_active,Y_active,1);
            patch(X_active(b),Y_active(b),'g','FaceAlpha',0.5,'LineSmoothing','on')
            plot(X_active,Y_active,'o',...
                'MarkerSize',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g');
            plot(X_ref,Y_ref,'o',...
                'MarkerSize',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r')
%             rectangle('Position',[0 0 220 120]);
%             ax = gca;
%             ax.XLim = [-10 230];
%             ax.XTick = 0:20:220;
%             ax.YLim = [-10 130];
%             ax.YTick = 0:20:120;
%             grid minor
            axis tight
        end
        
        function h = plotConnectivity(obj)
            h = figure;
            nw = obj.Network;
            con = isoutlier(obj.ConVector);
            con_v = [1 obj.ConVector(con)];
            ids = nw.ActiveChannels(con);
            ids = [obj.ID ids];
            for i = 1:numel(ids)
               Z(i) = length(nw.Templates(ids(i)).SpikeTimes(1,:));
               el_id = nw.Templates(ids(i)).Electrodes;
               X = nw.XYElectrodes(el_id,1);
               Y = nw.XYElectrodes(el_id,2);
               X = mean(X);
               Y = mean(Y);
               x(i) = X;
               y(i) = Y;
               
               plot([x(1) x(i)],[y(1) y(i)],'LineWidth',con_v(i)*5,'Color',[0,0,0,0.5],'AlignVertexCenters','on')
               hold on
            end
            Z = Z/nw.Duration;
            num_el = {nw.Templates.Electrodes};
            num_el = cellfun('length',num_el(ids));
            sz = 10*num_el;
            f = scatter(x,y,sz,Z,'filled');
            f.MarkerEdgeAlpha = 0.8;
            f.MarkerFaceAlpha = 0.5;
            colormap(jet);c = colorbar;
            ll = ceil(abs(log10(min(Z))));
            c.Ticks = [10.^(-ll:0) ceil(max(Z))];
            c.Label.String = 'Frequency per Unit [Hz]';
            set(gca,'ColorScale','log','CLim',[10^-ll ceil(max(Z))])
        end
        
        function spt = getSpikes(obj,st)
            ind = st(2,:)==obj.ID;
            spt = uint32(st(1,ind)*obj.Network.SamplingRate); %Spike times in sec
            obj.Activity = length(spt);
        end
        
        function freq = calcFrequency(obj,nw)
            freq = obj.Activity/nw.Duration;
        end        
       
        function nc = normCon(obj,cc,nw)
            id = obj.ID;
            if ~isempty(id)
                act = nw.ActiveChannels;
                act = act(act~=id);
                con = cc(id,act);
                spk = nw.Templates(act);
                sz = size(obj.SpikeTimes,2);
                for i = 1:length(con)
                    nc(i) = con(i)/sqrt(size(spk(i).SpikeTimes,2)*sz);
                end
                obj.ConVector = single(nc);
                obj.MeanCon = mean(nc);
            end
        end
        
        function inferFeatures(obj)
            [m,~] = min(obj.Waveform');
            [~,n] = min(m);
            thiswave = obj.Waveform(n,:)';
            obj.MaxWf = thiswave';
            obj.Amplitude = abs(min(obj.MaxWf));
            table_wfs_features = wave_form_features(thiswave);
            idx = find(sum(obj.Waveform'));
            wfs_table.single_template_spatial = obj.Waveform(idx,:)';
            wfs_table.noise_thresh = 0.5;
            table_risefall = fit_risefall(wfs_table,table_wfs_features);
            obj.Halfwidth = table_wfs_features.HalfWidth;
            obj.Asymmetry = table_wfs_features.assy;
            obj.T2Pratio = table_wfs_features.t2pratio;
            obj.T2Pdelay = table_wfs_features.t2pdelay;
            obj.AUCpeak1 = table_wfs_features.AUCpeak1;
            obj.WidthPeak1 = length(table_wfs_features.WidthPeak1);
            obj.AUCpeak2 = table_wfs_features.AUCpeak2;
            obj.WidthPeak2 = length(table_wfs_features.WidthPeak2);
            obj.AUCtrough = table_wfs_features.AUCtrough;
            obj.WidthTrough = length(table_wfs_features.WidthTrough);
            obj.Rise = table_risefall.lm_rising;
            obj.Decay = table_risefall.lm_decay;
                        
            % Activity-based features
            df = diff(obj.SpikeTimes);
            bl = df(df>0.2);
            dis = 1./bl-obj.Network.BaselineFrequency;
            obj.BlDev = mean(abs(dis));
            obj.Burstiness = sum(df<0.2)/length(df);
            obj.ISIMean = mean(df);
            obj.ISIVar = var(df);
            obj.Fano = var(df)/mean(df);
            pacf = parcorr(df,1);
            obj.PACF = pacf(2);
            
            if ~isempty(obj.Network.Bursts)
            burst_idx = obj.Network.Active(3,obj.Network.Active(2,:) == obj.ID);
            burst_start_idx = find(diff(burst_idx)==1)+1;
            burst_end_idx = find(diff(burst_idx)==-1);
            if length(burst_start_idx) < length(burst_end_idx)
                burst_start_idx = [1 burst_start_idx];
            elseif length(burst_start_idx) > length(burst_end_idx)
                burst_end_idx = [burst_end_idx length(burst_idx)];
            end
            diff_inburst = [];
            diff_outburst = [];
            for i = 1:length(burst_start_idx)
                inburst = obj.SpikeTimes(burst_start_idx(i):burst_end_idx(i));
                diff_inburst = [diff_inburst diff(inburst)];
                if i < length(burst_start_idx)
                    outburst = obj.SpikeTimes(burst_end_idx(i)+1:burst_start_idx(i+1)-1);
                    diff_outburst = [diff_outburst diff(outburst)];
                end
            end
            if length(diff_inburst)>2
                obj.ISIMeanInburst = mean(diff_inburst);
                obj.ISIVarInburst = var(diff_inburst);
                obj.FanoInburst = obj.ISIVarInburst/obj.ISIMeanInburst;
                pacf_inburst = parcorr(diff_inburst,1);
                obj.PACFInburst = pacf_inburst(2);
            end
            if length(diff_outburst)>2
                obj.ISIMeanOutburst = mean(diff_outburst);
                obj.ISIVarOutburst = var(diff_outburst);
                obj.FanoOutburst = obj.ISIVarOutburst/obj.ISIMeanOutburst;
                pacf_outburst = parcorr(diff_outburst,1);
                obj.PACFOutburst = pacf_outburst(2);
            end
            end
            % Multichannel features
            [obj.Sholl,obj.Spread,obj.CritV,obj.DendMax] = obj.ShollAnalysis;
            [obj.ResonanceMagnitude,obj.ResonanceFrequency,obj.ResonanceFit] = obj.calcRegularity;
        end
                
        
        
        
        
        
        function [mag,freq,exp_fit] = calcRegularity(obj)
           binning = obj.Network.Params.Regularity.binning;
           Fs = 1/binning;
           spk = histc(obj.SpikeTimes,0:binning:obj.SpikeTimes(end));
           NFFT = length(spk);
           F = (0 : 1/NFFT : 1/2-1/NFFT)*Fs;
           TEMP = fft(spk,NFFT);
           TEMP(1) = 0;
           if length(abs(TEMP(1:round(NFFT/2))))>2
               [mag,idx] = findpeaks(abs(TEMP(1:round(NFFT/2))),'NPeaks',1,'SortStr','descend');
               freq = F(idx);
               
               binned = abs(TEMP(1:round(NFFT/2)));
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
               if cutoff<4 && length(x)>3
                   cutoff = 4;
               end
               x = x(1:cutoff);
               p = p(1:cutoff);
               log_p = log10(p)-min(log10(p));
               try
                   f = fit(x,log_p','exp1');
                   exp_fit = f.b;
               catch
                   exp_fit = NaN;
               end
           else
               mag = NaN;
               freq = NaN;
               exp_fit = NaN;
           end
        end
        
        function d = calcDistance(obj,ids)
            XY = obj.Network.XYElectrodes(obj.RefElectrode,:);
            refs = [obj.Network.Templates(ids).RefElectrode];
            coor = obj.Network.XYElectrodes(refs,:);
            d = sqrt(sum((XY-coor).^2,2));            
        end
        
        function g = calcSpread(obj,plot_ind)
            xy = obj.Network.XYElectrodes(obj.Electrodes,:);
            wf = obj.Waveform(obj.Electrodes,:);
            [vmin,tmin] = min(wf,[],2);
            v_idx = vmin~=0;
            vmin = vmin(v_idx);
            tmin = tmin(v_idx);
            xy = xy(v_idx,:);
            wf = wf(v_idx,:);
            r = floor(tiedrank(tmin));
            d = squareform(pdist([xy(:,1), xy(:,2)]));
            d(d==0) = NaN; %Avoid zero as minimum
            sources = find(r==min(r))';
            targets = 1:length(vmin);
            targets(sources) = [];
            w = zeros(1,length(targets));
            src = zeros(1,length(targets));
            tgt = zeros(1,length(targets));
            i = 0;
            while ~isempty(targets)
                i = i+1;
                dist_mat = d(sources,targets);
                sorted = sort(reshape(dist_mat,[1,numel(dist_mat)]),'ascend');
                ii = 1;
                [s,t] = find(dist_mat==sorted(ii));
                rd = r(targets(t)) - r(sources(s));
                [~,tind] = min(rd(rd>0));
                while isempty(tind)
                    ii = ii+1;
                    [s,t] = find(dist_mat==sorted(ii));
                    rd = r(targets(t)) - r(sources(s));
                    [~,tind] = min(rd(rd>0));
                end
                w(i) = sorted(ii);
                src(i) = sources(s(tind));
                tgt(i) = targets(t(tind));
                sources = [sources targets(t(tind))];
                targets(targets==targets(t(tind))) =[];
            end
            g = digraph(src,tgt,w);
            %             g.Nodes.X = xy(:,1);
            %             g.Nodes.Y = xy(:,2);
            if plot_ind
%                 figure('color','w');
%                 subplot(1,2,1)
%                 plot(g,'XData',xy(:,1),'YData',xy(:,2),'MarkerSize',abs(vmin));
                min_x = min(abs(diff(unique(xy(:,1)))));
                min_y = min(abs(diff(unique(xy(:,2)))));
                x_offset = (min_x*0.9)/2;
                y_scale = 0.9*(min_y/abs(min(min(wf))))/2;
                c = colormap('hot'); c = c(1:200,:);
                c = flipud(c(round(linspace(1,length(c),ceil(abs(min(vmin))*y_scale))),:));
                x_l = xy(:,1)-min_x/2;
                x_r = xy(:,1)+min_x/2;
                y_l = xy(:,2)-min_y/2;
                y_r = xy(:,2)+min_y/2;
                x = [x_l' x_r' x_l' x_r'];
                y = [y_l' y_l' y_r' y_r'];
                for e = 1:size(xy,1)
                    plot(linspace(xy(e,1)-x_offset,xy(e,1)+x_offset,size(wf,2)),y_scale*wf(e,:)+xy(e,2),'Color',c(ceil(abs(vmin(e))*y_scale),:),'LineWidth',1)
                    
                    hold on
                end
                
                b = boundary(x',y',0.8);
                smooth_window = 9;
                smooth_order = 2;
                interp_factor = 2;
                interp_x = interp(x(b),interp_factor);
                interp_y = interp(y(b),interp_factor);
                smooth_x = sgolayfilt(interp_x,smooth_order,smooth_window);
                smooth_y = sgolayfilt(interp_y,smooth_order,smooth_window);
                b = boundary(smooth_x',smooth_y',0.5);
%                 patch(smooth_x(b),smooth_y(b),'b','FaceColor',c(end,:),'FaceAlpha',0.1,'EdgeAlpha',0.1,'EdgeColor',c(end,:))
                colormap(c)
                cb = colorbar;
                cb.Ticks = [0 1];
                cb.TickLabels = [0,round(min(min(wf))*6.2)];
                cb.Label.String = 'Amplitude';
                cb.Location = 'southoutside';
                cb.Position = [0.33 0.1 0.33 0.02];
                set(gca,'visible','off')
                ref = obj.Network.XYElectrodes(obj.RefElectrode,:);
                xlim([ref(1)-100 ref(1)+100])
                ylim([ref(2)-100 ref(2)+100])
            end
        end
            
            function [sh,sp,crit,dm] = ShollAnalysis(obj)
                if length(obj.Electrodes)>5
                    g = obj.calcSpread(0);
                    P = obj.Network.XYElectrodes(obj.Electrodes,:);
                    C = P(obj.Electrodes==obj.RefElectrode,:);
                    ic = zeros(1,5);
                    d = pdist2(P,C);
                    for k = 1:5
                        R = k*35;
                        itsc = 0;
                        for i = 1:height(g.Edges)
                            P1 = d(g.Edges.EndNodes(i,1));
                            P2 = d(g.Edges.EndNodes(i,2));
                            if sum([P1 P2]<R)==1
                                itsc = itsc+1;
                            end
                        end
                        A = pi*R^2;
                        lg = log10(itsc/A);
                        lgic(k) = lg;
                        ic(k) = itsc;
                    end
                    lm = fitlm(35:35:175,lgic);
                    if ~isinf(lm.Coefficients.Estimate(2))
                        sh = lm.Coefficients.Estimate(2);
                    else
                        sh= [];
                    end
                    sp = sum(g.Edges.Weight);
                    [mic,ind] = max(ic);
                    crit = ind*35;
                    dm = mic;
                else
                    sh = [];
                    sp = [];
                    crit = [];
                    dm = [];
                end
            end
        
        function wfd = WfDecay(obj)
            if length(obj.Electrodes) > 3
                XY = obj.Network.XYElectrodes(obj.Electrodes,:);
                d = sqrt(sum((XY-obj.Network.XYElectrodes(obj.RefElectrode,:)).^2,2));
                wf = permute(obj.Network.TemplateMatrix(obj.ID,:,obj.Electrodes),[3 2 1]);
                amp = min(wf,[],2);
                wfd = fit(d,amp,'exp1');
            else
                wfd = [];
            end
        end
        
        function plotWf(obj)
            wf_mat = obj.Waveform(obj.Electrodes,:)'*6.2;
            min_wf = min(wf_mat);
            max_wf = max(wf_mat);
            peak_wfs = nan(1,size(wf_mat,2));
            for i = 1:length(peak_wfs)
                if abs(min_wf(i))>abs(max_wf(i))
                    peak_wfs(i) = min_wf(i);
                else
                    peak_wfs(i) = max_wf(i);
                end
            end
            c_idx = round(peak_wfs-min(peak_wfs))+1;
            cmap = colormap('hot');
            colors = cmap(round(linspace(1,200,max(c_idx))),:);
            p = plot(wf_mat);
            arrayfun(@(x) set(p(x),'Color',colors(c_idx(x),:)),1:length(c_idx))
            axis tight
            xlabel('Time [samples]')
            ylabel('Amplitude (\muV)')
        end
        
        function prepareSave(obj)
           obj.SparseMatrix = [];
           obj.SpikeFrames = [];
           obj.ConVector = [];
        end
    end
end

