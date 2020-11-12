classdef Template < handle

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
        Velocity
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
        function tmp = Template(network,id,st,w)
            
            if nargin > 0
                tmp.Network = network;
                tmp.ID = id;
                if ~isempty(network.PrefElectrodes)
                tmp.RefElectrode = network.PrefElectrodes(id);
                else
                    tmp.RefElectrode = [];
                end
                tmp.SpikeFrames= tmp.getSpikes(st);
                tmp.Frequency = tmp.calcFrequency(network);
                tmp.SparseMatrix = w;
                tmp.Electrodes = tmp.findElectrodes();
                tmp.TemplateArea = tmp.calcTemplateArea();
            end
        end
        
        function st = get.SpikeTimes(obj)
           st = double(obj.SpikeFrames(1,:))/20000; 
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
            [X_all,Y_all] = pos2cor(pos(:,1),pos(:,2));
            X_active = X_all(obj.Electrodes);
            Y_active = Y_all(obj.Electrodes);
            X_ref = X_all(obj.RefElectrode);
            Y_ref = Y_all(obj.RefElectrode);
            name = sprintf('Active Electrodes Template %d',obj.ID);
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
            rectangle('Position',[0 0 220 120]);
            ax = gca;
            ax.XLim = [-10 230];
            ax.XTick = 0:20:220;
            ax.YLim = [-10 130];
            ax.YTick = 0:20:120;
            grid minor
        end
        
        function h = plotConnectivity(obj,nw,threshold)
            if nargin ==1
                threshold = 0.1;
            end
            h = figure;
            sd = std([nw.Templates.NormCon]);
            %nw.Synchronicity+20*sd;
            con = obj.NormCon>threshold;
            con_v = [1 obj.NormCon(con)];
            ids = nw.ActiveChannels(con);
            ids = [obj.ID; ids];
            for i = 1:numel(ids)
               Z(i) = length(nw.Templates(ids(i)).Spikes(1,:));
               el_id = nw.Templates(ids(i)).Electrodes;
               X = nw.XYElectrodes(el_id,1);
               Y = nw.XYElectrodes(el_id,2);
               X = mean(X);
               Y = mean(Y);
               [X,Y] = pos2cor(X,Y);
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
%             plot([repmat(cpx,length(x),1) x']', [repmat(cpy,length(y),1) y']')
            xlim([0 220])
            ylim([0 120])
        end
        
        function spt = getSpikes(obj,st)
            ind = st(2,:)==obj.ID;
            spt = uint32(st(1,ind)*20000); %Spike times in sec
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
%             wfs_table.x_inc_adj = temp_x(idx);
%             wfs_table.y_inc_adj = temp_y(idx);
            table_risefall = fit_risefall(wfs_table,table_wfs_features);
%             try
%             [lm,~,~] = fit_speed_3_circus_template(wfs_table);
%             catch
                lm = NaN;
%             end
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
            obj.Velocity = lm;
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
%             [~, rate, regularity] = BayesRR(double(obj.SpikeFrames(1,:)));
%             obj.Rate = mean(log(accumarray(ceil((1:numel(rate))/20)',rate(:),[],@mean)));
%             obj.Regularity = mean(log(accumarray(ceil((1:numel(regularity))/20)',regularity(:),[],@mean)));
        end
                
        
        
        
        
        
        function [mag,freq,exp_fit] = calcRegularity(obj)
           binning = 0.1;
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
                
                colorshift = 10;
                c = othercolor('YlGnBu3',ceil(abs(min(vmin))*y_scale)+colorshift);
                x_l = xy(:,1)-min_x/2;
                x_r = xy(:,1)+min_x/2;
                y_l = xy(:,2)-min_y/2;
                y_r = xy(:,2)+min_y/2;
                x = [x_l' x_r' x_l' x_r'];
                y = [y_l' y_l' y_r' y_r'];
%                 subplot(1,2,2)
                for e = 1:size(xy,1)
                    plot(linspace(xy(e,1)-x_offset,xy(e,1)+x_offset,size(wf,2)),y_scale*wf(e,:)+xy(e,2),'Color',c(ceil(abs(vmin(e))*y_scale)+colorshift,:),'LineWidth',1)
                    
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
                patch(smooth_x(b),smooth_y(b),'b','FaceColor',c(end,:),'FaceAlpha',0.1,'EdgeAlpha',0.1,'EdgeColor',c(end,:))
                colormap(c)
                cb = colorbar;
                cb.TickLabels = {};
                cb.Label.String = 'Amplitude';
                cb.Location = 'southoutside';
                cb.Position = [0.33 0.1 0.33 0.05];
                set(gca,'visible','off')
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
        
        function plotTemplateSpread(obj,rawFile,n_spk,window)
            arguments
                obj (1,1) Template
                rawFile (1,1) = []
                n_spk (1,1) {isnumeric} = 1000
                window (1,1) {isnumeric} = 30
            end
            xy = obj.Network.XYElectrodes(obj.Electrodes,:);
            if isempty(rawFile) 
                wf = obj.Waveform(obj.Electrodes,:);
            else
                %                 spkt = rawFile.getSpikeTimes;
                %                 spk_t = spk_t(spkt.channel==obj.RefElectrode)-min(spkt.frameno);
                %                 wf = arrayfun(@(x) rawFile.getCutoutsOneChannel(x,double(spk_t(n_spk))+5,window,window),...
                %                     obj.RefElectrode,'UniformOutput',0);%
                spk_t = double(obj.SpikeFrames);
                wf = arrayfun(@(x) rawFile.getCutoutsOneChannel(x,40001,window,window),...
                    obj.RefElectrode,'UniformOutput',0);
                wf = [wf{:}];
                %
                %                 spk_t = spkt.frameno-rawFile.getFrameNoAt(1);
                %                 wf = rawFile.getData(double(spk_t(n_spk))-window+1720,2*window);
                wf = wf-mean(wf,1);
                [b1, a1] = butter(3,[300/20000,6000/20000]*2,'bandpass');
% [b1, a1] = butter(3,300/20000*2,'high');
                wf_filt = filtfilt(b1,a1,wf);
                
                figure;plot(wf_filt,'LineWidth',0.5)%linspace(0,2*window,size(wf_filt,1)),
                xlabel('Time [ms]')
                hold on
                yline(-7*mad(wf_filt,0),'r--')
                xline(window)
%                 ylim([-20 10])
                wf = wf_filt';
            end
            [vmin,~] = min(wf,[],2);
            v_idx = vmin~=0;
            vmin = vmin(v_idx);
            xy = xy(v_idx,:);
            wf = wf(v_idx,:);
            min_x = min(abs(diff(unique(xy(:,1)))));
            min_y = min(abs(diff(unique(xy(:,2)))));
            x_offset = (min_x*0.9)/2;
            y_scale = 0.9*(min_y/abs(min(min(wf))))/2;
            
            colorshift = 10;
            c = othercolor('YlGnBu3',ceil(abs(min(vmin))*y_scale)+colorshift);
            x_l = xy(:,1)-min_x/2;
            x_r = xy(:,1)+min_x/2;
            y_l = xy(:,2)-min_y/2;
            y_r = xy(:,2)+min_y/2;
            x = [x_l' x_r' x_l' x_r'];
            y = [y_l' y_l' y_r' y_r'];
            %                 subplot(1,2,2)
            for e = 1:size(xy,1)
                plot(linspace(xy(e,1)-x_offset,xy(e,1)+x_offset,size(wf,2)),y_scale*wf(e,:)+xy(e,2),'Color',c(ceil(abs(vmin(e))*y_scale)+colorshift,:),'LineWidth',1)
                
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
            patch(smooth_x(b),smooth_y(b),'b','FaceColor',c(end,:),'FaceAlpha',0.1,'EdgeAlpha',0.1,'EdgeColor',c(end,:))
            colormap(c)
            cb = colorbar;
            cb.TickLabels = {};
            cb.Label.String = 'Amplitude';
            cb.Location = 'southoutside';
            cb.Position = [0.33 0.1 0.33 0.05];
            set(gca,'visible','off')
        end
        function prepareSave(obj)
           obj.SparseMatrix = [];
           obj.SpikeFrames = [];
           obj.ConVector = [];
        end
    end
end

