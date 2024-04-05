%% Read info from sorted neurons extracted with MountainSort and create 1 file per spike.
%-
%- Extracted neurons are saved in I:\Havoc_REC\*session*
%- Include an info matrix for each file include electrode depth and isolation measure for that neuron
%- 
%- step already performed
%-
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY
%- Date: 2023

clear
session = input('Session to convert into .mat files : ','s');
if strcmp(session(1),'G')
   mk = 'GENGER';
else
    mk = 'Havoc_REC'
end
cd(['I:\' mk  '\' session '\']);

list = dir('*firings_curated*.mda');

for i = 1 : length(list)
    
    ch = list(i).name(end-6:end-4);
    
    %- look for the cluster metrics files
    cluster_metrics = dir([list(i).name(1:9) 'clusters_metrics_' list(i).name(end-6:end-4) '.json']);
    [metrics] = textread(cluster_metrics(1).name,'%q');
    
    spikes=readmda(list(i).name); %- load spike times (It's in SAMPLE in this matrix)
    Clust = unique(spikes(3,:));
    
    info.ch = ch;
    info.session = session;
    
    %- load each isolated neuron
    for j = 1 : length(Clust)
        temp = spikes(2,spikes(3,:)==Clust(j));
        
        ClustID = Clust(j);
        timestamps = (temp * (1/40000)); %- convert from sample to sec
        
        info.ClustID = ClustID;
        
        %- extract info from cluster metric file
        label_loc = [find(strcmp(metrics,'label:')==1)+1 ; length(metrics)]; %- look for the different cluster labels
        dumm = label_loc(strcmp(metrics(label_loc),[num2str(ClustID) ','])); % find the location of the one currently loaded
        metric_clust = metrics(dumm:label_loc(find(label_loc==dumm)+1)); %- extract only the value for the currently loaded spike
        
        dumm = metric_clust( find(strcmp(metric_clust,'firing_rate:')==1)+1 ); %- extract Firing rate
        info.fr = str2num(dumm{1}(1:end-1));
        dumm = metric_clust( find(strcmp(metric_clust,'isolation:')==1)+1 ); %- extract Isolation
        info.isolation = str2num(dumm{1}(1:end-1));
        dumm = metric_clust( find(strcmp(metric_clust,'noise_overlap:')==1)+1 ); %- extract Noise Overlap
        info.noise_overlap = str2num(dumm{1}(1:end-1));
        dumm = metric_clust( find(strcmp(metric_clust,'peak_amp:')==1)+1 ); %- extract Peak Amplitude
        info.peak_amp = str2num(dumm{1}(1:end-1));
        dumm = metric_clust( find(strcmp(metric_clust,'peak_noise:')==1)+1 ); %- extract Peak Noise
        info.peak_noise = str2num(dumm{1}(1:end-1));
        dumm = metric_clust( find(strcmp(metric_clust,'peak_snr:')==1)+1 ); %- extract Peak SNR
        info.peak_snr = str2num(dumm{1}(1:end-1));
        dumm1 = metric_clust( find(strcmp(metric_clust,'t1_sec:')==1)+1 ); %- extract tstart
        dumm2 = metric_clust( find(strcmp(metric_clust,'t2_sec:')==1)+1 ); %- extract tend
        info.start_end = [str2num(dumm1{1}(1:end-1)) str2num(dumm2{1}(1:end-1))];
        clear dumm dumm1 dumm2
        
        if ClustID<10 ; Clustname=['0' num2str(ClustID)];
        else Clustname=[num2str(ClustID)];
        end
        
        %- save 1 file per neuron
        save(['I:\\' mk '\SORTED NEURONS\' session '_Ch' ch '_Clus' Clustname '.mat'],'timestamps','info')
        
        clear timestamps
    end
    
end
clear info

%- Extract waveform for every neuron
pre = 30;
post = 50;

cd(['I:\\' mk '\SORTED NEURONS\'])
list_spk = dir(['I:\\' mk '\SORTED NEURONS\' session '*Ch*Clus*'])
rawpath = ['I:\\' mk '\' session '\'];
disp('#############################################################')


Ch = cat(1,list_spk(:).name);
Ch = str2num(Ch(:,12:14));
uniqueCh = unique(Ch);

for c = 1 : length(uniqueCh) %- for each ch with a spike at least
    spks = find(Ch==uniqueCh(c));
    temp = list_spk(spks(1)).name;
    
    rawfile = [rawpath temp(1:8) '_spl_WB' temp(12:14) '.mda'];
    disp(['Extracting Raw signal : ' rawfile])
    raw=readmda(rawfile);
    
    %- filter the raw signal
    raw = double(raw);
    Wp = [ 700 8000] * 2 / 40000; % pass band for filtering
    Ws = [ 500 10000] * 2 / 40000; % transition zone
    [N,Wn] = buttord( Wp, Ws, 3, 20); % determine filter parameters
    [B,A] = butter(N,Wn); % builds filter
    rawfilt = filtfilt( B, A, raw ); % runs filter
    raw = single(rawfilt); clear rawfilt
    
    for s = 1 : length(spks) %- for each spike on a given channel
        spkname = list_spk(spks(s)).name;
        spike = load(spkname);
        
        samples = round(spike.timestamps*40000); %- go back in samples instead of sec
        
        %- remove spikes that are too close to recording boundaries
        while samples(1)-pre<=0
            samples = samples(2:end);
            spike.timestamps = spike.timestamps(2:end);
            disp('Removed First Spike')
        end
        while samples(end)+post>length(raw)
            samples = samples(1:end-1);
            spike.timestamps = spike.timestamps(1:end-1);
            disp('Removed Last Spike')
        end
        
        %- extract waveforms
        waveforms = zeros(length(samples),length(-pre:post));
        for sp = 1 : length(samples)
            waveforms(sp,:) = raw(samples(sp)-pre:samples(sp)+post);
        end
        
        %- update matrix with waveforms now
        info = spike.info;
        timestamps = spike.timestamps;
        save(spkname,'timestamps','info','waveforms')
        
    end
end

clear raw
