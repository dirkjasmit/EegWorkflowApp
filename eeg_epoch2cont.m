function EEG = eeg_epoch2cont(EEG)

if EEG.trials == 1
    return
end

np = EEG.pnts;
nt = EEG.trials;

% channels x pnts x trials  ->  channels x (pnts*trials)
EEG.data = reshape(EEG.data, EEG.nbchan, np*nt);

% boundary between each pair of adjacent epochs
for k = 1:nt-1
    n = numel(EEG.event) + 1;
    EEG.event(n).type     = 'boundary';
    EEG.event(n).latency  = k*np + 0.5;   % half-sample: sits between epochs
    EEG.event(n).duration = 1;
end

if isfield(EEG.event, 'epoch')
    EEG.event = rmfield(EEG.event, 'epoch');
end

EEG.pnts   = np * nt;
EEG.trials = 1;
EEG.xmin   = 0;
EEG.xmax   = (EEG.pnts - 1) / EEG.srate;
EEG.times  = (0:EEG.pnts-1) / EEG.srate * 1000;

EEG.epoch  = [];
EEG.icaact = [];      % force recompute from icaweights
EEG.reject = [];
EEG.stats  = [];

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG);