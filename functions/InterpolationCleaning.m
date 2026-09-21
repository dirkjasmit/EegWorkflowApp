function [BadChans, BadEpochs, Mask, R, SD, ImputedData, RealData] = InterpolationCleaning(EEGIn, Rcrit, SDcrit, doDetrend, UseParallel)

% InterpolationCleaning is an automated EEG cleaning script. Chunks into 2s
% epochs and compares the actual data to imputed data by comparing the
% spectra. Probably works best with high density EEG only (>32)
%
% Requires channel information in EEG struct!
%
% Input: EEG struct (EEGLAB)
% Output: 
% - suggested bad channels
% - suggested bad epochs (bad channels not considered). NOTE bad channels
%   are not removed!
% - full mask for epoch/channel combinations marked bad
if nargin<5
    UseParallel = false;
end
if nargin<4
    doDetrend = true;
end
if nargin<3
    SDcrit = 10;
end
if nargin<2
    Rcrit = 0;
end

RealData = EEGIn.data;
%ImputedData = RealData; % preallocation, to be overwritten

tmp = InterpolationReplace(EEGIn, UseParallel);
ImputedData = tmp.data;

%{
% parallel version!! impute each channel using all other channels
if UseParallel
    Temp = {};
    parfor ch=1:EEGIn.nbchan
        Temp{ch} = pop_interp(EEGIn,ch,'spherical');
    end

    for ch=1:EEGIn.nbchan    
        try
            ImputedData(ch,:) = Temp{ch}.data(ch,:);
        catch E
            warning('Interpolation failed, using the original data')
            ImputedData(ch,:) = RealData(ch,:);
        end
    end
    
else
    % impute each channel using all other channels
    for ch=1:EEGIn.nbchan
        Temp = pop_interp(EEGIn,ch,'spherical');
        try
            ImputedData(ch,:) = Temp.data(ch,:);
        catch E
            warning('Interpolation failed, using the original data')
            ImputedData(ch,:) = RealData(ch,:);
        end
    end

end
%}

if EEGIn.trials==1
    % chunk the data in 1 s epochs, then compare two ways: (i) correlation
    % between real and imputed data, and (ii) stdev scaled to median stdev.
    len = EEGIn.srate; % length of 1 s
    nepochs = floor(size(RealData,2)/len); % only full epochs
    FullLen = nepochs*len;
    % cut the data into epochs
    RealData = reshape(RealData(:,1:FullLen),[],len,nepochs);
    ImputedData = reshape(ImputedData(:,1:FullLen),[],len,nepochs);
else
    len = size(EEGIn.data,2);
    nepochs = EEGIn.trials;
    FullLen = nepochs*len;
end


% compute spectra
%[P_imp,fs,AllP_imp] = pfft(ImputedData,EEGIn.srate,hanning(len),.5);
%[P_real,fs,AllP_real] = pfft(RealData,EEGIn.srate,hanning(len),.5);

% calculate R and SD, determine artifact mask (nchans x nepochs)
if doDetrend
    R = cell2mat(arrayfun(@(x)diag(corr(detrend(ImputedData(:,:,x)'),detrend(RealData(:,:,x)'))),1:size(RealData,3),'uni',0));
    SD = cell2mat(arrayfun(@(x)std(detrend(ImputedData(:,:,x)')-detrend(RealData(:,:,x)')),1:size(RealData,3),'uni',0)')';
else
    R = cell2mat(arrayfun(@(x)diag(corr(ImputedData(:,:,x)',RealData(:,:,x)')),1:size(RealData,3),'uni',0));
    SD = cell2mat(arrayfun(@(x)std(ImputedData(:,:,x)'-RealData(:,:,x)'),1:size(RealData,3),'uni',0)')';
end
SD = SD./median(SD(:)); % scaling of SD
Mask = R<Rcrit | SD>SDcrit;
BadChans = (nepochs-sum(Mask,2))<0;
BadEpochs = sum(Mask(~BadChans,:))>0;






