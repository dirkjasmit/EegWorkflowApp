function REST_EEG = eeg_REST_reref(EEG, chanlist, leadfieldfile)

% wrapper to run the REST rereference (infinite ref) on EEG data. This is
% an EEGLAB plugin that cannot be called from the command line, so this
% wrapper was made.

if nargin<3
    leadfieldfile = '';
end
if nargin<2 || isempty(chanlist)
    chanlist = 1:EEG.nbchan;
end

if isempty(EEG.data)
    errordlg('EEG is empty!!!!','Data Error');
    return
end

% ----------------------------------------------
% Load lead field --> G

if ~isempty(leadfieldfile)
    disp(sprintf('Loading Lead Field... %s',leadfieldfile));

    [~, ~, ext] = fileparts(leadfieldfile);
    switch ext
        case '.xlsx'
            G = xlsread(leadfieldfile);
        case '.xls'
            G = xlsread(leadfieldfile);
        case '.dat'
            G = load(leadfieldfile);
        case '.txt'
            G = load(leadfieldfile,'-ascii');
    end
    if sum(isnan(G(:)))>0
        errordlg('NaN is contained in lead field matrix!!!!','Data Error');
        return
    end
    
else
    % calculating leadfield at once
    disp('Calculating leadfield based on 3-concentric spheres headmodel at once...');

    % -------------
    % use xyz coordinates in the EEG.chanlocs.
    if isfield(EEG.chanlocs,'X') && isfield(EEG.chanlocs,'Y') && isfield(EEG.chanlocs,'Z')
        if ~isempty(EEG.chanlocs(1).X) && ~isempty(EEG.chanlocs(1).Y) && ~isempty(EEG.chanlocs(1).Z) && std([EEG.chanlocs.Z])
            channs = chanlist; % selected channs
            xyz_elec = zeros(length(channs),3);
            for nc = 1:length(channs)
                xyz_elec(nc,1) = EEG.chanlocs(channs(nc)).X;
                xyz_elec(nc,2) = EEG.chanlocs(channs(nc)).Y;
                xyz_elec(nc,3) = EEG.chanlocs(channs(nc)).Z;
            end
        else
            errordlg('EEG coordinates (EEG.chanlocs.X/Y/Z) are empty / zeroed, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
            return
        end
    else
        errordlg('EEG coordinates (EEG.chanlocs.X/Y/Z) are empty, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
        return
    end
    
    % -------------------
    % load fixed dipoles and define their oritations. It can be defined by
    % a file with dipole coordinates.
    [ProgramPath, ~, ~] = fileparts(which('eeg_REST_reref.m'));
    xyz_dipoles = load([ProgramPath,filesep,'corti869-3000dipoles.dat']);
    
    % Calculate the dipole orientations.
    xyz_dipOri = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
    xyz_dipOri ( 2601: 3000, 1 ) = 0;
    xyz_dipOri ( 2601: 3000, 2 ) = 0;
    xyz_dipOri ( 2601: 3000, 3 ) = 1;
    
    % ------------------
    % define headmodel
    headmodel        = [];
    headmodel.type   = 'concentricspheres';
    headmodel.o      = [ 0.0000 0.0000 0.0000 ];
    headmodel.r      = [ 0.8700,0.9200,1];
    headmodel.cond   = [ 1.0000,0.0125,1];
    headmodel.tissue = { 'brain' 'skull' 'scalp' };
    
    % -------------------
    % calculate leadfield
    [G,~] = dong_calc_leadfield3(xyz_elec,xyz_dipoles,xyz_dipOri,headmodel);
    G = G';
end

disp(['Lead Field Matrix: ',num2str(size(G,1)),' sources X ',num2str(size(G,2)),' channels']);

% ----------------------------------------------
% Load EEG data
disp('Loading EEG data...');
try disp(['Current data set: ',EEG.setname]);catch;end;

if isempty(chanlist)
    chanlist = 1:EEG.nbchan;
end
if iscell(chanlist)
    chanlist = find(ismember(chanlist,{EEG.chanlist.labels}));
end

channs = chanlist;
if length(size(EEG.data)) == 3
    OrigData = EEG.data(channs,:);
    disp('********EEG.data is 3D epoched data!!!! Default of data demension is channels X timepoints X epochs!!!');
    disp('********Reshape to channels X timepoints');
else
    OrigData = EEG.data(channs,:);
end
disp(['EEG data: ',num2str(size(OrigData,1)),' channels X ',num2str(size(OrigData,2)),' time points'])

% ----------------------------------------------
if size(OrigData,1) == size(G,2)
    disp('Start from average reference...');
    OrigData = OrigData - repmat(mean(OrigData),size(OrigData,1),1);
    % refer to REST
    disp('Re-referencing to REST...');
    REST_EEG = EEG;
    REST_EEG.data = rest_refer(OrigData,G);
    if length(size(EEG.data)) == 3
        REST_EEG.data = reshape(REST_EEG.data,size(REST_EEG.data,1),size(EEG.data,2),size(EEG.data,3));
        disp('********Reshape to channels X timepoints X epochs!!!!');
    end
    disp('Completed...');
else
    errordlg('No. of Channels of lead field matrix and data are NOT equal!!!','Data Error');
    return;
end



