%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing fMRI data - Multiple subjects   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% General
owd = '/imaging/henson/users/sl01/FacesData/';
if ~exist(owd); mkdir(owd); end
cd(owd)
sessions = [3:11];
subs = {'subject_01','subject_02','subject_03','subject_05','subject_06','subject_08','subject_09','subject_10','subject_12','subject_14','subject_15','subject_16','subject_17','subject_18','subject_19','subject_23','subject_24','subject_25'};

%% Realignment
% remove dummy scans
for s = 1:length(subs)
    sub = subs{1,s};
    mkdir([owd,sub,'/bold/unused_scans'])
    cd([sub,'/bold/'])
    for i = 1:length(sessions)
        cd(num2str(sessions(i), '%03d')) 
        scan1 = dir('fMR*00001-*.nii');
        scan2 = dir('fMR*00002-*.nii');
        if ~isempty(scan1)
            movefile([scan1.name],[owd,sub,'/bold/unused_scans'])
        end;
        if ~isempty(scan2)
            movefile(scan2.name,[owd,sub,'/bold/unused_scans'])
        end;
        cd .. 
    end
    cd(owd)
end;

% realign
for s = 1:length(subs)
    sub = subs(s);
    P = cell(length(sessions),1);
    for i=1:length(sessions)
        P{i,:} = spm_select('FPList',fullfile(owd,sub,'/bold',num2str(sessions(i), '%03d')),'^f.*nii');
    end
    flags = struct('rtn',1);
    %flags.rtm = 1;
    spm_realign(P, flags);
    spm_reslice(P, flags);
end;
sub='subject_01';

%% Slice-Timing Correction
for s=1:length(subs)
    sub = subs(s);
    PF = cell(length(sessions),1);
    for i=1:length(sessions)
        PF{i,:} = spm_select('FPList',fullfile(owd,sub,'/bold',num2str(sessions(i), '%03d')),'^rfMR.*nii');
    end
    nslices = 33; TR = 2;
    timing = [TR / nslices, TR / nslices];
    refslice = 2;
    sliceorder = [1:2:33 2:2:32];
    spm_slice_timing(PF, sliceorder, refslice, timing)
end;
sub = 'subject_01';

clearvars -except owd sub sessions subs

%% Coregistration
for s=11:length(subs)
    sub = subs(s);
    [PS,flags] = spm_select('FPList',fullfile(owd,sub),'^mprage.nii');
    [PF,flags] = spm_select('FPList',fullfile(owd,sub,'/bold/003'),'^meanfMR.*nii');
    out = spm_coreg(PS, PF, flags);
    
    M  = inv(spm_matrix(out));    % convert 6-body params to a matrix
    MM = spm_get_space(PS);  % get vox-2-world matrix for structural
    spm_get_space(PS, M*MM); % Update position of structural (in hdr)
end;
sub = 'subject_01';
flags = [];

%% Segmentation
for s = 1:length(subs)
    sub = subs(s);
    % segment
    [PS] = spm_select('FPList',fullfile(owd,sub),'^mprage.nii');
    out = spm_preproc(PS);

    % convert output to normalisation params
    [sn,isn] = spm_prep2sn(out); 
    save(sprintf('%s_seg_sn.mat',spm_str_manip(PS,'sd')),'sn')
    save(sprintf('%s_seg_inv_sn.mat',spm_str_manip(PS,'sd')),'isn')

    % write converted data
    opts.biascor = 1;
    opts.GM  = [1 0 1];	
    opts.WM  = [1 0 0];
    opts.CSF = [1 0 0];
    opts.cleanup = 1;
    spm_preproc_write(sn,opts);

%% Normalisation and Smoothing 
    PS = spm_select('FPList',fullfile(owd,sub),'^mmprage.nii'); 
    PF = spm_select('FPList',fullfile(owd,sub,'/bold/003'),'^meanfMR.*nii');

    % normalise structural image
    flags.interp = 7;
    flags.wrap = [0 1 0];
    flags.vox = [1 1 1]; 
    spm_write_sn(PS,sn,flags);

    % normalise functional mean
    flags.vox = [3 3 3];
    spm_write_sn(PF,sn,flags);

    % normalise and smooth all images
    sessions = [3:11];
    parfor sess=1:length(sessions)
        PF = spm_select('FPList',fullfile(owd,sub,'/bold',num2str(sessions(sess), '%03d')),'^arf.*nii');
        for m=1:size(PF,1)
            VO = spm_write_sn(PF(m,:),sn,flags); % or just use spm_write_sn(PF(m,:),sn,flags), to keep warf files (for Hunar)
            [pth,nam,ext] = fileparts(PF(m,:));
            spm_smooth(VO,fullfile(pth,['sw' nam ext]),[8 8 8]);
        end   
    end
end;

clearvars -except owd sub sessions subs

%% Extract trigger codes
for s = 1:length(subs)
    sub = subs{1,s};
    cd(owd)
    sessions = [3:11];
    cd([owd,sub, '/bold/'])
    for i=1:length(sessions)
        fname = [owd,sub,'/bold/',num2str(sessions(i), '%03d'),'/',sub,'_',num2str(i),'.txt_output.txt'];
        trial_onsets{1,i} = get_ons(fname);
    end;

    save trial_onsets trial_onsets 
end;

%% 1st level analysis (using batch)
deb = 1; %1: corrected behav data
if deb == 1
    stat_pth = 'stat_deb';
    ons = 'trial_onsets_deb';
else
    stat_pth = 'stat';
    ons = 'trial_onsets';
end

%% Concatenate Sessions
for s = 1:length(subs)  
    sub = subs{1,s};
    cd(owd)
    mkdir([owd,sub,'/bold/stat_concat_deb'])
    conds = {'F_Init','F_Im','F_L','U_Init','U_Im','U_L','S_Init','S_Im','S_L'};
    load([owd,sub,'/bold/trial_onsets_deb']);

    % create and estimate design matric (can use spm_fMRI_design instead of batch)
    matlabbatch{1}.spm.stats.fmri_spec.dir = {[owd,sub,'/bold/stat_concat_deb']};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    onsets = cell(1,length(sessions)); R = [];
    scans = {}; time_so_far = 0;
    for i = 1:length(sessions)
        scans{i} = spm_select('FPList',fullfile(owd,sub,'/bold',num2str(sessions(i), '%03d')),'^swarf.*nii');
        
        rp_file = spm_select('FPList',fullfile(owd,sub,'/bold/',num2str(sessions(i), '%03d')),'^rp.*txt');
        R = [R; load(rp_file)];
        
        for c = 1:9
            onsets{c} = [onsets{c}; trial_onsets{1,i}{1,c}' + time_so_far];
        end
        nscan(i) = length(scans{i});
        time_so_far = time_so_far + nscan(i)*2;
    end
    move_file = fullfile(owd,sub,'bold','rp_all_sess.mat');
    save(move_file,'R');
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(strvcat(scans{:}));
    
    for c = 1:9
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).name = conds{1,c};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).onset = onsets{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(c).orth = 1;
    end;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = cellstr(move_file);
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
    
    % concatenate
    spm_fmri_concatenate(fullfile(owd,sub,'bold','stat_concat_deb','SPM.mat'),nscan)
end;

%% Estimate
est_dir = 'stat_concat_deb'; sessrep = 'none';

for s = 1:length(subs)
    sub = subs{1,s};
    cd([owd,sub,'/bold/',est_dir])
    load SPM
    spm_spm(SPM)
end;

clearvars -except owd sub sessions subs concat sessrep

%% Define Contrasts
contrastsF.names = {'effects of interest',...
                    'effects vs baseline',...
                    'main effect - faces',...
                    'main effect - repetition'};
contrastsF.weights = {detrend(eye(9),0)',...
                      eye(9)',... 
                      [1 0.5 0.5 -1 -0.5 -0.5 0 0 0; 1 0.5 0.5 0 0 0 -1 -0.5 -0.5 ;0 0 0 1 0.5 0.5 -1 -0.5 -0.5],...
                      [1 -1 0 1 -1 0 1 -1 0; 1 0 -1 1 0 -1 1 0 -1; 0 1 -1 0 1 -1 0 1 -1]};
                
contrastsT.names = {'F_Init',...
                    'F_Im',...
                    'F_L',...
                    'U_Init',...
                    'U_Im',...
                    'U_L',...
                    'S_Init',...
                    'S_Im',...
                    'S_L',...
                    'faces (F+U) > scrambled (S)',...
                    'main effects of RS of faces'};
contrastsT.weights = {[1 0 0 0 0 0 0 0 0],...
                      [0 1 0 0 0 0 0 0 0],...
                      [0 0 1 0 0 0 0 0 0],...
                      [0 0 0 1 0 0 0 0 0],...
                      [0 0 0 0 1 0 0 0 0],...
                      [0 0 0 0 0 1 0 0 0],...
                      [0 0 0 0 0 0 1 0 0],...
                      [0 0 0 0 0 0 0 1 0],...
                      [0 0 0 0 0 0 0 0 1],...
                      [0.25 0.125 0.125 0.25 0.125 0.125 -0.5 -0.25 -0.25],...
                      [0.5 -0.25 -0.25 0.5 -0.25 -0.25 0 0 0]};
                  
for s = 1:length(subs)
    sub = subs{1,s};
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_select('FPList',fullfile(owd,sub,'/bold/',est_dir),'^SPM.mat'));
    for i =1:length(contrastsF.names) % create f contrasts
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = contrastsF.names{1,i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.convec = contrastsF.weights{1,i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.sessrep = sessrep;
    end;
    for i = 1:length(contrastsT.names) % create t contrasts
        matlabbatch{1}.spm.stats.con.consess{i+length(contrastsF.names)}.tcon.name = contrastsT.names{1,i};
        matlabbatch{1}.spm.stats.con.consess{i+length(contrastsF.names)}.tcon.convec = contrastsT.weights{1,i};
        matlabbatch{1}.spm.stats.con.consess{i+length(contrastsF.names)}.tcon.sessrep = sessrep;
    end;
    matlabbatch{1}.spm.stats.con.delete = 1;
    spm('defaults', 'FMRI');
    spm_jobman('run',matlabbatch);
end;
clear matlabbatch

% display contrasts
figure
for s = 1:length(subs) 
    sub = subs{1,s};
    cd([owd,sub,'/bold/stat'])
    
    subplot(4,5,s)
    V=spm_vol('ess_0008.nii'); % select contast
    Y=spm_read_vols(V);
    YF = Y;
    YF(isnan(YF)) = 30;
    imagesc(squeeze(YF(:,:,40))); % display the wanted axial slice
    colormap('bone')
    title(sub)
    set(gca, 'Clim', [-20 20])
end

%% 2nd level analysis
matlabbatch{1}.spm.stats.factorial_design.dir = {[owd,'GroupStats/Univariate_concat_deb']};
for s = 1:length(subs)
    contrasts = cellstr(spm_select('FPList',fullfile(owd,subs(s),'/bold/stat_concat_deb'),'^con.*nii'));
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).scans = contrasts(1:9,1); % All 9 conditions
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).conds = [1 2 3 4 5 6 7 8 9];
end;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);

%% Estimate design matrix
cd([owd,'/GroupStats/Univariate_concat_deb'])
load SPM
spm_spm(SPM) % estimate

%% Define contrasts
clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_select('FPList',fullfile(owd,'/GroupStats/Univariate_concat_deb'),'^SPM.*mat'));

contrastsF.names = {'effects of interest',...
                    'main effect - faces',...
                    'main effect - repetition'};
contrastsF.weights = {detrend(eye(9),0)',...
                      [1 0.5 0.5 -1 -0.5 -0.5 0 0 0; 1 0.5 0.5 0 0 0 -1 -0.5 -0.5 ;0 0 0 1 0.5 0.5 -1 -0.5 -0.5],...
                      [1 -1 0 1 -1 0 1 -1 0; 1 0 -1 1 0 -1 1 0 -1; 0 1 -1 0 1 -1 0 1 -1]};
                  
contrastsT.names = {'faces (F+U) > scrambled (S)',...
                    'main effects of RS of faces'};
contrastsT.weights = {[0.25 0.125 0.125 0.25 0.125 0.125 -0.5 -0.25 -0.25],...
                      [0.5 -0.25 -0.25 0.5 -0.25 -0.25 0 0 0]};
                  
for i =1:length(contrastsF.names) % create f contrasts
    matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = contrastsF.names{1,i};
    matlabbatch{1}.spm.stats.con.consess{i}.fcon.convec = contrastsF.weights{1,i};
    matlabbatch{1}.spm.stats.con.consess{i}.fcon.sessrep = sessrep;
end;

for i = 1:length(contrastsT.names) % create t contrasts
    matlabbatch{1}.spm.stats.con.consess{i+length(contrastsF.names)}.tcon.name = contrastsT.names{1,i};
    matlabbatch{1}.spm.stats.con.consess{i+length(contrastsF.names)}.tcon.convec = contrastsT.weights{1,i};
    matlabbatch{1}.spm.stats.con.consess{i+length(contrastsF.names)}.tcon.sessrep = sessrep;
end;

matlabbatch{1}.spm.stats.con.delete = 1;
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);

clear matlabbatch

%% Display contrasts
matlabbatch{1}.spm.stats.results.spmmat = cellstr(fullfile(owd,sub,'/bold/stat_concat_deb/SPM.mat'));
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = true;
spm('defaults', 'FMRI');
spm_jobman('run',matlabbatch);