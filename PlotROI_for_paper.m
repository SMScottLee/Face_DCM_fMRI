Nvoi = 6;
cols = [1:9]; %9con
betas = []; dcm_betas = []; pvar = [];

for s = 1:length(dosubs)
    for v = 1:Nvoi
        sub = dosubs(s);
        swd = sprintf('subject_%02d',sub);
        
        %% Get GLM betas (Might already have betas from above!)
        load(fullfile(rwd,swd,'bold',stat_dir,'SPM.mat'));
        VOI = load(fullfile(owd,swd,stat_dir,VOInames{v}));
        tmp = pinv(SPM.xX.xKXs.X)*VOI.Y;
        betas(s,v,:) = tmp(cols);
        mbeta = mean(squeeze(betas(:,v,:)))';
        mbeta = mbeta/mean(mbeta);
        if v == 1 rEVC(1,1:9) = mbeta;
        elseif v == 2 rOFA(1,1:9) = mbeta;
        elseif v == 3 rFFA(1,1:9) = mbeta;
        elseif v == 4 lEVC(1,1:9) = mbeta;
        elseif v == 5 lOFA(1,1:9) = mbeta;
        elseif v == 6 lFFA(1,1:9) = mbeta;
        end
                
    end
end

PEB_dir = {'Right2';'Right3';'Bilateral6';'Left2';'Left3'};
for N = 1:5
    if N == 1 | N == 4
        load(fullfile(rwd,'PEB',PEB_dir{N},'GCM_peb_fitABC.mat'));
    elseif N == 2 | N == 3 | N == 5
        load(fullfile(rwd,'PEB',PEB_dir{N},'GCM_peb_fitAB.mat'));
    end
    if N == 1 | N == 4 Nvoi = 2;
    elseif N == 2 || N == 5 Nvoi = 3;
    elseif N == 3 Nvoi = 6;
    end
        for s = 1:length(dosubs)
            for v = 1:Nvoi
                sub = dosubs(s);
                swd = sprintf('subject_%02d',sub);
                load(fullfile(rwd,swd,'bold',stat_dir,'SPM.mat'));
                %% Get DCM betas
                Y2 = GCM{s,1}.y(:,v); % fitted data
                Y2 = Y2/GCM{s,1}.Y.scale;
                tmp = pinv(SPM.xX.xKXs.X)*Y2;
                dcm_betas(s,v,:) = tmp(cols);
            end
        end
        if N == 1 
            mdcm_beta = mean(squeeze(dcm_betas(:,1,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rOFA(2,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,2,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rFFA(2,:) = mdcm_beta;
        elseif N == 2 
            mdcm_beta = mean(squeeze(dcm_betas(:,1,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rEVC(3,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,2,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rOFA(3,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,3,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rFFA(3,:) = mdcm_beta;
        elseif N == 3 
            mdcm_beta = mean(squeeze(dcm_betas(:,1,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rEVC(4,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,2,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rOFA(4,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,3,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            rFFA(4,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,4,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lEVC(4,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,5,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lOFA(4,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,6,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lFFA(4,:) = mdcm_beta;
        elseif N == 4 
            mdcm_beta = mean(squeeze(dcm_betas(:,1,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lOFA(2,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,2,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lFFA(2,:) = mdcm_beta;
        elseif N == 5 
            mdcm_beta = mean(squeeze(dcm_betas(:,1,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lEVC(3,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,2,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lOFA(3,:) = mdcm_beta;
            mdcm_beta = mean(squeeze(dcm_betas(:,3,:)))';
            mdcm_beta = mdcm_beta/mean(mdcm_beta);
            lFFA(3,:) = mdcm_beta;
        end    
end

figure; b = bar(rEVC'); 
b(1).FaceColor='b'; b(2).FaceColor='r'; b(3).FaceColor=[1,0.3,0.3]; b(4).FaceColor=[1,0.6,0.6];
legend('real data','fitted data in 2-ROI network','fitted data in 3-ROI network','fitted data in 6-ROI network');
title('rEVC');
xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});

figure; b = bar(rOFA'); 
b(1).FaceColor='b'; b(2).FaceColor='r'; b(3).FaceColor=[1,0.3,0.3]; b(4).FaceColor=[1,0.6,0.6];
legend('real data','fitted data in 2-ROI network','fitted data in 3-ROI network','fitted data in 6-ROI network');
title('rOFA');
xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});

figure; b = bar(rFFA'); 
b(1).FaceColor='b'; b(2).FaceColor='r'; b(3).FaceColor=[1,0.3,0.3]; b(4).FaceColor=[1,0.6,0.6];
legend('real data','fitted data in 2-ROI network','fitted data in 3-ROI network','fitted data in 6-ROI network');
title('rFFA');
xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});

figure; b = bar(lEVC'); 
b(1).FaceColor='b'; b(2).FaceColor='r'; b(3).FaceColor=[1,0.3,0.3]; b(4).FaceColor=[1,0.6,0.6];
legend('real data','fitted data in 2-ROI network','fitted data in 3-ROI network','fitted data in 6-ROI network');
title('lEVC');
xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});

figure; b = bar(lOFA'); 
b(1).FaceColor='b'; b(2).FaceColor='r'; b(3).FaceColor=[1,0.3,0.3]; b(4).FaceColor=[1,0.6,0.6];
legend('real data','fitted data in 2-ROI network','fitted data in 3-ROI network','fitted data in 6-ROI network');
title('lOFA');
xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});

figure; b = bar(lFFA'); 
b(1).FaceColor='b'; b(2).FaceColor='r'; b(3).FaceColor=[1,0.3,0.3]; b(4).FaceColor=[1,0.6,0.6];
legend('real data','fitted data in 2-ROI network','fitted data in 3-ROI network','fitted data in 6-ROI network');
title('lFFA');
xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});


cw = [
     1/6 1/6 1/6         1/6 1/6 1/6         -1/3 -1/3 -1/3; % Perception - faces > scrambled
     1 -1 0              1 -1 0               1 -1 0; % Imm Rep
     1 0 -1              1 0 -1               1 0 -1; % Del Rep
     1/3 1/3 1/3         -1/3 -1/3 -1/3       0 0 0; % Recognition - famous > unfamous     
     1 -1/2 -1/2         1 -1/2 -1/2          0 0 0; % Imm vs Del Rep
     1/2 -1/2 0          1/2 -1/2 0           -1 1 0; % Perception X Imm Rep
     1/2 0 -1/2          1/2 0 -1/2           -1 0 1; % Perception X Del Rep
     1 -1 0              -1 1 0                0 0 0; % Recognition X Imm Rep
     1 0 -1              -1 0 1                0 0 0; % Recognition X Del Rep
    ];

T=[]; p=[]; figure; clf
sp = [2 4 6 1 3 5];
for v=1:Nvoi
    for c=1:size(cw,1)
         dc(:,c) = squeeze(betas(:,v,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    if Nvoi == 2 && H == 1
        VOItitle = VOInames{v+1};
    elseif Nvoi == 3 && H == 1
        VOItitle = VOInames{v};
    elseif Nvoi == 2 && H == 2
        VOItitle = VOInames{v+4};
    elseif Nvoi == 3 && H == 2
        VOItitle = VOInames{v+3};
    elseif Nvoi == 6
        VOItitle = VOInames{v};
    end
    fprintf('%s\n',VOItitle)
    disp([cw T' p'])
    
    for c=1:size(cw,1)
         dc(:,c) = squeeze(dcm_betas(:,v,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    fprintf('DCM: %s\n',VOItitle)
    disp([cw T' p'])
    
    subplot(3,2,sp(v)),
    mbeta = mean(squeeze(betas(:,v,:)))';
    mbeta = mbeta/mean(mbeta)
    mdcm_beta = mean(squeeze(dcm_betas(:,v,:)))';  
    mdcm_beta = mdcm_beta/mean(mdcm_beta)
    bar([mbeta mdcm_beta])
    legend('real data','DCM fitted');
    Title = {'rEVC','rOFA','rFFA','lEVC','lOFA','lFFA'};
    if Nvoi == 2 && H == 1
        title(Title{v+1},'Interpreter','none')
    elseif Nvoi == 3 && H == 1
        title(Title{v},'Interpreter','none')
    elseif Nvoi == 2 && H == 2
        title(Title{v+4},'Interpreter','none')
    elseif Nvoi == 3 && H == 2
        title(Title{v+3},'Interpreter','none')
    elseif Nvoi == 6
        title(Title{v},'Interpreter','none')
    end
    xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});
end
