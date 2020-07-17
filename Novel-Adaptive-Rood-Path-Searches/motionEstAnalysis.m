clear
clc
mbSize = 16;
p = 15;

Final_result(1,1) = {'Name'};
Final_result(1,2) = {'PSNR'};
Final_result(1,3) = {'Time(sec)'};

datt = {'Foreman_CIF\','Tennis_CIF\','Garden_CIF\'}; % CIF
% datt = {'City\','Hall\','Container\'}; % QCIF
for iiiiiii = 1:size(datt,2)
dat_path = datt{iiiiiii};
dat_path = ['..\VideoSeq\' dat_path];
dat = dir([dat_path '*.png']);

NumOfFrame = size(dat,1); % Number of video frames
clearvars Seq Seq_r MVpre1 MVpre2 motionVect imgComp MSE computations SSIM_out GBIQA_out;


NOF = min(300,NumOfFrame);
for K = 1:NOF %length(dat)
    
    filename2 = strcat(dat_path,dat(K).name(1:end));
    Seq(:,:,:,K) = double(rgb2ycbcr(imread(filename2))); 
%     Seq_o(:,:,:,K) = double(imread(filename2)); 
end
Seq_r(:,:,:,1) = Seq(:,:,:,1); 
[row,col,c] = size(Seq(:,:,:,1)); 
MVpre1 = zeros(2,row*col/mbSize^2);
MVpre2 = zeros(2,row*col/mbSize^2);
MVpre3 = zeros(2,row*col/mbSize^2);
tic
for K = 2:NOF
    kk = K-1;
    imgI = Seq(:,:,:,kk);
    imgP = Seq(:,:,:,K);
    [sa,sb] = size(imgP(:,:,1));

%% Exhaustive Search
    [motionVect, computationsES, Rooding] = motionEstES(imgP(:,:,1),imgI(:,:,1),mbSize,p);
    imgComp = motionComp(imgI, motionVect, mbSize);

    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue only Y
    MSE_ES(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).
    computations_ES(kk) = computationsES;

%% Diamond Search
    [motionVect, computat] = motionEstDS(imgP(:,:,1),imgI(:,:,1),mbSize,p);
    computations_DS(kk) = computat;
%     Seq_r(:,:,:,K) = imgComp;
    imgComp = motionComp(imgI, motionVect, mbSize);
    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue
    MSE_DS(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).

%%    Adaptive Rood Patern Search with ZMP
    [motionVect, ZMPcomputat] = motionEstARPSZMP(imgP(:,:,1),imgI(:,:,1),mbSize,p);
    computations_ARPS(kk) = ZMPcomputat;
    imgComp = motionComp(imgI, motionVect, mbSize);
%     Seq_r(:,:,:,K) = imgComp;
    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue
    MSE_ARPS(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).
%     computations(kk) = computat;
% %     
%% FDGDS (2010)
    [motionVect, FDGDScomputations] = MEFDGDS(imgP(:,:,1), imgI(:,:,1), mbSize, p);
    imgComp = motionComp(imgI, motionVect, mbSize);
    computations_FDGDS(kk) = FDGDScomputations;
    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue
    MSE_FDGDS(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).

%% MDPS (2013)
    [motionVect, MDPScomputations] = motionEstDPS(imgP(:,:,1), imgI(:,:,1), mbSize, p, MVpre1);
    computations_MDPS(kk) = MDPScomputations;
    MVpre1 = motionVect;
    imgComp = motionComp(imgI, motionVect, mbSize);
    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue only Y
    MSE_MDPS(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).
    
%% DSP- MV (2016)
    [motionVect, DPSMVcomputations] = motionEstDSPMV(imgP(:,:,1), imgI(:,:,1), mbSize, p, MVpre2);
    computations_DPSMV(kk) = DPSMVcomputations;
    MVpre2 = motionVect;
    imgComp = motionComp(imgI, motionVect, mbSize);
%     Seq_r(:,:,:,K) = imgComp;
    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue
    MSE_DPSMV(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).
%% Proposed CHEN
    [motionVect, Chencomputations] = MEChenARPS2(imgP(:,:,1), imgI(:,:,1), mbSize, p, MVpre3);
    computations_Pro(kk) = Chencomputations;
    MVpre3 = motionVect;
    
    imgComp = motionComp(imgI, motionVect, mbSize);
    Seq_r(:,:,:,K) = imgComp;
    Res = imgComp(:,:,1)-imgP(:,:,1); % Residue only Y
    MSE_Pro(kk)  = norm(Res(:), 'fro')^2/numel(Res); % sqrt(sum(diag(X'*X))).

end
t = toc;
%% Ploting  
string1 = [dat_path(1:end-5),' CIF MSE'];
% MSE
figure;
plot(MSE_ES, 'k--o');

hold on;
plot(MSE_DS, 'r-+');
plot(MSE_ARPS, 'y-x');
plot(MSE_FDGDS, 'c->');
plot(MSE_MDPS, 'g-diamond'); 
plot(MSE_DPSMV, 'm-^');
plot(MSE_Pro, 'b-pentagram');
hold off;

xlabel('Frame number');
ylabel('M.S.E');
set(gca,'fontsize',20)
h = legend('ES','DS[18]','ARPS[19]','FDGDS[38]','MDPS[37]','DPSMV[26]','Proposed',1);
set(h,'Fontsize',25);       
title(string1,'fontsize',25,'fontweight','bold');

% % PSNR
figure;
plot(computations_DS, 'r-+');
hold on;
plot(computations_ARPS, 'y-x');
plot(computations_FDGDS, 'c->');
plot(computations_MDPS, 'g-diamond'); 
plot(computations_DPSMV, 'm-^');
plot(computations_Pro, 'b-pentagram');
hold off;

xlabel('Frame number');
ylabel('Search Points / Bblock)');
set(gca,'fontsize',20)
h = legend('DS[18]','ARPS[19]','FDGDS[38]','MDPS[37]','DPSMV[26]','Proposed',1);
set(h,'Fontsize',25);       
% title(string2,'fontsize',25,'fontweight','bold');
Final_result
end