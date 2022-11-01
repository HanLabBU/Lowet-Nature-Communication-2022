clear all

if 0
   %%%%%%%%%%%%%%%% 
allsamp=[];
 Z1=[]; Z2=[];Z3=[];
 BRATE=[];SRATE=[];PRATE=[];sesid=[];
   allangs8=[];  allangs140=[];  allangs40=[];
    allangs8B=[];  allangs140B=[];  allangs40B=[];
 %%%%
stim_freq=[  0]

if stim_freq==1
cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DATA\optoDBS\140\')
cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\140\')
else
  cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DATA\optoDBS\40\')
  cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\40\')
end


ses=dir('*.mat')
mj=0;mr=0;
for ind=[1:length(ses)]
    
    load(ses(ind).name)
    
    
       trial_numb=unique(result.trial_vec);
    if length(find(result.trial_vec==trial_numb(1) )) <4000 & length(trial_numb)>2
 mr=mr+1;
    %%%
     FS=828;
    Fn = FS/2;FB=[ 86.4 88.7];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    LFPg= ((filtfilt(B,A,    result.traces(:,1))));
    result.traces(:,1)= result.traces(:,1)-LFPg;
    %%%
   ind
   tr=0;
    lfp=[];  asA=[];clear aV allV8 allV40 allV140 allS alls1
     for  ne=unique(result.trial_vec)
         tr=tr+1;
  v= result.traces(result.trial_vec==ne,1);%./result.tracesB(result.trial_vec==ne)   ;

      Fn = FS/2;FB=[ 87 88];
                         
                                           
 % v=(v-fastsmooth(v,2300,1,1));
  
       % exponential fitting to remove photobleaching
     [ fitbaseline, coeff]=exp_fit_Fx(v,FS);
  
   v=(v-fitbaseline');
  v=v(1:2500);
%   try
%       vsub= result.lfp.trial{ne}(1,:);%result.resultS{ne}.trace_ws  ;
%   end
      vsub=result.resultS{ne}.trace_ws  ;

        vsub=vsub(1:2500)';
        vsub=zscore(vsub-mean(vsub));
        
               Fn = FS/2;FB=[ 7 9];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    Vsub8= angle(hilbert(filtfilt(B,A,     zscore(vsub))));
            Fn = FS/2;FB=[ 35 45];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    Vsub40= angle(hilbert(filtfilt(B,A,     vsub)));
       Fn = FS/2;FB=[ 135 145];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    Vsub140= angle(hilbert(filtfilt(B,A,     vsub)));
    
     allV8(1:length(Vsub8),tr)=Vsub8;
        allV40(1:length(Vsub40),tr)=Vsub40;
         allV140(1:length(Vsub140),tr)=Vsub140;
       
       %%spikes
          strain= result.resultS{ne}.roaster(1:2500);
          allS(1:length(strain),tr)=strain;
           sid=result.resultS{ne}.spike_idx{1};sid=sid-5;sid(sid<1)=[];
           vect=zeros(1,size(allS,1));vect(sid)=1;
          alls1( :,tr)=vect(1:2500);
          
          %%%%
          lfp.trial{tr}(1,:)= vsub';
aV(:,tr)=v; % raw Vm
  lfp.trial{tr}(2,:)= v'; % raw Vm
  lfp.time{tr}= (1:size(v,1))./828;
    
  spx= result.resultS{ne}.spike_idx{1}   ;
         samp= result.resultS{ne}.spike_amplitude{1};%./mean(result.tracesB(result.trial_vec==ne));
         samp(spx <10)=NaN;
         asA=[asA, samp];
     end
  for id=1:2 %:size(vsig,2)
lfp.label{id}= num2str(id);
  end  
 
allV2(:,mr)= nanmean(aV,2);%
 allV2(:,mr)=allV2(:,mr);%
 
            warning off
  cfg = []; %block_type == cfg.blk
    cfg.method ='wavelet'; %'mvar';
    cfg.output ='fourier';
     cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
     cfg.channel= 'all'%; %chans=cfg.channel;
    cfg.foi= [4:1:160];
     cfg.toi=lfp.time{1}(1:1:end) ;
     cfg.width =5;
    cfg.t_ftimwin =[ones(1,length(cfg.foi))*0.4];
freq2 = ft_freqanalysis(cfg, lfp);

 wavD = abs(squeeze(freq2.fourierspctrm(:,1,:,:)));
wavA = angle(squeeze(freq2.fourierspctrm(:,1,:,:)));
fsel=freq2.freq>=7 & freq2.freq<=9;

%allV8= squeeze(nanmean(wavA(:,fsel,:),2))';
 
 slim=5;
   %% PLV with filtered signals
   allS1= alls1;
   allS1([1:828 1650:size(allS1,1)],:)=0;
%    allangs140= [allangs140; allV140(allS1==1)];
%      allangs40= [allangs40; allV40(allS1==1)];
%        allangs8= [allangs8; allV8(allS1==1)];
      Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
      NT=sum(sum(allS1==1));
         T=   Z.^2;
%     Z= (((1/(NT-1))*((T.*NT-1))));
     if NT>=slim;  allPLVs(:,1,ind)= Z;else; allPLVs(:,1,ind)= NaN;end
       
        allS1= alls1;
   allS1([ 800:size(allS1,1) ],:)=0;
%    allangs140B= [allangs140B; allV140(allS1==1)];
%    allangs40B= [allangs40B; allV40(allS1==1)];
%     allangs8B= [allangs8B; allV8(allS1==1)];
        Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
        
      NT=sum(sum(allS1==1));
         T=   Z.^2;
%       Z= (((1/(NT-1))*((T.*NT-1))));
     if NT>=slim;  allPLVs(:,2,ind)= Z;else; allPLVs(:,2,ind)= NaN;end
     
            allS1= alls1;
   allS1([ 1:1660 ],:)=0; allS1([2340:size(allS1,1) ],:)=0;
     Z=abs(nanmean(exp(1i.*allV8((allS1==1)))));
     
      NT=sum(sum(allS1==1));
          T=   Z.^2;
 %      Z= (((1/(NT-1))*((T.*NT-1))));
     if NT>=slim;  allPLVs(:,3,ind)= Z;else; allPLVs(:,3,ind)= NaN;end
    
    %%%%%%%%%%%%%%%%%
    
    
 %%%%%%
 

for fg=1:size(wavD,1)
mj=mj+1;

allP(:,1:size(wavD,3),mj)= squeeze(nanmean(wavD(fg,:,:),1));

allV(1:length(lfp.trial{fg}),mj)= lfp.trial{fg}(2,:);


sesid=[sesid; ind];
end
allsamp=[allsamp, nanmean(asA,2)];
    end
allIND2{ind}= mr;
end

allP(allP==0)=NaN;
clear allPS
sesU=unique(sesid)
for ih=1:length(unique(sesid))
    allPS(:,:,ih)=nanmean(allP(:,:,sesid==ih),3);
end


%%%%
%%%%%%%%%%%%%%%%%%%%%%

%save('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\optoDBS_40_corr2','allPS', 'freq2','allV2','allPLVs','allsamp')

%%%%%%%%%%%%%%%%%%
end

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2\'
pheight=150;

cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\')
load('optoDBS_40_corr2.mat')

allPS2=allPS;
freq22=freq2;
allV22=allV2;
allPLVs2=allPLVs;
samp2=nanmean(allsamp);
load('optoDBS_140_corr2.mat')
samp=nanmean(allsamp);
%%%%%
%% POW bar plot and quant
ftim=freq2.time< 0.9 %&freq2.time>0.2;
fsel=freq2.freq>=7 & freq2.freq<=9;
F1=squeeze(nanmean(nanmean(allPS(fsel,ftim,:),1),2));
ftim=freq2.time>1.1 &freq2.time<2;
F2=squeeze(nanmean(nanmean(allPS(fsel,ftim,:),1),2));
%%%
ftim=freq22.time< 1;
fsel=freq22.freq>=7 & freq22.freq<=9;
F12=squeeze(nanmean(nanmean(allPS2(fsel,ftim,:),1),2));
ftim=freq22.time>1.1 &freq22.time<2 ;
F22=squeeze(nanmean(nanmean(allPS2(fsel,ftim,:),1),2));

%%%%%
%% Membrane depolarization  %%%
FS=828; %sampling rate
pre_tim=FS./2:FS;
stim_time= FS:FS*2;
X2=mean(allV2(stim_time,:),1);
X1=mean(allV2(pre_tim,:),1);
X22=mean(allV22(stim_time,:),1);
X12=mean(allV22(pre_tim,:),1);
%%%%


%%%%%
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot([X2-X1 ]./samp   , [1-F2./F1]  ,'.','COlor', [ 0.7 0.3 0.3],'Markersize',15); hold on,
plot([X22-X12 ]./samp2    , [1-F22./F12]  ,'.','COlor', [ 0.3 0.3 0.7],'Markersize',15)
M1=[(X2-X1)./samp (X22-X12)./samp2]';M=[1-F2./F1 ;1-F22./F12];
fitResults1 = polyfit(M1,M,1);
yplot1 = polyval(fitResults1,M1);
plot(M1,yplot1,'k-')
axis tight
 [B,BINT,R,RINT,stats] = regress(zscore(M),[ones(length(M),1),M1(~isnan(M))]);
 stats
length(find(~isnan(M)))


Z1=squeeze(allPLVs(:,2,:)-allPLVs(:,1,:));
 Z2=squeeze(allPLVs2(:,2,:)-allPLVs2(:,1,:));
   figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot([X2-X1 ]./samp   , [Z1]  ,'.','COlor', [ 0.7 0.3 0.3],'Markersize',15); hold on,
plot([X22-X12 ]./samp2   , [Z2]  ,'.','COlor', [ 0.3 0.3 0.7],'Markersize',15)
M1=[(X2-X1)./samp (X22-X12)./samp2]';;M=[Z1 ;Z2];
fitResults1 = polyfit(M1(~isnan(M)),M(~isnan(M)),1);
yplot1 = polyval(fitResults1,M1(~isnan(M)));
plot(M1(~isnan(M)),yplot1,'k-')
axis tight
 [B,BINT,R,RINT,stats] = regress(zscore(M(~isnan(M))),[ones(length(M(~isnan(M))),1),M1(~isnan(M))]);
 stats
length(find(~isnan(M)))
 
 
 
