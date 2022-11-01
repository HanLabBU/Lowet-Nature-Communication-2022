clear all


allsamp=[];
 Z1=[]; Z2=[];Z3=[];
 BRATE=[];SRATE=[];PRATE=[];sesid=[];
   allangs8=[];  allangs140=[];  allangs40=[];
    allangs8B=[];  allangs140B=[];  allangs40B=[];
 %%%%
stim_freq=[ 1]

if stim_freq==1
cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\140\')
else
  cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\40\')
end


ses=dir('*.mat')
mj=0;
for ind=[1:length(ses)]
    
    load(ses(ind).name)
     trial_numb=unique(result.trial_vec);
    if length(find(result.trial_vec==trial_numb(1) )) <4000 & length(trial_numb)>2
    
    %%%
     FS=828;
    Fn = FS/2;FB=[ 86.5 88.5];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    LFPg= ((filtfilt(B,A,    result.traces(:,1))));
    result.traces(:,1)= result.traces(:,1)-LFPg;
    %%%
   ind
    lfp=[];  asA=[];clear aV allV8 allV40 allV140 allS alls1
     for  ne=unique(result.trial_vec)
         
  v= result.traces(result.trial_vec==ne,1);%



       % exponential fitting to remove photobleaching
     [ fitbaseline, coeff]=exp_fit_Fx(v,FS);
  
   v=(v-fitbaseline');
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
    
     allV8(1:length(v),ne)=Vsub8;
        allV40(1:length(v),ne)=Vsub40;
         allV140(1:length(v),ne)=Vsub140;
       
       %%spikes
          strain= result.resultS{ne}.roaster(1:2500);
          allS(1:length(strain),ne)=strain;
           sid=result.resultS{ne}.spike_idx{1};sid=sid-5;sid(sid<1)=[];
           vect=zeros(1,size(allS,1));vect(sid)=1;
          alls1( :,ne)=vect(1:2500);
          
          %%%%
          lfp.trial{ne}(1,:)= vsub';
aV(:,ne)=v; % raw Vm
  lfp.trial{ne}(2,:)= v'; % raw Vm
  lfp.time{ne}= (1:size(v,1))./828;
    
  spx= result.resultS{ne}.spike_idx{1}   ;
         samp= result.resultS{ne}.spike_amplitude{1}./mean(result.tracesB(result.trial_vec==ne));
         samp(spx <10)=NaN;
         asA=[asA, samp];
     end
  for id=1:2 %:size(vsig,2)
lfp.label{id}= num2str(id);
  end  
 
allV2(:,ind)= nanmean(aV,2);%./nanmean(asA);
 allV2(:,ind)=allV2(:,ind);%-nanmean(allV2(500:800,ind)); 
 
            warning off
  cfg = []; %block_type == cfg.blk
    cfg.method ='wavelet'; %'mvar';
    cfg.output ='fourier';
     cfg.taper='hanning';
    cfg.keeptapers ='yes';
    cfg.keeptrials ='yes';
    cfg.trials='all';cfg.tapsmofrq =5;%
     cfg.channel= 'all'%; %chans=cfg.channel;
    cfg.foi= [4:1:80];
     cfg.toi=lfp.time{1}(1:1:end) ;
     cfg.width =5;
    cfg.t_ftimwin =[ones(1,length(cfg.foi))*0.4];
freq2 = ft_freqanalysis(cfg, lfp);

 wavD = abs(squeeze(freq2.fourierspctrm(:,1,:,:)));
wavA = angle(squeeze(freq2.fourierspctrm(:,1,:,:)));
fsel=freq2.freq>=7 & freq2.freq<=9;

 slim=5;
   %% PLV with filtered signals
   allS1= alls1;
   allS1([1:FS (2*FS):size(allS1,1)],:)=0;
   allangs140= [allangs140; allV140(allS1==1)];
     allangs40= [allangs40; allV40(allS1==1)];
       allangs8= [allangs8; allV8(allS1==1)];
      Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
      NT=sum(sum(allS1==1));
         T=   Z.^2;
%     Z= (((1/(NT-1))*((T.*NT-1))));
     if NT>=slim;  allPLVs(:,1,ind)= Z;else; allPLVs(:,1,ind)= NaN;end
       
        allS1= alls1;
   allS1([ 1:50 FS:size(allS1,1) ],:)=0;
   allangs140B= [allangs140B; allV140(allS1==1)];
   allangs40B= [allangs40B; allV40(allS1==1)];
    allangs8B= [allangs8B; allV8(allS1==1)];
        Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
        
      NT=sum(sum(allS1==1));
         T=   Z.^2;
%       Z= (((1/(NT-1))*((T.*NT-1))));
     if NT>=slim;  allPLVs(:,2,ind)= Z;else; allPLVs(:,2,ind)= NaN;end
     
            allS1= alls1;
   allS1([ 1:(2*FS) ],:)=0; allS1([size(allS1,1)-50 :size(allS1,1) ],:)=0;
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

spx= result.resultS{fg}.spike_idx{1}   ;
         samp= result.resultS{fg}.spike_amplitude{1};
         samp(spx <10)=NaN;
     
allsamp=[allsamp ,samp];
spt=result.resultS{fg}.spike_idx{1}-7;
spt(spt<0)=[];
ssel=spt(spt>3 &spt<828);
if length(ssel)>1
Z1=[Z1,squeeze(wavA(fg,:,spt(spt>5 &spt<828)))];
else;Z1=[Z1,squeeze(wavA(fg,:,spt(spt>5 &spt<828)))'];end
BRATE=[BRATE,(length(spt(spt>3 &spt<825))./822).*1000];
ssel=spt(spt>828 &spt<1628);
if length(ssel)>1
Z2=[Z2,squeeze(wavA(fg,:,ssel))];
else;Z2=[Z2,squeeze(wavA(fg,:,ssel))'];end
SRATE=[SRATE,(length(ssel)./800).*1000];
sesid=[sesid; ind];
ssel=spt(spt>1650 & spt< 2400);
if length(ssel)>1
Z3=[Z3,squeeze(wavA(fg,:,ssel))];
else;Z3=[Z3,squeeze(wavA(fg,:,ssel))'];end
PRATE=[PRATE,(length(ssel)./750).*1000];

end

    end
end


%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2\'
pheight=150;
%%%%%%%%%%%%%%%%%%%%%%

allP(allP==0)=NaN;
clear allPS
sesU=unique(sesid)
for ih=1:length(unique(sesid))
    allPS(:,:,ih)=nanmean(allP(:,:,sesid==ih),3);
end
clear BRATES
for ih=1:length(unique(sesid))
    BRATES(ih)=nanmean( BRATE(sesid==ih));
end
clear SRATES
for ih=1:length(unique(sesid))
    SRATES(ih)=nanmean( SRATE(sesid==ih));
end
clear PRATES
for ih=1:length(unique(sesid))
   PRATES(ih)=nanmean( PRATE(sesid==ih));
end


%%%%%
%% POW bar plot and quant
ftim=freq2.time< 1;
fsel=freq2.freq>=7 & freq2.freq<=9;
F1=squeeze(nanmean(nanmean(allPS(fsel,ftim,:),1),2));
ftim=freq2.time>1.1 &freq2.time<2 ;
F2=squeeze(nanmean(nanmean(allPS(fsel,ftim,:),1),2));
ftim=freq2.time>2.1  ;
F3=squeeze(nanmean(nanmean(allPS(fsel,ftim,:),1),2));

  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
V1=F1;V2=F2;V3=F3;
  M=[V1,V2, V3];
b1=bar(1,nanmean(V1),'Facecolor',[ 0.3 0.3 0.3]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(2,nanmean(V2),'Facecolor',[ 0.7 0.6 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(V3),'Facecolor',[ 0.3 0.3 0.3]);hold on,
set(b1,'FaceAlpha',0.7)
errorbar([1 2 3],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k');ylim([1.5 5])




%%%%%%%%

 
V1=squeeze(allPLVs(:,2,:));
V2=squeeze(allPLVs(:,1,:));
V3=squeeze(allPLVs(:,3,:));
%%%%%%
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
  M=[V1,V2, V3];
b1=bar(1,nanmean(V1),'Facecolor',[ 0.3 0.3 0.3]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(2,nanmean(V2),'Facecolor',[ 0.7 0.6 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(V3),'Facecolor',[ 0.3 0.3 0.3]);hold on,
set(b1,'FaceAlpha',0.7)
errorbar([1 2 3],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
ylim([ 0 0.8])
%%%%
%% FIring RATE

