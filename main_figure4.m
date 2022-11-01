clear all

waveP=[]; wavePp=[]; wavePb=[];
allsamp=[];
 Z1=[]; Z2=[];Z3=[];
 BRATE=[];SRATE=[];PRATE=[];sesid=[];
   allangs8=[];  allangs140=[];  allangs40=[];
    allangs8B=[];  allangs140B=[];  allangs40B=[];
    allwS=[];   allwSp=[];   allwSb=[];
    M=[];MS=[];O=[];VS=[]; VS1=[];
 %%%%
stim_freq=1
FS=828;
if stim_freq==0
cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\140\')  
else
   cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\40\')
end

ses=dir('*.mat')
mj=0;
for ind=[1:length(ses)]
    sesid=ind;
    try
    Cpath=ses(ind).name;
    load(Cpath)
     try
         lfp2=[];
    lfp2=result.lfp;end

    Cpath=ses(ind).name;
    load(Cpath)

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
         v=v(1:2500);
         
         vsub= result.resultS{ne}.trace_ws  ;
         vsub=vsub(1:2500)';vsub(1:10)=vsub(11);vsub(end-10:end)=vsub(end-11);
         
         [ fitbaseline, coeff]=exp_fit_Fx(vsub,FS);
         
         vsub=(vsub-fitbaseline');   
         vsub=zscore(vsub-mean(vsub));
         
         filt_freqs=[3:2:200]
         for citer=1:length(filt_freqs)
             Fn = FS/2;FB=[  filt_freqs(citer)-2  filt_freqs(citer)+2];
             [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
             Vsub8= angle(hilbert(filtfilt(B,A,     vsub))); % Filtering and Hilbert
             allV8(1:length(Vsub8),ne,citer)=Vsub8;
             
         end
         
         allV(1:length(vsub),ne)=vsub;
         %%spikes
         strain= result.resultS{ne}.roaster(1:2500);
         allS(1:length(strain),ne)=strain;
         sid=result.resultS{ne}.spike_idx{1};sid=sid-0;sid(sid<1)=[];
         vect=zeros(1,size(allS,1));vect(sid)=1;
         alls1( :,ne)=vect(1:2500);
         
         %%%%
         
         aV(:,ne)=v;
         
         
         spx= result.resultS{ne}.spike_idx{1}   ;
         samp= result.resultS{ne}.spike_amplitude{1}./mean(result.tracesB(result.trial_vec==ne));
         samp(spx <10)=NaN;
         asA=[asA, samp];
     end

 
   %% PLV with filtered signals  
      allS1= alls1;
   allS1([1:10 FS:size(allS1,1)],:)=0;   
  allMs=[];
   for id=1:size(allV8,3)
   filts= squeeze(allV8(:,:,id));
  allMs=[allMs, filts(allS1==1)];
   end
   
   wavePb=[wavePb;allMs];
   
   bSp=size(allMs,1);
   if bSp>5
       allwSb= [allwSb;abs(nanmean(exp(1i.*allMs),1))];
   end
   
   allS1= alls1;
   allS1([1:FS (2*FS):size(allS1,1)],:)=0;
   
   allMs=[];
   for id=1:size(allV8,3)
       filts= squeeze(allV8(:,:,id));
       allMs=[allMs, filts(allS1==1)];
   end
   
   
   waveP=[waveP;allMs];
   bSp=size(allMs,1);
   if bSp>5
       allwS= [allwS;abs(nanmean(exp(1i.*allMs),1))];
   end

 
   %%%%%%%%%%%%%%%%%
   allS1= alls1;
   allS1([1:(2*FS) size(allS1,1)-20:size(allS1,1)],:)=0;
   
   allMs=[];
   for id=1:size(allV8,3)
       filts= squeeze(allV8(:,:,id));
       allMs=[allMs, filts(allS1==1)];
   end
   
   wavePp=[wavePp;allMs];
   
   bSp=size(allMs,1);
   if bSp>5
       allwSp= [allwSp;abs(nanmean(exp(1i.*allMs),1))];
   end
   

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2\'
pheight=150;

 allV2(:,mj)=allV2(:,mj)-nanmean(allV2(500:800,mj));  
allses(mj)=     sesid;

    end

end



savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\Vmcoh\'
pheight=160

Z= abs(nanmean(exp(1i.*waveP(:,:)),1));
Zb= abs(nanmean(exp(1i.*wavePb(:,:)),1));
Zp= abs(nanmean(exp(1i.*wavePp(:,:)),1));


shuffle_size=150;itermax=1000;
[n1 fn]=size(waveP);
clear Z
for iter=1:itermax
fx=randperm(n1);
Z(:,iter)= abs(nanmean(exp(1i.*waveP(fx(1:shuffle_size),:)),1));
end
[n1 fn]=size(wavePp);
clear Zp
for iter=1:itermax
fx=randperm(n1);
Zp(:,iter)= abs(nanmean(exp(1i.*wavePp(fx(1:shuffle_size),:)),1));
end
[n1 fn]=size(wavePb);
clear Zb
for iter=1:itermax
fx=randperm(n1);
Zb(:,iter)= abs(nanmean(exp(1i.*wavePb(fx(1:shuffle_size),:)),1));
end

  figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
plot(filt_freqs,nanmean(Z,2))
fill_error_area2(filt_freqs,nanmean(Z,2),nanstd(Z,[],2),[ 0.7 0 0 ])
fill_error_area2(filt_freqs,nanmean(Zb,2),nanstd(Zb,[],2),[ 0. 0 0.7 ])
%fill_error_area2(filt_freqs,nanmean(Zp,2),nanstd(Zp,[],2),[ 0.2 0.6 0 ])
plot(filt_freqs,nanmean(Z,2),'Linewidth',1,'Color', [ 0.7 0 0 ])
plot(filt_freqs,nanmean(Zb,2),'Linewidth',1,'Color', [ 0. 0 0.7 ])
%plot(filt_freqs,nanmean(Zp,2),'Linewidth',1,'Color', [ 0.2 0.6 0 ])
axis tight
set(gca,'Xscale','log')
ylim([0 0.8])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vmcoh_base_stim_' num2str(stim_freq) '.pdf'])

fsel=filt_freqs> 3 &filt_freqs< 12
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
  M=[nanmean(Zb(fsel,:),1)',nanmean(Z(fsel,:),1)',nanmean(Zp(fsel,:),1)'];
  b1=bar(1,mean(nanmean(Zb(fsel,:),1)),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(2,mean(nanmean(Z(fsel,:),1)),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(3,mean(nanmean(Zp(fsel,:),1)),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.3)
errorbar([1 2 3 ],nanmean(M,1), nanstd(M),'.k')
%ylim([0 1.2])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vmcoh_base_stimBARTHETA_' num2str(stim_freq) '.pdf'])

%%%%%%%%%%%%%%%%%

fsel=filt_freqs> 138 &filt_freqs< 142
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
  M=[nanmean(Zb(fsel,:),1)',nanmean(Z(fsel,:),1)',nanmean(Zp(fsel,:),1)'];
  b1=bar(1,mean(nanmean(Zb(fsel,:),1)),'Facecolor',[0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(2,mean(nanmean(Z(fsel,:),1)),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(3,mean(nanmean(Zp(fsel,:),1)),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.3)
errorbar([1 2 3 ],nanmean(M,1), nanstd(M),'.k')
ylim([0 0.8])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vmcoh_base_stimBAR140_' num2str(stim_freq) '.pdf'])


fsel=filt_freqs> 38 &filt_freqs< 42
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
  M=[nanmean(Zb(fsel,:),1)',nanmean(Z(fsel,:),1)',nanmean(Zp(fsel,:),1)'];
  b1=bar(1,mean(nanmean(Zb(fsel,:),1)),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(2,mean(nanmean(Z(fsel,:),1)),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(3,mean(nanmean(Zp(fsel,:),1)),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.3)
errorbar([1 2 3 ],nanmean(M,1), nanstd(M),'.k')
%ylim([0 1.2])
ylim([0 0.8])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vmcoh_base_stimBAR40_' num2str(stim_freq) '.pdf'])



%%%%%%%%%%%%%%%%%%%%%


