clear all


allsamp=[];
 Z1=[]; Z2=[];Z3=[];
 BRATE=[];SRATE=[];PRATE=[];sesid=[];
   allangs8=[];  allangs140=[];  allangs40=[];
    allangs8B=[];  allangs140B=[];  allangs40B=[];
 %%%%
stim_freq=[0]
M=[];MS=[];FS=828;O=[];VS=[];V=[];VS2=[];MS2=[];
if stim_freq==0
cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DATA\optoDBS\140\')
else
  cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DATA\optoDBS\40\')
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
  v= result.traces(result.trial_vec==ne,1)./result.tracesB(result.trial_vec==ne)   ;
v=v(1:2500);
 
       % exponential fitting to remove photobleaching
     [ fitbaseline, coeff]=exp_fit_Fx(v,FS);
  
   v=(v-fitbaseline');
  
  
      vsub= result.resultS{ne}.trace_ws  ;
        vsub=vsub(1:2500)';
           vsub=(vsub-fitbaseline');
        vsub=zscore(vsub-mean(vsub));
        
               Fn = FS/2;FB=[ 7  9];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    Vsub8= angle(hilbert(filtfilt(B,A,     vsub)));
            Fn = FS/2;FB=[ 35 45];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    Vsub40= angle(hilbert(filtfilt(B,A,     vsub)));
       Fn = FS/2;FB=[ 135 145];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    Vsub140= angle(hilbert(filtfilt(B,A,     vsub)));
    
     allV8(1:length(Vsub8),ne)=Vsub8;
     allV40(1:length(Vsub40),ne)=Vsub40;
       allV140(1:length(Vsub140),ne)=Vsub140;
          allV(1:length(vsub),ne)=vsub;
       %%spikes
          strain= result.resultS{ne}.roaster(1:2500);
          allS(1:length(strain),ne)=strain;
           sid=result.resultS{ne}.spike_idx{1};sid=sid-0;sid(sid<1)=[];
           vect=zeros(1,size(allS,1));vect(sid)=1;
          alls1( :,ne)=vect(1:2500);
          
          %%%%
          lfp.trial{ne}(1,:)= vsub';
aV(:,ne)=v;
  lfp.trial{ne}(2,:)= v';
  lfp.time{ne}= (1:size(v,1))./828;
    
  spx= result.resultS{ne}.spike_idx{1}   ;
         samp= result.resultS{ne}.spike_amplitude{1}./mean(result.tracesB(result.trial_vec==ne));
         samp(spx <10)=NaN;
         asA=[asA, samp];
     end
  for id=1:2 %:size(vsig,2)
lfp.label{id}= num2str(id);
  end  
 
 
 
   %% PLV with filtered signals
   allS1= alls1;
   allS1([1:825 1650:size(allS1,1)],:)=0;
   allangs140= [allangs140; allV140(allS1==1)];
     allangs40= [allangs40; allV40(allS1==1)];
       allangs8= [allangs8; allV8(allS1==1)];
      Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
      NT=sum(sum(allS1==1));
     if NT>7;  allPLVs(:,1,ind)= Z;else; allPLVs(:,1,ind)= NaN;end
       
        allS1= alls1;
   allS1([1:20 830:size(allS1,1) ],:)=0;
   allangs140B= [allangs140B; allV140(allS1==1)];
   allangs40B= [allangs40B; allV40(allS1==1)];
    allangs8B= [allangs8B; allV8(allS1==1)];
        Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
      NT=sum(sum(allS1==1));
     if NT>7;  allPLVs(:,2,ind)= Z;else; allPLVs(:,2,ind)= NaN;end
     
            allS1= alls1;
   allS1([ 1:1650 ],:)=0;
     Z=abs(nanmean(exp(1i.*allV8(allS1==1))));
      NT=sum(sum(allS1==1));
     if NT>7;  allPLVs(:,3,ind)= Z;else; allPLVs(:,3,ind)= NaN;end
    
    %%%%%%%%%%%%%%%%%
    
    
 %%%%%%



%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2\'
pheight=150;



%%%%5%%%
   MS1=[];   VS1=[];
for tr=1:size( alls1,2)
    try
opto=lfp2.trial{tr}(1,1:end);
windS=230;
trigs=find(diff(opto)>0.2);
for ind=1:length(trigs)-1
    if trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5& trigs(ind)<820& trigs(ind)>50
     M=[M,  fastsmooth(alls1( trigs(ind)-windS:trigs(ind)+windS,tr),5,1,1)];
    V=[V,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
    end
    
     if trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5 & trigs(ind)>900  &trigs(ind)<1650
     MS=[MS,  fastsmooth(alls1( trigs(ind)-windS:trigs(ind)+windS,tr),5,1,1)];
        VS=[VS,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
             MS1=[MS1,  alls1( trigs(ind)-windS:trigs(ind)+windS,tr)];
        VS1=[VS1,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
        O=[O,  opto( trigs(ind)-windS:trigs(ind)+windS)'];
     elseif trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5 & trigs(ind)>1700
     MS2=[MS2,  fastsmooth(alls1( trigs(ind)-windS:trigs(ind)+windS,tr),5,1,1)];
        VS2=[VS2,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
   
    end
end;
    end
end

if ~isempty(MS1)

mj=mj+1;

allmod_S(:,mj)=nanmean(MS1,2);
allmod_V(:,mj)=nanmean(VS1,2);allmod_V(:,mj)=allmod_V(:,mj)-nanmean(allmod_V(:,mj));

end


allV2(:,mj)= nanmean(aV,2)./nanmean(asA);
 allV2(:,mj)=allV2(:,mj)-nanmean(allV2(500:800,mj)); 

 
allses(mj)=     sesid;

    end

end



savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\revision\'
pheight=160
timax=([-windS:windS]-1).*1.2;
selT=find(timax>-50 & timax<0);
msig=fastsmooth((nanmean(MS,2)),10,1,1).*FS;
FS=828;Fn = FS/2;FB=[ 36 44];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]); filt_T=((filtfilt(B,A, msig))); % zero-phasefilter
msig=msig-filt_T;
FB=[ 76 84];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]); filt_T=((filtfilt(B,A, msig))); % zero-phasefilter
msig=msig-filt_T;

pheight=160
timax=([-windS:windS]-1).*1.2;
  figure('COlor','w','Position', [ 300 200 250 160],'Renderer', 'painters')
plot(timax,msig-nanmean(msig(selT)),'r');hold on,
fill_error_area2(timax, msig-nanmean(msig(selT)) ,   (nanstd(MS,[],2).*FS)./sqrt(size(MS,2)), [ 0.5 0.5 .5])
plot(timax,fastsmooth(nanmean(M,2),10,1,1).*FS - mean(fastsmooth(nanmean(M(selT,:),2),1,1,1).*FS) ,'k');hold on,
fill_error_area2(timax, fastsmooth(nanmean(M,2).*FS,10,1,1)- mean(fastsmooth(nanmean(M(selT,:),2),1,1,1).*FS),   (nanstd(M,[],2).*FS)./sqrt(size(M,2)), [ 0.5 0.5 .5])
%plot(timax,fastsmooth((nanmean(MS2,2)).*FS,10,1,1)- mean(fastsmooth(nanmean(MS2(selT,:),2),1,1,1).*FS),'m');hold on,
%fill_error_area2(timax, fastsmooth(nanmean(MS2,2).*FS,10,1,1)- mean(fastsmooth(nanmean(MS2(selT,:),2),1,1,1).*FS),   (nanstd(MS2,[],2).*FS)./sqrt(size(MS2,2)), [ 0.5 0.5 .5])
plot(timax,msig-nanmean(msig(selT)),'r');hold on,
plot(timax,fastsmooth(nanmean(M,2),10,1,1).*FS - mean(fastsmooth(nanmean(M(selT,:),2),1,1,1).*FS) ,'k');hold on,
%plot(timax,fastsmooth((nanmean(MS2,2)).*FS,10,1,1)- mean(fastsmooth(nanmean(MS2(selT,:),2),1,1,1).*FS),'m');hold on,
axis tight
xlim([-50 200])
ylim([ -2 10])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OPTOPuls_trig_FR_base_stim_post' num2str(stim_freq) '.pdf'])



pheight=160
timax=([-windS:windS]-1).*1.2;
  figure('COlor','w','Position', [ 300 200 200 400],'Renderer', 'painters')
  nx=0;
for ind=1:size(M,2)
    n=M(:,ind);n(n==0)=NaN;
    if nansum(n)>0
        nx=nx+1;
    plot(timax,n+nx,'.b'); hold on,    end
end;nx1=0;
for ind=1:size(M,2)
    n=MS(:,ind);n(n==0)=NaN; if nansum(n)>0
          nx1=nx1+1;
    plot(timax,n+nx1+nx,'.k'); hold on,  end  
end;nx2=0;
for ind=1:size(M,2)
    n=MS2(:,ind);n(n==0)=NaN;if nansum(n)>0
          nx2=nx2+1;
    plot(timax,n+nx2+nx+nx1,'.m'); hold on, end   
end
% line([ 25 25],[ 0 30],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
axis tight
xlim([-80 100])
%xlim([-10 28])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_FR_pop40hzOPTO.pdf'])
selT=find(timax>-50 & timax<0);
msig=fastsmooth((nanmean(VS,2)),1,1,1);
FS=828;Fn = FS/2;FB=[ 37 43];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]); filt_T=((filtfilt(B,A, msig))); % zero-phasefilter
msig=msig-filt_T;
FB=[ 77 83];
[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]); filt_T=((filtfilt(B,A, msig))); % zero-phasefilter
msig=msig-filt_T;
timax=([-windS:windS]-1).*1.2;
  figure('COlor','w','Position', [ 300 200 250 160],'Renderer', 'painters')
plot(timax,msig-nanmean(msig(selT)),'r');hold on,
fill_error_area2(timax,msig-nanmean(msig(selT)),   (nanstd((VS),[],2))./sqrt(size(VS,2)), [ 0.5 0.5 .5])
plot(timax,(nanmean((V),2))-mean(nanmean((V(selT,:)),2)),'k');hold on,
fill_error_area2(timax, nanmean((V),2)-mean(nanmean((V(selT,:)),2)),   (nanstd((V),[],2))./sqrt(size(V,2)), [ 0.5 0.5 .5])
plot(timax,(nanmean((V),2))-mean(nanmean((V(selT,:)),2)),'k');hold on,
plot(timax,msig-nanmean(msig(selT)),'r');hold on,

axis tight
xlim([-50 200])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'OPTOPuls_trig_VM_base_stim_post' num2str(stim_freq) '.pdf'])

%xlim([-10 28])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_VM_pop40hzOPTO.pdf'])

%xlim([-10 10])
if 0
timax=([-windS:windS]-1).*1.2;
timL=find(timax>-10  & timax<27);
MM=zscore(allmod_S(timL,:));
MM1=zscore(allmod_V(timL,:));
[AA BB]= max(MM1(:,:),[],1);
  [BB1 BB2]=sort(BB);
  figure('COlor','w','Position', [ 300 400 200 300],'Renderer', 'painters')
 imagesc(timax(timL),[],smoothn(MM1(:,BB2),0.4)')
 axis xy;hold on,
 %plot(timax,[],nanmean( O(timL,:) ,2).*4,'k');axis tight
minM=-0; maxM=18;
line([ 0 0],[ minM maxM],'Color',[ 1 1 1],'Linestyle','--')
line([25 25],[ minM maxM],'Color',[ 1 1 1],'Linestyle','--')
colormap(jet)
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_imag_40Hz_VmOPTO.pdf'])
timax=([-windS:windS]-1).*1.2;
timL=find(timax>-10  & timax<27);
MM=zscore(allmod_S(timL,:));
MM1=zscore(allmod_V(timL,:));
[AA BB]= max(MM(:,:),[],1);
  [BB1 BB2]=sort(BB);
  figure('COlor','w','Position', [ 300 400 200 300],'Renderer', 'painters')
 imagesc(timax(timL),[],smoothn(MM(:,BB2),0.4)')
 axis xy;hold on,
 %plot(timax,[],nanmean( O(timL,:) ,2).*4,'k');axis tight
minM=-0; maxM=18;
line([ 0 0],[ minM maxM],'Color',[ 1 1 1],'Linestyle','--')
line([25 25],[ minM maxM],'Color',[ 1 1 1],'Linestyle','--')
colormap(jet)
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_imag_40Hz_FROPTO.pdf'])
 
 
 
  figure('COlor','w','Position', [ 300 400 200 300],'Renderer', 'painters')
  for id=1:size(allmod_V,2)
plot([-windS:windS]-1,zscore(allmod_V(:,id))+id*4,'r');hold on,
plot([-windS:windS]-1,zscore(allmod_S(:,id))+id*4,'k');hold on,
  end
plot([-windS:windS]-1,nanmean( O,2).*10,'b');axis tight
%xlim([-10 10])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_Vm_FR_shifted40HzOPTO.pdf'])


 st=std(allmod_S);col=colormap(jet(size(allmod_S,2)));
timax=([-windS:windS]-1).*1.2;[Aa BB]= sort(st);
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
  for ind=1:size(allmod_S,2)
plot(timax,allmod_S(:,BB(ind)).*FS,'Color',col(ind,:));hold on,end
minM=-0; maxM=240;
line([ 0 0],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
line([ 25 25],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
axis tight
xlim([-10 28])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_FR_colorplot40HzOPTO.pdf'])


 st=std(allmod_V);col=colormap(jet(size(allmod_V,2)));
timax=([-windS:windS]-1).*1.2;[Aa BB]= sort(st);
  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
  for ind=1:size(allmod_S,2)
plot(timax,allmod_V(:,BB(ind)),'Color',col(ind,:));hold on,end
minM=-0.4; maxM=1;
line([ 0 0],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
line([ 25 25],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
xlim([-10 28])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_Vm_colorplot40HzOPTO.pdf'])


%   figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
% for ind=1:size(MS,2)
%  Sspikes=(find(MS(:,ind)')./FS)-0.0;
% line([ Sspikes; Sspikes],[  ones(1,length(Sspikes))+ind; ones(1,length(Sspikes))+ind+2.8],'COlor', [ 0.2 0.2 0.2],'Linewidth',2)
% hold on,
% end;plot((1:size(O,1))./FS- 0.00,nanmean(O,2).*323,'b','Linewidth',1)
% axis tight;%xlim([-0.05 0.05])
% %print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath ['example_neuron_optotrigStim_' num2str(stim_freq) '.pdf']])
% %savefig(gcf, [ savepath ['example_neuron_optotrigStim_' num2str(stim_freq) '.fig']])
% 
%%%%%%%%%
end
