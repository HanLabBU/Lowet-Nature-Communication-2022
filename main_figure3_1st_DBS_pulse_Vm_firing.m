clear all


allsamp=[];
 Z1=[]; Z2=[];Z3=[];
 BRATE=[];SRATE=[];PRATE=[];sesid=[];
   allangs8=[];  allangs140=[];  allangs40=[];
    allangs8B=[];  allangs140B=[];  allangs40B=[];
 %%%%
stim_freq=[0]
M=[];MS=[];FS=828;O=[];VS=[];
if stim_freq==1
    cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\140\')
else
      cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\40\')  
end

% 140 ses= 2  very strong, 1
ses=dir('*.mat')
mj=0;
% 3,4 , no LFP
% 3 good
for ind=[1:length(ses)]
    sesid=ind;
    try
    Cpath=ses(ind).name;
    load(Cpath)
     try
         lfp2=[];
    lfp2=result.lfp;end
% 
% if isempty(lfp2)
%     Cpath=ses(3).name;
%     load(Cpath)
%      try
%          lfp2=[];
%     lfp2=result.lfp;end
% end

    Cpath=ses(ind).name;
    load(Cpath)
    
    
   
    %%%
     FS=828;
    Fn = FS/2;FB=[ 86.4 88.7];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    LFPg= ((filtfilt(B,A,    result.traces(:,1))));
    result.traces(:,1)= result.traces(:,1)-LFPg;
    %%%
   ind
    lfp=[];  asA=[];clear aV allV8 allV40 allV140 allS alls1
     for  ne=unique(result.trial_vec)
         v= result.traces(result.trial_vec==ne,1)./result.tracesB(result.trial_vec==ne)   ;
         v=v(1:2500);
         Fn = FS/2;FB=[ 86.5 88.5];
         
         
         vsub= result.resultS{ne}.trace_ws  ;
         
         [ fitbaseline, coeff]=exp_fit_Fx(vsub',FS);
         
         vsub=(vsub-fitbaseline);

        vsub=zscore(vsub(1:2500))';
  
        
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
          alls1( :,ne)=fastsmooth(vect(1:2500),2,1,1);
          
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


%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2\'
pheight=150;



%%%%5%%%
   MS1=[];   VS1=[];
for tr=1:size( alls1,2)
    try
opto=lfp2.trial{tr}(2,1:end);
windS=500;
trigs=find(diff(opto)>0.2);
for ind=1%:length(trigs)-1
    if trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>6& trigs(ind)<820
     M=[M,  alls1( trigs(ind)-windS:trigs(ind)+windS,tr)];
   %   O=[O,  opto( trigs(ind)-windS:trigs(ind)+windS)'];
    end
    
     if trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5 %& trigs(ind)>818  &trigs(ind)<1650
     MS=[MS,  alls1( trigs(ind)-windS:trigs(ind)+windS,tr)];
        VS=[VS,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
             MS1=[MS1,  alls1( trigs(ind)-windS:trigs(ind)+windS,tr)];
        VS1=[VS1,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
        O=[O,  opto( trigs(ind)-windS:trigs(ind)+windS)'];
        
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
%savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\revision\'

ff= MS(490,:); % reference to 10time points before DBS onset
clear allT
for id=1:size(VS,1)
    [h,p,ci,stats] =  ttest(MS(id,:),ff);
    if p<0.05
        allT(id,:)=stats.tstat;end
end
pheight=160
timax=([-windS:windS]-1).*1.2;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(timax,(nanmean(MS,2)).*FS,'k');hold on,
fill_error_area2(timax, nanmean(MS,2).*FS,   (nanstd(MS,[],2).*FS)./sqrt(size(MS,2)), [ 0.5 0.5 .5])
line([ 0 0],[ 0 100],'Color',[ 0.8 0.3 0.3],'Linestyle','-','Linewidth',1)
for dg=1:250
    if stim_freq==1
        line([ 7.2 7.2].*dg,[ 0 100],'Color',[ 0.8 0.3 0.3],'Linestyle','-','Linewidth',1);else
        line([ 25 25].*dg,[ 0 100],'Color',[ 0.8 0.3 0.3],'Linestyle','-','Linewidth',1);end
end
axis tight
xlim([-100 500])
plot(timax(find(allT>0)),ones(1,length( find(allT>0) ))+180,'.y','Markersize',15)
xlim([-20 40])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff= VS(490,:);
clear allT
for id=1:size(VS,1)
    [h,p,ci,stats] =  ttest(VS(id,:),ff);
    if p<0.05
        allT(id,:)=stats.tstat;end
end
timax=([-windS:windS]-1).*1.2;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(timax,(nanmean(zscore(VS,[],1),2)),'k');hold on,
fill_error_area2(timax, nanmean(zscore(VS,[],1),2),   (nanstd(zscore(VS,[],1),[],2))./sqrt(size(VS,2)), [ 0.5 0.5 .5])
minM=-0.5; maxM=0.6;
line([ 0 0],[ minM maxM],'Color',[ 0.8 0.3 0.3],'Linestyle','-','Linewidth',1)
for dg=1:250
    if stim_freq==1
        line([ 7.2 7.2].*dg,[minM maxM],'Color',[ 0.8 0.3 0.3],'Linestyle','-','Linewidth',1);else
        line([ 25 25].*dg,[minM maxM],'Color',[ 0.8 0.3 0.3],'Linestyle','-','Linewidth',1);end
end
axis tight
xlim([-100 500])
plot(timax(find(allT>0)),ones(1,length( find(allT>0) )),'.y','Markersize',15)
xlim([-20 40])

