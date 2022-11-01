clear all


allsamp=[];
 Z1=[]; Z2=[];Z3=[];
 BRATE=[];SRATE=[];PRATE=[];sesid=[];
   allangs8=[];  allangs140=[];  allangs40=[];
    allangs8B=[];  allangs140B=[];  allangs40B=[];
    M=[];MS=[];FS=828;O=[];VS=[];VS1=[];
 %%%%
stim_freq=0; % Select frequency of DBS , 0=40Hz, 1=140Hz
if stim_freq==1
  cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\140\')
else
  cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\40\')
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
% 
if isempty(lfp2)
    Cpath=ses(5).name;
    load(Cpath)
     try
         lfp2=[];
    lfp2=result.lfp;end
end

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
  
  % exponential fitting to remove photobleaching
  [ fitbaseline, coeff]=exp_fit_Fx(v,FS);
  
  v=(v-fitbaseline');
  
  vsub= result.resultS{ne}.trace_ws  ;
  vsub=vsub(1:2500)';
  vsub=zscore(vsub-mean(vsub));
  allV(1:length(vsub),ne)=vsub;
  %%spikes
  sid=result.resultS{ne}.spike_idx{1};sid=sid-0;sid(sid<1)=[];
  vect=zeros(1,size(allV,1));vect(sid)=1;
  if stim_freq==0
      alls1( :,ne)=fastsmooth(vect(1:2500),5,1,1);
  else
      alls1( :,ne)=fastsmooth(vect(1:2500),3,1,1);
  end
  %%%%
  lfp.trial{ne}(1,:)= vsub';
  aV(:,ne)=v;
  lfp.trial{ne}(2,:)= v';
  lfp.time{ne}= (1:size(v,1))./FS;
  
  spx= result.resultS{ne}.spike_idx{1}   ;
         samp= result.resultS{ne}.spike_amplitude{1}./mean(result.tracesB(result.trial_vec==ne));
         samp(spx <10)=NaN;
         asA=[asA, samp];
     end
  for id=1:2 %:size(vsig,2)
lfp.label{id}= num2str(id);
  end  
 
 

    
    %%%%%%%%%%%%%%%%%
    
    
 %%%%%%



%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2\'
pheight=150;



%%%%5%%%
   MS1=[];  M1=[];
    VSI1=[]; VSI=[];
for tr=1:size( alls1,2)
    try
opto=lfp2.trial{tr}(2,1:end);optoTrace=lfp2.trial{tr}(1,1:end);
if stim_freq==1
windS=15;windSL=15;else
windS=25;windSL=25;end


trigs=find(diff(opto)>0.2);
for ind=1:length(trigs)-1
    if trigs(ind)-windSL>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5 & trigs(ind)>FS &trigs(ind)<(2*FS) & mean(optoTrace(trigs(ind)-15:trigs(ind)+20)) <0.1
     M=[M,  alls1( trigs(ind)-windSL:trigs(ind)+windS,tr)];
         VS1=[VS1,  allV( trigs(ind)-windSL:trigs(ind)+windS,tr)-mean( allV( trigs(ind)-windSL:trigs(ind)+windS,tr))];
            VSI1=[VSI1,  allV( trigs(ind)-windSL:trigs(ind)+windS,tr)-mean( allV( trigs(ind)-windSL:trigs(ind)+windS,tr))];
        M1=[M1,  alls1( trigs(ind)-windSL:trigs(ind)+windS,tr)];

    end
    
    if trigs(ind)-windSL>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5 & trigs(ind)>FS  &trigs(ind)<(2*FS) & mean(optoTrace(trigs(ind)-15:trigs(ind)+20)) >0.1
     MS=[MS,  alls1( trigs(ind)-windSL:trigs(ind)+windS,tr)];
       
     VS=[VS,  allV( trigs(ind)-windSL:trigs(ind)+windS,tr)-mean(allV( trigs(ind)-windSL:trigs(ind)+windS,tr))   ];
          VSI=[VSI,  allV( trigs(ind)-windSL:trigs(ind)+windS,tr)-mean(allV( trigs(ind)-windSL:trigs(ind)+windS,tr))   ];
     
        MS1=[MS1,  alls1( trigs(ind)-windSL:trigs(ind)+windS,tr)];
      % VS1=[VS1,  allV( trigs(ind)-windS:trigs(ind)+windS,tr)];
        O=[O,  opto( trigs(ind)-windS:trigs(ind)+windS)'];
        
    end
end;
    end
end

if ~isempty(MS1)

mj=mj+1
if stim_freq==1
allmod_S(:,mj)=nanmean(MS1,2);allmod_S(:,mj)=allmod_S(:,mj)-nanmean(allmod_S(10:15,mj));
allmod_S1(:,mj)=nanmean(M1,2);allmod_S1(:,mj)=allmod_S1(:,mj)-nanmean(allmod_S1(10:15,mj));
allmod_V(:,mj)=nanmean(VSI,2);allmod_V(:,mj)=allmod_V(:,mj)-nanmean(allmod_V(10:15,mj));
allmod_V1(:,mj)=nanmean(VSI1,2);allmod_V1(:,mj)=allmod_V1(:,mj)-nanmean(allmod_V1(10:15,mj));
else
 allmod_S(:,mj)=nanmean(MS1,2);allmod_S(:,mj)=allmod_S(:,mj)-nanmean(allmod_S(20:25,mj));
allmod_S1(:,mj)=nanmean(M1,2);allmod_S1(:,mj)=allmod_S1(:,mj)-nanmean(allmod_S1(20:25,mj));
allmod_V(:,mj)=nanmean(VSI,2);allmod_V(:,mj)=allmod_V(:,mj)-nanmean(allmod_V(20:25,mj));
allmod_V1(:,mj)=nanmean(VSI1,2);allmod_V1(:,mj)=allmod_V1(:,mj)-nanmean(allmod_V1(20:25,mj));   
end

end


allV2(:,mj)= nanmean(aV,2)./nanmean(asA);
 allV2(:,mj)=allV2(:,mj)-nanmean(allV2(500:800,mj)); 

 
allses(mj)=     sesid;

    end

end



savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\com_puls_opto\'


pheight=160;
timax=([-windSL:windS]-1).*1.2;
if stim_freq== 1
tsel=timax>-3 & timax<7;else;
tsel=timax>-10 & timax<25;end;
DD=minmax(allmod_V(tsel,:)');MM=DD(:,2)-DD(:,1);
DD1=minmax(allmod_V1(tsel,:)');;MM1=DD1(:,2)-DD1(:,1);

[h,p,ci,stats] = ttest(MM,MM1)
%   0.0878
% df=20, t=-4.8
  figure('COlor','w','Position', [ 300 400 160 pheight],'Renderer', 'painters') 
  M=[MM1,MM];
  b1=bar(1,mean(MM1),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(2,mean(MM),'Facecolor',[ 0 0. 0.9]);hold on,
set(b1,'FaceAlpha',0.7)
errorbar([1 2  ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
if stim_freq== 1
ylim([0 .1]);else
ylim([0.2 0.8]);end
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_VM_14OBAR.pdf'])


if stim_freq== 1
tsel=timax>-3 & timax<7;else;
tsel=timax>-10 & timax<25;end;
DD=minmax(allmod_S(tsel,:)');MM=(DD(:,2)-DD(:,1)).*FS;
DD1=minmax(allmod_S1(tsel,:)');;MM1=(DD1(:,2)-DD1(:,1)).*FS;
[h,p,ci,stats] = ttest(MM,MM1)


  figure('COlor','w','Position', [ 300 400 160 pheight],'Renderer', 'painters') 
  M=[MM1,MM];
  b1=bar(1,mean(MM1),'Facecolor',[ 0.5 0.5 0.5]);hold on,
set(b1,'FaceAlpha',0.7)
  b1=bar(2,mean(MM),'Facecolor',[ 0 0. 0.9]);hold on,
set(b1,'FaceAlpha',0.7)
errorbar([1 2  ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
if stim_freq== 1
ylim([0 5]); else;
ylim([0 40]); end