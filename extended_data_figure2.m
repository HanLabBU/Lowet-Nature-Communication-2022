addpath('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_Scripts\')


%%
clear all

stim_type=[];
% 1= 140Hz, other 40Hz
allsp=[]; allVm=[]; firing_conc=[];Vm_conc=[];
allsamp=[];stim_type_tr=[];
allangs40=[];allangs140=[];stim_type_sp=[];
mr=0;
miceNames= {'C00014169', 'C00014168','619706', '617432','603841' ,'603831', '615783','603809', '615784'}

for stim_freq=[0 1 3 4]
    
    if stim_freq==1
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\140\')
    elseif  stim_freq==0
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\40\')
    elseif  stim_freq==2
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\40\')
    elseif  stim_freq==3
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\140\')
        
    end
    
    
    ses=dir('*.mat');
    
    for ind=1:length(ses)
        Cpath= ses(ind).name;
        
        
        load(Cpath)  %loading
        
        trial_numb=unique(result.trial_vec);
        if length(find(result.trial_vec==trial_numb(1) )) <4000 & length(trial_numb)>2
            
            
            %%  denoise
            FS=828;
          %  Fn = FS/2;FB=[ 86.4 88.7];
          %  [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
         %   LFPg= ((filtfilt(B,A,    result.traces(:,1))));
            %result.traces(:,1)= result.traces(:,1)-LFPg;
            %%
            lfp=[];
            clear allV allS alls1  allV40 allV140
            tr=0;
            for  ne=unique(result.trial_vec) % trials
                tr=tr+1;
                vsub= result.resultS{ne}.trace_ws(1,:)'  ;
                vsub=vsub(10:2500);
                
                [ fitbaseline, coeff]=exp_fit_Fx(vsub,FS);                
                vsub=(vsub-fitbaseline');
                
                lfp.trial{tr}= ( vsub)';
                allV(1:length(vsub),ne)=vsub;
                spx= result.resultS{ne}.spike_idx{1}   ;
                samp= result.resultS{ne}.spike_amplitude{1};
                samp(spx <10)=NaN;
                
                allsamp=[allsamp ,samp];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            clear allUL
            for iter=1:length(Cpath)
                F2= Cpath(iter);
                F1='_';
                allUL(iter)= strcmp(F2,F1);
            end
            
            clear allp
            for iter=1:length(Cpath)
                F2= Cpath(iter);
                F1='.';
                allp(iter)= strcmp(F2,F1);
            end
            
            %%%%%%%%%%%%%%%%%
            UL=find(allUL); p=find(allp);
            if length(UL)>3
                
                
                nameD= Cpath(UL(4)+1: UL(4)+2) ;
                if ~isempty(str2num(nameD))
                    mr=mr+1
                    data(1,mr)= str2num(nameD);
                    
                    nameD= Cpath(1:UL(1)-1) ;
                    data(2,mr)= find(strcmp(nameD,miceNames));
                    
                    samp=      nanmean(allsamp);
                    sm=50 ;%smoothing parameter (rectangular window)
                    allVmSM(:,mr)= nanfastsmooth(nanmean(allV,2),sm,1,1)./nanmean(samp);
                    allVmSM(:,mr)= allVmSM(:,mr)-nanmean(allVmSM(50:750,mr));
                    
                end
            end
        end
    end % sessions
    
end % stimulation type


%%%%%%%%%%%%%%
mean(data(1,:))
std(data(1,:))


allIntra =[];
allInter=[];clear allN
for gh=unique(data(2,:))
    
    varA= data(1,data(2,:)==gh);
    allInter(gh)=nanmean(varA);
    allIntra(gh)=std(varA);
    allN(gh)= length(varA);
end

nanmean(allIntra(allN>1))
std(allInter)


V1=nanmean(allVmSM(828+350:1656,:),1);
V2=nanmean(allVmSM(300:800,:),1);

depolarization=V1-V2;
[r, p]=corrcoef(data(1,:),depolarization)


figure('COlor','w','Position', [ 300 400 300 180],'Renderer', 'painters')
fitResults1 = polyfit(data(1,:),depolarization,1);
% Evaluate polynomial
yplot1 = polyval(fitResults1,5:5:65);
plot(5:5:65,yplot1,'k','Linewidth',1)
hold on,
f=scatter(data(1,:),depolarization,14,[0.3 0.3 0.8],'filled');
alpha(f,0.5)
xlim([0 70])
fitResults1 = polyfit(data(1,:),depolarization,1);




