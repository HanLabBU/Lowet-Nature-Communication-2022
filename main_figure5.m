addpath('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_Scripts\')

%%
clear all

stim_type=[];
% 1= 140Hz, other 40Hz
allsp=[]; allVm=[]; firing_conc=[];Vm_conc=[];
allsamp=[];stim_type_tr=[];
allangs40=[];allangs140=[];stim_type_sp=[];
mr=0;

for stim_freq=[ 0 1]
    
    if stim_freq==1
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\140\')
    else
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_volt_data_github\optoDBS\40\')
    end
    
    
    ses=dir('*.mat');
    
    for ind=1:length(ses)
        Cpath= ses(ind).name;
        load(Cpath)  %loading
        
        trial_numb=unique(result.trial_vec);
        if length(find(result.trial_vec==trial_numb(1) )) <3500 & length(trial_numb)>2
            
            %%  denoise
            FS=828;
            Fn = FS/2;FB=[ 86.5 88.5];
            [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
            LFPg= ((filtfilt(B,A,    result.traces(:,1))));
           % result.traces(:,1)= result.traces(:,1)-LFPg;
            %%
            lfp=[];
            clear allV allS alls1  allV40 allV140
            tr=0;
            for  ne=unique(result.trial_vec) % trials
                tr=tr+1;
                %% Vm
                v= result.traces(result.trial_vec==ne,1)  ;
                v=v(10:2500);
                
                vsub= result.resultS{ne}.trace_ws(1,:)'  ;
                vsub=vsub(10:2500);
                
                
                
                
                % exponential fitting to remove photobleaching
                [ fitbaseline, coeff]=exp_fit_Fx(v,FS);
                
                v=(v-fitbaseline');
                % v(1:10)=NaN;
                lfp.trial{tr}= ( vsub)';
                lfp.time{tr}= (1:size(v,1))./FS;
                allV(1:length(v),ne)=v;
                %%spikes
                strain= result.resultS{ne}.roaster(10:2500);
                allS(1:length(strain),tr)=strain;
                sid=result.resultS{ne}.spike_idx{1};sid=sid-15;sid(sid<1)=[];
                vect=zeros(1,size(allS,1));vect(sid)=1;
                alls1( :,tr)=vect(1:2491);
                spx= result.resultS{ne}.spike_idx{1}   ;
                samp= result.resultS{ne}.spike_amplitude{1};
                samp(spx <10)=NaN;
                
                allsamp=[allsamp ,samp];
                
            end
            for id=1%:size(vsig,2)
                lfp.label{id}= num2str(id);
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            mr=mr+1;
            allname{mr}=  Cpath;
            
            
            allfiring(:,mr)= nanmean(allS.*FS,2);
            sm=10 ;%smoothing parameter (rectangular window)
            allfiringSM(:,mr)= nanfastsmooth(nanmean(allS.*FS,2),sm,1,1);
            %  allfiringSM(:,mr)= allfiringSM(:,mr)-nanmean(allfiringSM(5:800,mr));
            allVm(:,mr)= nanmean(allV,2)./nanmean(     allsamp);
            sm=50 ;%smoothing parameter (rectangular window)
            allVmSM(:,mr)= nanfastsmooth(nanmean(allV,2),sm,1,1)./nanmean(     allsamp);
            allVmSM(:,mr)= allVmSM(:,mr)-nanmean(allVmSM(50:750,mr));
            firing_conc=[ firing_conc,allS];
            Vm_conc=[Vm_conc,bsxfun(@rdivide,allV, nanmean(     allsamp))];
            stim_type=[stim_type, stim_freq];
            stim_type_tr=[ stim_type_tr, ones(1,size(allS,2)).*stim_freq];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            baseTimsel= [100:800];
            StimTimsel= [830:1600];
            PostTimsel= [1690:2490];
            
            ;
            
            
            %%%%%%%%%%%%%%
            %% PLV with filtered signals
            allS1= alls1;
            allS1([1:800 1650:size(allS1,2)],:)=0;
            
       %     stim_type_sp=[ stim_type_sp, ones(1,length(find(allS1==1))).*stim_freq];
            
        end
    end % sessions
    
end % stimulation type
FS=828;
baseTimsel= round([FS.*0.4:FS.*0.9]);
StimTimsel= round([FS:FS*2]);
StimTimselTR= round([FS:FS+FS.*0.2]);
StimTimselSU= round([FS+FS.*0.2:FS*2]);
PostTimsel= round([FS*2:2500-10]);

tim_axis=([1:size(allV,1)]-(FS-10))./FS;

%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\figure2_suppl\'
pheight=150;
%%%%%%%%%%%%%%%%%%
V1=nanmean(allVm(:, stim_type==0),2);
V1s=nanstd(allVm(:, stim_type==0),[],2)./sqrt(length(find(stim_type==0)));
V2=nanmean(allVm(:, stim_type==1),2);
V2s=nanstd(allVm(:, stim_type==1),[],2)./sqrt(length(find(stim_type==1)));
%% NOn smoothed Vm , overlay 40 and 140
figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
plot(tim_axis,V1,'b','Linewidth',1.5)
hold on,plot(tim_axis,V2,'r','Linewidth',1.5)
axis tight;xlim([-0.7 1.7])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Nonsmooth_Vm_av_40_140opto.pdf'])
%savefig(gcf, [ savepath 'Nonsmooth_Vm_av_40_140opto.fig'])

%% Vm average plots (smoothed)
V1=nanmean(allVmSM(:, stim_type==0),2);
V1s=nanstd(allVmSM(:, stim_type==0),[],2)./sqrt(length(find(stim_type==0)));
V2=nanmean(allVmSM(:, stim_type==1),2);
V2s=nanstd(allVmSM(:, stim_type==1),[],2)./sqrt(length(find(stim_type==1)));

figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
plot(tim_axis,V1,'k','Linewidth',1.5)
fill_error_area2(tim_axis,V1,V1s, [ 0.5 0.5 0.5]);
axis tight;; xlim([-0.7 1.7]);ylim([-0.6 1.2])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'smooth_Vm_av_40opto.pdf'])
%savefig(gcf, [ savepath 'smooth_Vm_av_40opto.fig'])

figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
plot(tim_axis,V2,'k','Linewidth',1.5)
fill_error_area2(tim_axis,V2,V2s, [ 0.5 0.5 0.5]);
axis tight; xlim([-0.7 1.7]);ylim([-0.6 1.2])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'smooth_Vm_av_140opto.pdf'])
%savefig(gcf, [ savepath 'smooth_Vm_av_140opto.fig'])

%% Vm bar plots /quantifications

V1b=nanmean(allVmSM(baseTimsel, stim_type==0),1);
V2b=nanmean(allVmSM(baseTimsel, stim_type==1),1);
V1s=nanmean(allVmSM(StimTimselTR, stim_type==0),1);
V2s=nanmean(allVmSM(StimTimselTR, stim_type==1),1);
V1s2=nanmean(allVmSM(StimTimselSU, stim_type==0),1);
V2s2=nanmean(allVmSM(StimTimselSU, stim_type==1),1);
V1p=nanmean(allVmSM(PostTimsel, stim_type==0),1);
V2p=nanmean(allVmSM(PostTimsel, stim_type==1),1);
%% STATS
%% Base vs trans
[h,p,ci,stats] = ttest(V1s) % 40Hz%
%df=19,  0.0059
[h,p,ci,stats] = ttest(V2s)% 140Hz
%df=20,     0.0093
%% Base vs sustained
[h,p,ci,stats] = ttest( V1s2) % 40Hz
%df=19,      0.0040
[h,p,ci,stats] = ttest( V2s2)% 140Hz
%df=20,   0.0024
%%


% %%
% %% Base vs post
% [h,p,ci,stats] = ttest( V1p) % 40Hz
% % df=20, 0.3230
% [h,p,ci,stats] = ttest( V2p)% 140Hz
% %df=22, 0.5219
% %% trans between DBS cond
% [h,p,ci,stats] = ttest2(V2s, V1s) %
% %df=42,  0.0082
% %% sust between DBS cond
% [h,p,ci,stats] = ttest2(V2s2, V1s2) %
% %df=42,  0.3584

figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
M=[V1s', V1s2'];M2=[ V2s', V2s2'];
b1=bar(2,nanmean(V1s),'Facecolor',[ 0 0 0.9]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(V1s2),'Facecolor',[ 0 0 0.9]);hold on,
set(b1,'FaceAlpha',0.4)
b1=bar(5,nanmean(V2s),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(6,nanmean(V2s2),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
errorbar([2 3 ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
errorbar([5 6],nanmean(M2,1), nanstd(M2)./sqrt(size(M2,1)),'.k')
ylim([0 1.2])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_av_barplotopto.pdf'])
%savefig(gcf, [ savepath 'Vm_av_barplotopto.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRING RATE %%


V1=nanmean(allfiringSM(:, stim_type==0),2);
V1s=nanstd(allfiringSM(:, stim_type==0),[],2)./sqrt(length(find(stim_type==0)));
V2=nanmean(allfiringSM(:, stim_type==1),2);
V2s=nanstd(allfiringSM(:, stim_type==1),[],2)./sqrt(length(find(stim_type==1)));
SPM=firing_conc(:,stim_type_tr==0)==1;
figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
fill_error_area2(tim_axis,fastsmooth(V1,10,1,1),V1s, [ 0.5 0.5 0.5]);
plot(tim_axis,fastsmooth(V1,10,1,1),'k','Linewidth',1)
plot(tim_axis,fastsmooth(V1,300,1,1),'k','Linewidth',1,'Color',[0.7 0.5 0.5])
axis tight;; xlim([-0.7 1.7]);%ylim([-0.3 1.2])
% for ind=1:size(SPM,2)
%    plot(tim_axis(SPM(:,ind)),  ones(1,length(find(SPM(:,ind))))*ind+32,'r.');  hold on,
% end
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_40_smopto.pdf'])
%savefig(gcf, [ savepath 'Firing_40_smopto.fig'])

figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
fill_error_area2(tim_axis,fastsmooth(V2,10,1,1),V2s, [ 0.5 0.5 0.5]);
plot(tim_axis,fastsmooth(V2,10,1,1),'k','Linewidth',1.5)
plot(tim_axis,fastsmooth(V2,300,1,1),'k','Linewidth',1,'Color',[0.7 0.5 0.5])
axis tight; xlim([-0.7 1.7]);%ylim([-0.3 1.2])
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_140_smopto.pdf'])
%savefig(gcf, [ savepath 'Firing_140_smopto.fig'])

%% BARPLOT FIRING
V1b=nanmean(allfiringSM(baseTimsel, stim_type==0),1);
V2b=nanmean(allfiringSM(baseTimsel, stim_type==1),1);

V1s=nanmean(allfiringSM(StimTimselTR, stim_type==0),1)-V1b;
V2s=nanmean(allfiringSM(StimTimselTR, stim_type==1),1)-V2b;
V1s2=nanmean(allfiringSM(StimTimselSU, stim_type==0),1)-V1b;
V2s2=nanmean(allfiringSM(StimTimselSU, stim_type==1),1)-V2b;

V1p=nanmean(allfiringSM(PostTimsel, stim_type==0),1);
V2p=nanmean(allfiringSM(PostTimsel, stim_type==1),1);
figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters')
M=[V1s', V1s2'];M2=[V2s', V2s2'];
b1=bar(2,nanmean(V1s),'Facecolor',[ 0 0 0.9]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(V1s2),'Facecolor',[ 0 0 0.9]);hold on,
set(b1,'FaceAlpha',0.4)
b1=bar(5,nanmean(V2s),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(6,nanmean(V2s2),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
errorbar([ 2 3 ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
errorbar([ 5 6],nanmean(M2,1), nanstd(M2)./sqrt(size(M2,1)),'.k')
ylim([-5 5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_barplotopto.pdf'])
savefig(gcf, [ savepath 'Firing_barplotopto.fig'])


%% Base vs trans
[h,p,ci,stats] = ttest(V1s) % 40Hz%
%df=19,   0.23
[h,p,ci,stats] = ttest(V2s)% 140Hz
%df=20, 0.77
%% Base vs sustained
[h,p,ci,stats] = ttest( V1s2) % 40Hz
%df=19, 0.25
[h,p,ci,stats] = ttest( V2s2)% 140Hz
%df=20,   0.16
