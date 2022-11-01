clear all


allsamp=[];M=[];MS=[];O=[];VS=[]; VS1=[];
FS=828;
stim_freq=0;
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
        load(Cpath) % Loading of the data file
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
            vsub=zscore(vsub-mean(vsub));
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
        %%%%
        %%%%5%%%
        MS1=[];  VS1=[];
        for tr=1:size( alls1,2)
            try
                opto=lfp2.trial{tr}(2,1:end);
                windS=25;
                trigs=find(diff(opto)>0.2);
                for ind=1:length(trigs)-1
                    if trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))> 5 & trigs(ind)<FS
                        M=[M,  alls1( trigs(ind)-windS:trigs(ind)+windS,tr)];
                        
                    end
                    
                    if trigs(ind)-windS>0  & trigs(ind) +windS < size(alls1,1) & (trigs(ind+1)-trigs(ind))>5 & trigs(ind)>FS  &trigs(ind)<FS*2
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

savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\revision\'
pheight=160

timax=([-windS:windS]-1).*1.2;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(timax,(nanmean(MS,2)).*FS,'k');hold on,
fill_error_area2(timax, nanmean(MS,2).*FS,   (nanstd(MS,[],2).*FS)./sqrt(size(MS,2)), [ 0.5 0.5 .5])
if stim_freq==1 
minM=0; maxM=80;
   line([ 0 0],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([ -25 -25],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([ 25 25],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    axis tight
    xlim([-10 28])
else
    minM=0; maxM=20;
   line([ 0 0],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([ -7.2 -7.2],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([7.2 7.2],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    axis tight ;xlim([-5 10])
end
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_FR_pop40hz.pdf'])

timax=([-windS:windS]-1).*1.2;
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters')
plot(timax,(nanmean(zscore(VS,[],1),2)),'k');hold on,
fill_error_area2(timax, nanmean(zscore(VS,[],1),2),   (nanstd(zscore(VS,[],1),[],2))./sqrt(size(VS,2)), [ 0.5 0.5 .5])
if stim_freq==1
    minM=-0.4; maxM=0.6;
line([ 0 0],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([ -25 -25],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([ 25 25],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    axis tight
    xlim([-10 28])
else;       minM=-0.2; maxM=0.3;
line([ 0 0],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([ -7.2 -7.2],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    line([7.2 7.2],[ minM maxM],'Color',[ 0.3 0.3 0.3],'Linestyle','--')
    axis tight ;xlim([-5 10])
end
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Puls_trig_VM_pop40hz.pdf'])

%xlim([-10 10])





