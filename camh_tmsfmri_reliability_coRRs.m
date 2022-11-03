
%Do ICC on TMS-fMRI data splitting by trial numbers
addpath('C:\Users\colin_hawco\Documents\MATLAB\SCN_Core_Support\Statistics_tools\')
addpath C:\Users\colin_hawco\Documents\MATLAB\NIFTI_toolbox-master



cd F:\CAMH_TMSfMRI\trialnums_paper\inter_permute
basedir=(pwd)

cort=32767*2; 
trials=[50 40 30 25 20 15 10]

tic
pdx =[102:105 107:110 112:114 116:121 123:126]
for trials=[ 15 10] %40 30 25 20
    for perm=1:100
        parfor idx=1:21
            camhtms_spm_permtrials(['p' num2str(pdx(idx))], perm, trials)
        end
    end
end
toc


tic
pdx =[102:105 107:110 112:114 116:121 123:126]
for perm=1:100
    parfor idx=1:21
        camhtms_spm_permtrials_50(['p' num2str(pdx(idx))], perm)
    end
end
toc


%%


cd F:\CAMH_TMSfMRI\trialnums_paper\inter_permute
basedir=(pwd)

at=load_nii('Glass_atlas_nifti.nii');
glasser=[at.img(:,1)' at.img(:,2)'];

pdx =[102:105 107:110 112:114 116:121 123:126];
trials=[50 40 30 25 20 15 10]

%55555 DO SPEARMANS NOW ICC NO GOOD
gcor=zeros(21,7,100,360);
gscor=zeros(21,7,100,360);
gcor_hi=zeros(21,7,100,360);
gcor_lo=zeros(21,7,100,360);
cors_hilo=zeros(7,21,100);
sp_cors_hilo=zeros(7,21,100);

for  tdx=1:length(trials)         %%%%%%%%1:length(trials)
    tn=trials(tdx)
    for idx=1:21
 
        for perm=1:100
            [tdx idx perm]
            cd([basedir '\trials_' num2str(tn) '\p' num2str(pdx(idx)) '\perm' num2str(perm)]);
            d=zeros(2,cort);
            
            for con=1:2
                ff=load_nii(['spmT_000' num2str(con) '.nii']);
                d(con,1:cort) = [ff.img(:,1)' ff.img(:,2)'];
            end
            
           % cors_hilo(tdx, idx,perm) =  corr(d(1,:)', d(2,:)','rows', 'pairwise');
            %sp_cors_hilo(tdx, idx,perm) = corr(d(1,:)', d(2,:)','type', 'spearman','rows', 'pairwise');
                        
            for rdx=1:360
               % gcor(idx,tdx,perm,rdx) = corr(d(1,glasser==rdx)', d(2,glasser==rdx)','type', 'pearson','rows', 'pairwise');
                gscor(idx,tdx,perm,rdx) = corr(d(1,glasser==rdx)', d(2,glasser==rdx)','type', 'spearman','rows', 'pairwise');
            end
        end
    end
end

cd F:\CAMH_TMSfMRI\trialnums_paper
save corrdat 


%% boxplots
figure; boxplot([mean(cors_hilo,3)'],  'notch', 'on', 'colors', 'k')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'YLim', [-.2 .9])

hold on; boxplot([mean(simcor,3)],  'notch', 'on', 'colors', [0.5 0.5 0.5])
set(findobj(gca,'type','line'),'linew',2)
set(gca,'YLim', [-.2 .9])
set(gca, 'FontSize', 16, 'box','off','XTickLabel',[],'XTick',[])


set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[], )


%50 trials
figure; boxplot([squeeze((cors_hilo(1,:,:)))'],  'notch', 'on', 'colors', 'k')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'YLim', [0 0.9])
set(gca, 'FontSize', 14, 'box','off','XTickLabel',[],'XTick',[])
%set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])

%30 trials
figure; boxplot([squeeze((cors_hilo(3,:,:)))'],  'notch', 'on', 'colors', 'k')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'YLim', [0 0.9])
set(gca, 'FontSize', 14, 'box','off','XTickLabel',[],'XTick',[])
%set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])


figure; boxplot([squeeze((cors_hi(1,:,:)))'],  'notch', 'on', 'colors', 'k',  'positions', [2:4:21*4])
set(findobj(gca,'type','line'),'linew',2)
hold on; boxplot([squeeze((cors_lo(1,:,:)))'],  'notch', 'on', 'colors', 'k','positions', [3:4:21*4]) 


set(gca,'YLim', [0 0.9])
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])


%%
% Write region specific ICCs on glasser atlas
% THERES A BETTER WAY TO DO THIS
out=ff;
out.img(:) = 0;
tdx=1



for idx=1:3
    dd=mean(squeeze(gcor(idx,tdx,:,:)));
    for rdx=1:360
        out.img(at.img==rdx) = dd(rdx);
    end
    out.fileprefix=['sub' num2str(idx) '_R50trials'];
    save_nii(out, ['sub' num2str(idx) '_R50trials.nii']);
end


dd=squeeze(mean(mean(squeeze(gcor(:,1,:,:)))));
for rdx=1:360
    out.img(at.img==rdx) = dd(rdx);
end
out.fileprefix=['all_R50trials'];
save_nii(out, ['all_R50trials.nii']);



% NOW FOR 30 TRIALS
out=ff;

tdx=3
for idx=1:3
    dd=mean(squeeze(gcor(idx,tdx,:,:)));
    out.img(:)=0;
    for rdx=1:360
        out.img(at.img==rdx) = dd(rdx);
    end
    out.fileprefix=['sub' num2str(idx) '_R30trials'];
    save_nii(out, ['sub' num2str(idx) '_R30trials.nii']);
end


dd=squeeze(mean(mean(squeeze(gcor(:,tdx,:,:)))));
for rdx=1:360
    out.img(at.img==rdx) = dd(rdx);
end

out.fileprefix=['all_R30trials'];
save_nii(out, ['all_R30trials.nii']);

%% histograms

tdx=1
t=gcor(:,tdx,:,:);
figure; hist(t(:),20)
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[], 'YLim', [0 165000])

tdx=3
t=gcor(:,tdx,:,:);
figure; hist(t(:),20)
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[], 'YLim', [0 165000])


%% load data from individaul tmaps
at=load_nii('F:\CAMH_TMSfMRI\trialnums_paper\inter_permute\Glass_atlas_nifti.nii');
glasser=[at.img(:,1)' at.img(:,2)'];

pdx =[102:105 107:110 112:114 116:121 123:126];
statdir = 'F:\CAMH_TMSfMRI\HCP_spm\';

for idx=1:21
    
    cd([statdir '\P'  num2str(pdx(idx))]);
    
    for con=1:5
        
        ff=load_nii(['spmT_000' num2str(con) '.nii']);
        d = [ff.img(:,1)' ff.img(:,2)'];
        p=d;
        p(d<0) = 0;
        n=d;
        n(d>0) = 0;
        
        
        for rdx=1:360
           tstats(idx,con,rdx) = mean(abs(d(glasser==rdx)));
           mtstats(idx,con,rdx) = mean((d(glasser==rdx)));
           
          
           ptstats(idx,con,rdx) = mean(nonzeros(p(glasser==rdx) ));
           ntstats(idx,con,rdx) = mean(nonzeros(n(glasser==rdx) ));
           
        end
        
    end
end


for idx=1:21
    for tdx=1:7
        [ticc_corr(idx,tdx)  picc_corr(idx,tdx) ] = corr(squeeze((tstats(idx,1,:))), squeeze(mean(gcor(idx,tdx,:,:),3)), 'rows', 'pairwise', 'type', 'pearson');
        
%         t = squeeze((mtstats(idx,1,:)));
%         g=squeeze(mean(gcor(idx,tdx,:,:),3));
%         [pticc_corr(idx,tdx)  ppicc_corr(idx,tdx) ] = corr(t(t>0), g(t>0), 'rows', 'pairwise', 'type', 'Spearman');
%         [nticc_corr(idx,tdx)  npicc_corr(idx,tdx) ] = corr(t(t<0), g(t<0), 'rows', 'pairwise', 'type', 'Spearman');
%         
        [pticc_corr(idx,tdx)  ppicc_corr(idx,tdx) ] = corr(squeeze((ptstats(idx,1,:))), squeeze(mean(gcor(idx,tdx,:,:),3)), 'rows', 'pairwise', 'type', 'pearson');
        [nticc_corr(idx,tdx)  npicc_corr(idx,tdx) ] = corr(squeeze((ntstats(idx,1,:))), squeeze(mean(gcor(idx,tdx,:,:),3)), 'rows', 'pairwise', 'type', 'pearson');
        
        %[ hi_corr(idx,tdx),  hi_p(idx,tdx) ] = corr(squeeze((tstats(idx,4,:))), squeeze(mean(gICC_hi(idx,1,:,:),3)), 'rows', 'pairwise', 'type', 'Spearman');
        %[ lo_corr(idx,tdx) lo_p(idx,tdx) ] = corr(squeeze((tstats(idx,5,:))), squeeze(mean(gICC_lo(idx,1,:,:),3)), 'rows', 'pairwise', 'type', 'Spearman');
    end
end

for idx=1:5; 
figure; plot(squeeze((tstats(idx,1,:))), squeeze(mean(gcor(idx,1,:,:),3)), 'ko')
end
[ ticc_corr; hi_corr; lo_corr]'

figure; plot(squeeze((tstats(:,1,:))), squeeze(mean(gcor(:,1,:,:),3)), 'ko')

a = reshape(squeeze(mean(gcor(:,1,:,:),3)), 21*360,1);
b = reshape(squeeze(tstats(:,1,:)), 21*360, 1);

figure; plot(b,a,'ko')

a = reshape(squeeze(mean(gcor(:,1,:,:),3)), 21*360,1);
b = reshape(squeeze(tstats(:,1,:)), 21*360, 1);
figure; plot(b,a,'ro')

b = reshape(squeeze(ntstats(:,1,:)), 21*360, 1);
figure; plot(b,a,'bo')


csvwrite('tstat_plots.csv', [b a ])

%%

cd D:\programs\workbench\cvs_avg35_inMNI152\MNINonLinear\fsaverage_LR32k\
atlas = ft_read_cifti('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
cd F:\CAMH_TMSfMRI\trialnums_paper\inter_permute
data = ft_read_cifti('sub1_30trials.dscalar.nii');

for idx=1:21
    idx
    dd=squeeze(tstats(idx,1, :));
    for rdx=1:360
        data.dscalar(atlas.indexmax==rdx) = dd(rdx); 
    end
    
   ft_write_cifti(['sub' num2str(idx) '_tstat.dscalar.nii'], data, 'parameter' , 'dscalar')
    
end





%%
% RANDOM PERMUTATION USING RESTING STATE RUN AND RANDOM TRIAL ONSETS
% this part ran analsysi and made the maps
pdx =[102:105 108:110 112:114 116:121 123:126]   % 107
for trials=[50 40 30 25 20 15 10]

        parfor idx=1:length(pdx)
            camhtms_spm_permtrials_restsim( pdx(idx), perm, trials)
        end
    end
end

% This part loads the data into matlab for further useage. 
cd F:\CAMH_TMSfMRI\trialnums_paper\inter_permute
at=load_nii('Glass_atlas_nifti.nii');
glasser=[at.img(:,1)' at.img(:,2)'];
cort=32767*2; 

pdx =[102:105 108:110 112:114 116:121 123:126]; %107

% gdat = zeros(21,length(trials),6,100,360);
simcor=zeros(length(pdx),length(trials),100);

trials=[50 40 30 25 20 15 10]
for tdx=1:length(trials)
    tn=trials(tdx)
    for idx=1:length(pdx)
        for perm=1:100
            [tdx idx perm]
            cd( ['D:\trialnums_perms\rest\trials' num2str(tn) '\pt' num2str(pdx(idx)) '\perm' num2str(perm) '\']);
            d=zeros(2,cort);
            
            for con=1:2
                ff=load_nii(['spmT_000' num2str(con) '.nii']);
                d(con,1:cort) = [ff.img(:,1)' ff.img(:,2)'];
            end
            
            simcor(idx,tdx,perm) = corr(d(1,glasser==rdx)', d(2,glasser==rdx)','type', 'pearson','rows', 'pairwise');
             
            for rdx=1:360
                gSimcor(idx,tdx,perm,rdx) = corr(d(1,glasser==rdx)', d(2,glasser==rdx)','type', 'pearson','rows', 'pairwise');
            end
             
        end
    end
end

% we need to make boxplots, as we did for the actual data. 
figure; boxplot([mean(simcor,3)],  'notch', 'on')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'YLim', [0.1 0.9])
 
cd F:\CAMH_TMSfMRI\trialnums_paper
save corrdat 


tdx=1
t=gSimcor(:,tdx,:,:);
figure; hist(t(:),20)
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[], 'Ylim', [0 165000])

tdx=3
t=gSimcor(:,tdx,:,:);
figure; hist(t(:),20)
set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[], 'Ylim', [0 165000])


%% check histographs for tstats

for idx=1:21
    
    cd([statdir '\P'  num2str(pdx(idx))]);
   
        ff=load_nii(['spmT_0001.nii']);
        d = [ff.img(:,1)' ff.img(:,2)'];
        
        figure; hist(d, 20)
        sk(idx)=skewness(d);
        
end

            
            
            
