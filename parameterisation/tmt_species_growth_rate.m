% Calculate the mean growth rates for different species at grid cell level, based on single-species dominated stands
% Uses NFI data from the TreeMort dataset
%
% T. Pugh
% 12.08.24

% module load slurm-interactive
% fisbatch_tmux --qos bbdefault -n 1 -t 2:0:0
% module load bear-apps/2022b
% module load MATLAB/2023b-r6
% matlab -nodisplay -nosplash

use_neighbourhood=false;

thres_spec=0.75; % Fraction of basal area from one species to be consisdered single species plot

specnames={'Abies alba','Betula pendula','Betula pubescens','Carpinus betulus','Fagus sylvatica','Fraxinus excelsior',...
    'Larix decidua','Picea abies','Pinus sylvestris','Pinus halpensis','Populus tremula','Quercus ilex','Quercus pubescens',...
    'Quercus robur'};
nspec=length(specnames);

%datasets={'NSW','NFG','NFL','NCZ','NNL','NPO','NSI','NSP','NFR'}; % MISSING WALLONIA AND FINLAND
datasets={'FUN','NCZ','NFG','NFR','NNL','NSI','NSW','NPO'}; % Susanne's list for evaluation dataset (minus NOR)
ndatasets=length(datasets);
census_to_use=[2,2,3,2,7,4,3,2]; % Same as used by Susanne

% Grid locations to make calculations for
gridlist=readtable('/Users/pughtam/Documents/LPJ-GUESS/Fix_Europe/gridlist_fix_europe_rand_plus_norway_26gc.txt');
ngrid=height(gridlist);

% Initialise the output arrays
WBgrowth=zeros(ndatasets,ngrid,nspec);
WBgrowth_nplot=zeros(ndatasets,ngrid,nspec);

for dd=1:ndatasets

    fprintf('dataset is %s\n',datasets{dd})

    % Read plot data
    %format_dir='/rds/projects/2018/pughtam-treeinv/03_TreeMort-database/01_Format-data/output/';
    %process_dir='/rds/projects/2018/pughtam-treeinv/03_TreeMort-database/03_Query-data/output/data/';
    format_dir='/Users/pughtam/data/Plots/TreeMort_data/';
    process_dir='/Users/pughtam/data/Plots/TreeMort_data/';
    treetable_file=[process_dir,'/',datasets{dd},'/01_treedata-biomass_TMt_',datasets{dd},'.csv'];
    plottable_file=[format_dir,'/',datasets{dd},'/04_plot-info_TMt_',datasets{dd},'.csv'];
    sdtable_file=[process_dir,'/',datasets{dd},'/02_stand-level-dynamics_TMt_',datasets{dd},'.csv'];

    TMtt=readtable(treetable_file,'TreatAsMissing','NA');
    TMtp=readtable(plottable_file,'TreatAsMissing','NA');
    TMts=readtable(sdtable_file,'TreatAsMissing','NA');

    nplot=height(TMtp);
    ncensus=max(unique(TMts.census_n));

    % Select only from most recent census
    %TMtt_c=TMtt(TMtt.census_n==ncensus,:);
    %TMts_c=TMts(TMts.census_n==ncensus,:);
    % Select from predetermined census
    TMtt_c=TMtt(TMtt.census_n==census_to_use(dd),:);
    TMts_c=TMts(TMts.census_n==census_to_use(dd),:);

    % Identify if plots have dominant species and if so which
    spec_ind=zeros(nplot,1);
    spec_n=zeros(nplot,1);
    for nn=1:nplot
        aa=find(strcmp(TMtt_c.tmt_plot_id,TMtp.tmt_plot_id{nn}));
        if isempty(aa)
            continue
        end
        ba_plot=sum(TMtt_c.ba(aa));

        % Now test if any of the species dominate the plot
        for ss=1:nspec
            bb=find(strcmp(TMtt_c.species_cor(aa),specnames{ss}));
            ba_spec=sum(TMtt_c.ba(aa(bb)));
            if ba_spec/ba_plot > thres_spec
                spec_ind(nn)=ss;
                spec_n(nn)=length(aa); % Total number of individuals on the plot
                break
            end
        end
        if mod(nn,1000)==0
            fprintf('Processed %d of %d\n',nn,nplot)
        end
    end
    clear nn aa bb ss ba_plot ba_spec

    % Now calculate the mean growth rate for each species for each grid cell

    TMtp.longitude_0p5=(floor(TMtp.longitude*2)/2)+0.25;
    TMtp.latitude_0p5=(floor(TMtp.latitude*2)/2)+0.25;

    % WBgrowth=zeros(ngrid,nspec);
    % WBgrowth_nplot=zeros(ngrid,nspec);
    for gg=1:ngrid
        aa=find(TMtp.longitude_0p5==gridlist.Var1(gg) & TMtp.latitude_0p5==gridlist.Var2(gg));
        if isempty(aa)
            continue
        end
        for bb=1:length(aa)
            for ss=1:nspec
                if spec_ind(aa(bb))==ss
                    tmp=TMts_c.NPP2(strcmp(TMts_c.tmt_plot_id,TMtp.tmt_plot_id(aa(bb))));
                    WBgrowth(dd,gg,ss)=WBgrowth(dd,gg,ss) + tmp(end); % Using the tmp variable is a work around to avoid instances where a plot comes up with two entries for the same census (seen in NSI)
                    WBgrowth_nplot(dd,gg,ss)=WBgrowth_nplot(dd,gg,ss)+1;
                end
            end
        end
        WBgrowth(dd,gg,:) = WBgrowth(dd,gg,:)./WBgrowth_nplot(dd,gg,:); % Mean across all the plots for each species
    end
    clear gg aa bb tmp

end
clear dd

% Convert units from ton DM ha-1 to kg C m-2
WBgrowth=WBgrowth/10/2;

save -v7.3 WBgrowth_per_spec.mat WBgrowth WBgrowth_nplot

% Now merge the different datasets into one map per species, combining where there is overlap on country borders

WBgrowth_overall=squeeze(sum(WBgrowth.*WBgrowth_nplot,1));
WBgrowth_nplot_overall=squeeze(sum(WBgrowth_nplot,1));
WBgrowth_overall=WBgrowth_overall./WBgrowth_nplot_overall;

save WBgrowth_per_spec_overall.mat WBgrowth_overall WBgrowth_nplot_overall

%% Calculate neighbourhood growth rate and neighbourhood nplot
% PERHAPS MAKE THIS OPTIONAL - SEE IF WE NEED IT

WBgrowth_overall_neighbourhood=NaN(size(WBgrowth_overall));
WBgrowth_nplot_overall_neighbourhood=NaN(size(WBgrowth_nplot_overall));
for gg=1:ngrid
    aa=find(gridlist.Var1>=(gridlist.Var1(gg)-0.5) & gridlist.Var1<=(gridlist.Var1(gg)+0.5) &...
        gridlist.Var2>=(gridlist.Var2(gg)-0.5) & gridlist.Var2<=(gridlist.Var2(gg)+0.5));
    WBgrowth_nplot_overall_neighbourhood(gg,:)=nansum(WBgrowth_nplot_overall(aa,:),1);
    WBgrowth_overall_neighbourhood(gg,:)=nansum(WBgrowth_overall(aa,:).*WBgrowth_nplot_overall(aa,:),1);
    WBgrowth_overall_neighbourhood(gg,:)=WBgrowth_overall_neighbourhood(gg,:)./WBgrowth_nplot_overall_neighbourhood(gg,:);
end

if use_neighbourhood
    WBgrowth_overall_select=WBgrowth_overall_neighbourhood;
    WBgrowth_nplot_overall_select=WBgrowth_nplot_overall_neighbourhood;
    fprintf('Using neighbourhood statistics\n')
else
    WBgrowth_overall_select=WBgrowth_overall;
    WBgrowth_nplot_overall_select=WBgrowth_nplot_overall;
    fprintf('Using single grid cell statistics\n')
end

%% Make some maps of the results by species

% Number of plots per gridcell
figure
for ss=1:nspec
    subplot(3,5,ss)
    s1=scatter(gridlist.Var1,gridlist.Var2,8,WBgrowth_nplot_overall_select(:,ss));
    set(s1,'MarkerFaceColor','flat')
    title(specnames{ss})
    colorbar
    caxis([0 10])
end
clear ss

% Growth rates
figure
for ss=1:nspec
    subplot(3,5,ss)
    s1=scatter(gridlist.Var1,gridlist.Var2,8,WBgrowth_overall_select(:,ss));
    set(s1,'MarkerFaceColor','flat')
    title(specnames{ss})
    colorbar
    caxis([0 0.5])
end
clear ss

% Growth rates filtered by having at least 10 underlying plots per grid cell
ind_plot=WBgrowth_nplot_overall_select>=10;
figure
for ss=1:nspec
    subplot(3,5,ss)
    s1=scatter(gridlist.Var1(ind_plot(:,ss)),gridlist.Var2(ind_plot(:,ss)),8,WBgrowth_overall_select(ind_plot(:,ss),ss));
    set(s1,'MarkerFaceColor','flat')
    title(specnames{ss})
    colorbar
    caxis([0 0.5])
    set(gca,'XLim',[-10 33])
    set(gca,'YLim',[30 70])
end
clear ind_plot ss
% NOTE: CONCLUSION FROM THE ANALYSIS AT THIS POINT IS THAT GRIDCELL LEVEL ALONE IS NOT ENOUGH FOR MOST SPECIES

% Plot neighbourhood and normal versions against each other
figure
for ss=1:nspec
    subplot(3,5,ss)
    hold on
    plot(WBgrowth_overall(:,ss),WBgrowth_overall_neighbourhood(:,ss),'.')
    ind_plot=WBgrowth_nplot_overall_select>=10;
    plot(WBgrowth_overall(ind_plot(:,ss),ss),WBgrowth_overall_neighbourhood(ind_plot(:,ss),ss),'.')
    xlabel('Standard growth rate')
    ylabel('Neighbour growth rate')
    plot([0 1],[0 1],'k')
    title(specnames{ss})
    set(gca,'XLim',[0 1])
    set(gca,'YLim',[0 1])
end
clear ss
% NOTE: CONCLUSION IS THAT NEIGHBOURHOOD GROUPING IS SIGNIFICANTLY REDUCING SCATTER

figure
for ss=1:nspec
    subplot(3,5,ss)
    hold on
    plot(WBgrowth_nplot_overall_select(:,ss),WBgrowth_overall_select(:,ss),'.')
    title(specnames{ss})
end
clear ss

% Compare against the overall growth rate from the gridcell, as calculated elsewhere
% Note: only really a fair test if compared with gridcell level stats
grid_stats=readtable('/Users/pughtam/Documents/LPJ-GUESS/Fix_Europe/tmt_grid_4evaluation_03-01-24/tmt_grid_4evaluation_03-01-24.txt');

npp2=NaN(ngrid,1);
for gg=1:ngrid
    aa=find(grid_stats.Lon==gridlist.Var1(gg) & grid_stats.Lat==gridlist.Var2(gg));
    if isempty(aa)
        continue
    end
    npp2(gg)=grid_stats.npp2_mean(aa);
end
clear gg aa

npp2=npp2/10/2; % Unit conversion

ind_plot=WBgrowth_nplot_overall_select>=10;
figure
for ss=1:nspec
    subplot(3,5,ss)
    hold on
    plot(WBgrowth_overall_select(ind_plot(:,ss),ss),npp2(ind_plot(:,ss)),'.')
    plot([0 1],[0 1],'k')
    title(specnames{ss})
    set(gca,'XLim',[0 1])
    set(gca,'YLim',[0 1])
    xlabel('Species')
    ylabel('Gridcell')
end
clear ss ind_plot
% NOTE: CONCLUSION IS THAT THERE IS A PRETTY GOOD CORRESPONDENCE WHICH GENERALLY DEVIATES IN EXPECTED WAYS (E.G. SPRUCE 
% HIGHER THAN GRIDCELL AVERAGE) AND THE CALCULATIONS IN THE SCRIPT ARE THEREFORE LIKELY SOLID
% HOWEVER THE ADDITIONAL DISASSOCIATION FROM THE GRIDCELL BY USING THE NEIGHBOURHOOD MEANS THAT IT MIGHT BE BETTER TO
% TAKE ONLY GRID CELL LEVEL RESULTS, EVEN WHERE THESE ARE FEWER. BECAUSE THEY GENERALLY SHOW CONSISTENT BEHAVIOUR
% RELEVATIVE TO THE GRIDCELL AVERAGE, THEY ARE LIKELY SOLID

% OVERALL CONCLUSION: DO NOT APPLY THE NEIGHBOURHOOD RULE.

%% Create a list of grid cells to run a calibration on

% Calculate the number of gridcells per species with the number of plots above a threshold

nsample_grids=96; % Total number of gridcells to sample

% Number of gridcells to sample purely based on latitude (one per degree of latitude)
minlat=prctile(gridlist.Var2,1);
maxlat=prctile(gridlist.Var2,99);
nsample_grids_lats=maxlat-minlat+1; 

thres_nplot=10;

ngrid_spec_thres=sum(WBgrowth_nplot_overall_select>=thres_nplot);
% Use this to calculate a fraction to sample which is relative to the spatial abundance of each species
ngrid_spec_thres_tot=sum(ngrid_spec_thres); % Total of all gridcell-species combinations above the threshold
ngrid_spec_thres_frac=ngrid_spec_thres/ngrid_spec_thres_tot;

% First set a minimum of 3 gridcells per species and exclude any species that don't meet this criteria
grid_spec_sample=zeros(1,nspec);
for ss=1:nspec
    if ngrid_spec_thres(ss)<3
        fprintf('Not including %s because less than 3 grid cells meet the threshold\n',specnames{ss});
        continue
    end
    grid_spec_sample(ss)=3;
end
clear ss
% Then allocate the remaining gridcells so as to get more samples for the most common species
left_to_allocate=nsample_grids-nsample_grids_lats-sum(grid_spec_sample);
grid_spec_sample=grid_spec_sample+round(ngrid_spec_thres_frac.*left_to_allocate);


% Randomly select gridcells to sample for each species, excluding the gridcells with more extreme values
ngrid_sel=false(ngrid,nspec+1);
samp_size=NaN(nspec,1);

% First sample based on species (it is possible that some gridcells will be selected for more than one species. This is OK)
for ss=1:nspec
    perc_10th=prctile(WBgrowth_overall_select(:,ss),10);
    perc_90th=prctile(WBgrowth_overall_select(:,ss),90);
    aa=find(WBgrowth_nplot_overall_select(:,ss)>=thres_nplot & WBgrowth_overall_select(:,ss)>=perc_10th...
        & WBgrowth_overall_select(:,ss)<=perc_90th);
    if isempty(aa)
        fprintf('Warning: no values for species %s\n',specnames{ss})
        continue
    end
    rand_ind=randsample(length(aa),grid_spec_sample(ss));
    ngrid_sel(aa(rand_ind),ss)=true;
end
clear ss aa

% Then based on latitude (referenced in final index in ngrid_sel)
for ll=1:nsample_grids_lats
    lat_curr=minlat+ll-1;
    aa=find(gridlist.Var2==lat_curr & sum(ngrid_sel,2)==0); % Last check here is to avoid double sampling.
    if isempty(aa)
        continue
    end
    rand_ind=randsample(length(aa),1);
    ngrid_sel(aa(rand_ind),end)=true;
end

% Clean up the arrays to only include the selected gridcells for output
ngrid_sel_clean=ngrid_sel;
ngrid_sel_clean(sum(ngrid_sel,2)==0,:)=[];
gridlist_clean=gridlist;
gridlist_clean(sum(ngrid_sel,2)==0,:)=[];

WBgrowth_overall_select_clean=zeros(ngrid,nspec);
for ss=1:nspec
    for gg=1:ngrid
        if ngrid_sel(gg,ss)
            WBgrowth_overall_select_clean(gg,ss)=WBgrowth_overall_select(gg,ss);
        end
    end
end
clear ss gg
WBgrowth_overall_select_clean(sum(ngrid_sel,2)==0,:)=[];

% Plot selected gridcells

figure
subplot(4,5,1)
hold on
s1=scatter(gridlist.Var1,gridlist.Var2,8,npp2); % Uses index aa from loop above
set(s1,'MarkerFaceColor','flat')
caxis([0 0.5])
plot(gridlist_clean.Var1,gridlist_clean.Var2,'rx')
title('Overall')
colorbar
XLim=get(gca,'XLim');
YLim=get(gca,'YLim');

subplot(4,5,2)
hold on
s1=scatter(gridlist.Var1,gridlist.Var2,8,npp2); % Uses index aa from loop above
set(s1,'MarkerFaceColor','flat')
caxis([0 0.5])
plot(gridlist_clean.Var1(ngrid_sel_clean(:,end)),gridlist_clean.Var2(ngrid_sel_clean(:,end)),'rx')
title('Latitude')
colorbar
XLim=get(gca,'XLim');
YLim=get(gca,'YLim');
for ss=1:nspec
    subplot(4,5,ss+2)
    hold on
    % Line below for aa needs to be the same as in loop above for selecting the gridcells
    aa=find(WBgrowth_nplot_overall_select(:,ss)>=thres_nplot);
    s1=scatter(gridlist.Var1(aa),gridlist.Var2(aa),8,WBgrowth_overall_select(aa,ss));
    plot(gridlist.Var1(ngrid_sel(:,ss)),gridlist.Var2(ngrid_sel(:,ss)),'rx')
    set(s1,'MarkerFaceColor','flat')
    title(specnames{ss})
    colorbar
    caxis([0 0.5])
    set(gca,'XLim',XLim,'YLim',YLim)
end
clear ss aa

% Write out selected gridcells to file
specnames_out=specnames; specnames_out{nspec+1}='Latitude';
ngrid_sel_clean_table_1=array2table(ngrid_sel_clean,'VariableNames',specnames_out);
ngrid_sel_clean_table_2=array2table([gridlist_clean.Var1,gridlist_clean.Var2],'VariableNames',{'Lon','Lat'});
ngrid_sel_clean_table=cat(2,ngrid_sel_clean_table_2,ngrid_sel_clean_table_1);
clear ngrid_sel_clean_table_1 ngrid_sel_clean_table_2
writetable(ngrid_sel_clean_table,'gridlist_calibration_species.csv')

WBgrowth_clean_table_1=array2table(WBgrowth_overall_select_clean,'VariableNames',specnames);
WBgrowth_clean_table_2=array2table([gridlist_clean.Var1,gridlist_clean.Var2],'VariableNames',{'Lon','Lat'});
WBgrowth_clean_table=cat(2,WBgrowth_clean_table_2,WBgrowth_clean_table_1);
clear WBgrowth_clean_table_1 WBgrowth_clean_table_2
writetable(WBgrowth_clean_table,'WBgrowth_calibration_species.csv')

% Save the file so that we can come back to the same outputs later (needed because of the random selection)
save -v7.3 calibration_gridlist_calculation.mat



