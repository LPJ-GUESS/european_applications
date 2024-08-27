% Update trait values for LPJ-GUESS European PFTS based on information from
% - TRY request 6033 for Specific Leaf Area (SLA)
% - The Global Wood Density Database (Zanne et al., 2009, https://datadryad.org/stash/dataset/doi:10.5061/dryad.234) for
%
% T. Pugh
% 30.07.24

datafol='/Users/pughtam/Documents/LPJ-GUESS/Fix_Europe/traits/';

pftname={'Abi alb','Bet pen','Bet pub','Car bet','Cor ave','Fag syl','Fra exc','Jun oxy',...
    'Lar dec','Pic abi','Pic sit','Pin syl','Pin hal','Pop tre','Que coc','Que ile',...
    'Que pub','Que rob','Til cor','Ulm gla'}; % Names of PFTs in LPJ-GUESS

specname={'Abies alba','Betula pendula','Betula pubescens','Carpinus betulus','Corylus avellana','Fagus sylvatica','Fraxinus excelsior','Juniperus oxycedrus',...
    'Larix decidua','Picea abies','Picea sitchensis','Pinus sylvestris','Pinus halepensis','Populus tremula','Quercus coccifera','Quercus ilex',...
    'Quercus pubescens','Quercus robur','Tilia cordata','Ulmus glabra'}; % Species names corresponding to PFTs
nspec=length(specname);

%% SLA

tt=readtable([datafol,'/6033.txt']);

% Filter for species
spec_sla_median=NaN(nspec,1);
spec_sla_mean=NaN(nspec,1);
ndata_sla=NaN(nspec,1);
for ss=1:nspec
    spec_ind=find(strcmp(tt.SpeciesName,specname{ss}));

    tt_spec=tt(spec_ind,:);

    % Filter for SLA - take all measurement types
    sla_ind=find(strcmp(tt_spec.TraitName,'Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included'));
    sla_ind2=find(strcmp(tt_spec.TraitName,'Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded'));
    sla_ind3=find(strcmp(tt_spec.TraitName,'Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded'));
    
    tt_spec_sla=tt_spec(cat(1,sla_ind,sla_ind2,sla_ind3),:);

    % Filter for correct units
    unit_ind=find(strcmp(tt_spec_sla.UnitName,'mm2 mg-1'));

    tt_spec_sla_unit=tt_spec_sla(unit_ind,:);

    spec_sla_median(ss)=median(tt_spec_sla_unit.StdValue(isfinite(tt_spec_sla_unit.StdValue)))*2; % *2 to convert to carbon units
    spec_sla_mean(ss)=mean(tt_spec_sla_unit.StdValue(isfinite(tt_spec_sla_unit.StdValue)))*2; % *2 to convert to carbon units
    ndata_sla(ss)=sum(isfinite(tt_spec_sla_unit.StdValue));

    fprintf('%s %7.3f %7.3f %d\n',pftname{ss},spec_sla_median(ss),spec_sla_mean(ss),ndata_sla(ss))
end
clear ss tt_spec tt_spec_sla tt_spec_sla_unit

%% Wood density
ww=readtable([datafol,'/GlobalWoodDensityDatabase_dotdecimal.csv']);

spec_wd_median=NaN(nspec,1);
spec_wd_mean=NaN(nspec,1);
ndata_wd=NaN(nspec,1);
for ss=1:nspec
    spec_ind=find(strcmp(ww.Binomial,specname{ss}));

    ww_spec=ww(spec_ind,:);

    spec_wd_median(ss)=median(ww_spec.WoodDensity_g_cm_3__OvenDryMass_freshVolume(isfinite(ww_spec.WoodDensity_g_cm_3__OvenDryMass_freshVolume)))/2; % /2 to convert to carbon units
    spec_wd_mean(ss)=mean(ww_spec.WoodDensity_g_cm_3__OvenDryMass_freshVolume(isfinite(ww_spec.WoodDensity_g_cm_3__OvenDryMass_freshVolume)))/2; % /2 to convert to carbon units
    ndata_wd(ss)=sum(isfinite(ww_spec.WoodDensity_g_cm_3__OvenDryMass_freshVolume));

    fprintf('%s %7.3f %7.3f %d\n',pftname{ss},spec_wd_median(ss),spec_wd_mean(ss),ndata_wd(ss))
end
clear ss spec_ind ww_spec