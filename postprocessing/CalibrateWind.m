%%% Script to read LPJ-GUESS output of wind damage in severe storm,
%%% combine it with wind damage probability, agregate it to country level,
%%% calibrate it against observations and present the results in table, figure
%%% and maps

clear
close all hidden

addpath('C:\Program Files\MATLAB\R2018b\resources\borders');
% https://climate.copernicus.eu/visualising-and-processing-climate-data-within-matlab

path = 'C:\Users\Fredrik.Lagergren\Documents\ClimbFor\european_applications\run_240809_sim4_age_probharv_satdist07\';
%path = 'C:\Users\Fredrik.Lagergren\Documents\ClimbFor\european_applications\run_240809_sim5_age_probharv_satdist07_inclcc\';


%%% READING LUND-USE FILE INTO MAPPING FORMAT %%%
lulist = dlmread('C:\Users\Fredrik.Lagergren\Documents\ClimbFor\european_applications\LC_europe_nat_for_1801_2010_Pucher_noNatural.txt','\t',1,0);
% Format "Lon	Lat	year	NATURAL	FOREST	BARREN" Fraction of land

lulist(:,1:2) = lulist(:,1:2) - 0.25; % Change from center to SW corner coordinates
years = 210; % 1801-2010
gridcells = size(lulist,1)/years; % Number of gridcells in LU list
lulist = reshape(lulist,[years gridcells 6]);

% Find extent of the LU list and create a complete map
minlon = min(permute(lulist(1,:,1),[2 3 1]));
maxlon = max(permute(lulist(1,:,1),[2 3 1]));
minlat = min(permute(lulist(1,:,2),[2 3 1]));
maxlat = max(permute(lulist(1,:,2),[2 3 1]));
lonsLU = (minlon:0.5:maxlon)';
latsLU = (minlat:0.5:maxlat)';
lumap = NaN(size(latsLU,1),size(lonsLU,1),4);
% (lat,lon,variable (1=NATURAL fraction, 2=FOREST frac., 3=BARREN frac, 4=Age

% Variables for calculating age in 2010 based on difference in FOREST fraction
agetemp = NaN(210,1); % Fraction of FOREST in age classes
ages = (210:-1:1)'; % Age in the age classes

% Populate the complete LU map file
for gridcell = 1:gridcells
    lonnr = (lulist(1,gridcell,1) - minlon) * 2 + 1; % Number in x for 0.5 degree grid
    latnr = (lulist(1,gridcell,2) - minlat) * 2 + 1; % Number in y for 0.5 degree grid
    lumap(latnr,lonnr,1:3) = lulist(210,gridcell,4:6); % Use data from 2010
    
    % Calculating age based on difference in managed forest area from year to year
    agetemp(1) = lulist(1,gridcell,5) / (1 - lulist(1,gridcell,6));
    agetemp(2:210) = (lulist(2:210,gridcell,5) - lulist(1:209,gridcell,5)) / (1- lulist(1,gridcell,6));
    lumap(latnr,lonnr,4) = sum(agetemp .* ages);
end
%%%


%%% READING OF WIND DAMAGE PROBABILITY INTO MAPPING FORMAT %%%
% Original file (WDP_10km.tif) has a resolution of 10 km and the projection
% is epsg:3035 (Lamberts Equal Area). It was converted to a 0.5 degree
% resolution txt-file in ArgGIS Pro 
wdplist = dlmread('C:\Users\Fredrik.Lagergren\Documents\FORECO\CorneliusWDPmap\wdp_05deg_from01deg_colinear.txt','\t',1,0);
% Format "Lon Lat WindDamageProbability" Annual probability of a severe storm based on wind data from 1986-2020

wdplist(:,1:2) = wdplist(:,1:2) - 0.25; % Change from center to SW corner coordinates
gridcells = size(wdplist,1); % Number of gridcells in WDP list

% Find extent of the WDP list and create a complete map
minlon = min(wdplist(:,1));
maxlon = max(wdplist(:,1));
minlat = min(wdplist(:,2));
maxlat = max(wdplist(:,2));
lonsWDP = (minlon:0.5:maxlon)';
latsWDP = (minlat:0.5:maxlat)';
wdpmap = NaN(size(latsWDP,1),size(lonsWDP,1));

for gridcell = 1:gridcells
    lonnr = (wdplist(gridcell,1) - minlon) * 2 + 1;
    latnr = (wdplist(gridcell,2) - minlat) * 2 + 1;
    wdpmap(latnr,lonnr) = wdplist(gridcell,3);
end
%%%


%%% READING OF STORM DAMAGE IN SEVERE STORM FROM LPJ-GUESS OUTPUT INTO MAP FORMAT
thefile = sprintf('%s',path,'storm.out');
stormlist = dlmread(thefile,'',1,0);
% Format "Lon      Lat    Year       SI  DamFrac DamWoodC"

years = 122; % 1901-2022
gridlist = stormlist(:,1:2);
gridcells = size(gridlist,1) / years;
gridlist = reshape(gridlist,[years gridcells 2]);
gridlist = permute(gridlist(1,:,:),[2 3 1]);
gridlist = gridlist -0.25; % Change from center to SW corner coordinates

% Find extent of the STORM list and create a complete map
minlon = min(gridlist(:,1));
maxlon = max(gridlist(:,1));
minlat = min(gridlist(:,2));
maxlat = max(gridlist(:,2));
lons = (minlon:0.5:maxlon)';
lats = (minlat:0.5:maxlat)';

stormdata = stormlist(:,4:6);
stormdata = reshape(stormdata,[years gridcells 3]);
stormmap = NaN(size(lats,1),size(lons,1),3);


for gridcell = 1:gridcells
    lonnr = (gridlist(gridcell,1) - minlon) * 2 + 1;
    latnr = (gridlist(gridcell,2) - minlat) * 2 + 1;
    stormmap(latnr,lonnr,1:3) = mean(stormdata(86:120,gridcell,:));
end
%%%


%%% CLIP THE MAPS TO HAVE THE SAME EXTENT
lonstart = 1 + (lons(1) - lonsLU(1)) * 2;
latstart = 1 + (lats(1) - latsLU(1)) * 2;
lumap = lumap(latstart:latstart+size(lats,1)-1,lonstart:lonstart+size(lons,1)-1,:);

lonstart = 1 + (lons(1) - lonsWDP(1)) * 2;
latstart = 1 + (lats(1) - latsWDP(1)) * 2;
wdpmap = wdpmap(latstart:latstart+size(lats,1)-1,lonstart:lonstart+size(lons,1)-1);

damagemap_uncall = stormmap(:,:,3) .* wdpmap; % Uncallibrated map of modelled damage (kg C m-2 yr-1) 
%%%


%%% CALLIBRATION
% Callibration of modelled damage summed up to country level 1986-2020 against data from the DFDE database
countrycode = [
    'ES'	%	Spain
    'FR'	%	France
    'BE'	%	Belgium
    'NL'	%	Netherlands
    'CH'	%	Switzerland
    'SE'	%	Sweden
    'CZ'	%	Czech Republic
    'DE'	%	Germany
    'PL'	%	Poland
    'FI'	%	Finland
%    'NO'	%	Norway
    ];
countries = size(countrycode,1);
areapath = 'C:\Users\Fredrik.Lagergren\Documents\ClimbFor\european_applications\GridcellFractionsEMEP\';

% Calculate a factor to convert carbon mass (kg C total / m2) to volume (m3 / ha)
kgtoton = 0.001; % Ton per kg (ton / kg)
stemfrac = 0.65; % Fraction of biomass that is stem (ton stem / ton total)
cfrac = 0.5; % Carbon fraction of biomass (ton C / ton total)
dens = 0.4; % Density of of wood (ton biomass / m3)
m2toha = 0.0001; % Hectar per m2 (ha / m2)
cmtovol = kgtoton * stemfrac / cfrac / dens / m2toha; % Conversion factor

%Variables to store country totals
mod_totdam_uncall = NaN(countries,1); %Totalled modelled damage 1986-2020 by country
forestarea = NaN(countries,1); %Total forest area by country

% Loop over countries to sum up the modelled total damage
for country = 1:countries
    
    % Read country list of fraction of 0.1 degree grids that is within the
    % country from EMEP
    % https://www.ceip.at/the-emep-grid/grid-definiton
    thefile = sprintf('%s',areapath,countrycode(country,:),'.csv');
    aflist = dlmread(thefile,';',1,2); % (long,lat,fraction)
    gridcells = size(aflist,1);
    
    % Change coordinates to SW corner in 0.5 degree grid
    aflist(:,1:2) = floor(aflist(:,1:2)*2)/2;
    
    % Create a complete area map with the extent of the WDP map
    areamap = zeros(size(wdpmap));
    
    for gridcell = 1:gridcells
        lonnr = 1 + (aflist(gridcell,1) - lons(1)) * 2;
        latnr = 1 + (aflist(gridcell,2) - lats(1)) * 2;
        if lonnr > 0 && latnr > 0 && lonnr <= size(lons,1) && latnr <= size(lats,1)
            latC = aflist(gridcell,2) + 0.25; % Center grid cell latitude
            acell = -28.562 * latC^2 - 1229.7 * latC + 331456; % Empirical area (ha) of 0.5 degree cell by latitude
            areamap(latnr,lonnr) = areamap(latnr,lonnr) + aflist(gridcell,3) * acell / 25;
        end
    end
    
    totdammap = stormmap(:,:,3)*cmtovol*35 .* wdpmap .* areamap .* lumap(:,:,2); % Total damage over 35 years (m3 per gridcell)
    mod_totdam_uncall(country,1) = sum(sum(totdammap,'omitnan'),'omitnan');
    forestarea(country,1) = sum(sum(areamap .* lumap(:,:,2),'omitnan'),'omitnan');
   
end

%Expert gapfilled total damage by wind (m3 per ha forest) by country 1986-2019 from the DFDE database
% as analazed by Patacca et al. (2022)
% https://doi.org/10.1038/s41586-021-03292-x
DFDE1986_2019 = [
	0.311712469	%	Spain
	15.92983022	%	France
	8.624930994	%	Belgium
	1.83575288	%	Netherlands
	18.95101885	%	Switzerland
	5.811679096	%	Sweden
	38.41414753	%	Czech Republic
	17.37980347	%	Germany
	6.684530283	%	Poland
	1.54078413	%	Finland
%	0.678861877	%	Norway
];

%Wind damage per country in 2020 directly taken from the DFDE database (milj m3)
%https://dfde.efi.int/db/dfde_app.php (accessed 7 Aug 2024, search terms: Period start 2020; Period end 2020; cause wind)
DFDE2020= [
	0.0000	%	Spain
	0.0000	%	France
	0.0000	%	Belgium
	0.0000	%	Netherlands
	0.5906	%	Switzerland
	0.0000	%	Sweden
	0.0000	%	Czechia
	1.7410	%	Germany
	0.0000	%	Poland
	0.0000	%	Finland
%	0.0000	%	Norway
];

DFDEdam = DFDE1986_2019 .* forestarea + DFDE2020 * 1000000; % Total DFDE damage 1986-2020 (m3) by country

% Linear regression through origin
xdata = mod_totdam_uncall / 1000000;
ydata = DFDEdam / 1000000;
dlm = fitlm(xdata,ydata,'Intercept',false);
slope = dlm.Coefficients.Estimate;
R2 = dlm.Rsquared.Adjusted;

damagemap_call = damagemap_uncall * slope; % Callibrated map of modelled damage (kg C m-2 yr-1)
stormmap_call = stormmap * slope; % All outputs in storm.out are linearly scaled by a calibration factor

% Make a table of total modelled and reported damage 1986-2020 for export to Excel or Word
country_table = NaN(countries,3);
country_table(:,1) = ydata; % Reported damage
country_table(:,2) = xdata; % Modelled damage before callibration
country_table(:,3) = xdata * slope; % Modelled damage after callibration
country_table = country_table'; % Transpose
%%%


%%% SCATTERPLOT OF MODELLED AND REPORTED DAMAGE
% Plot the modelled total damage 1986-2020 by country against DFDE reported
% values and the fitted equation with parameter and R2
figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 15.0 11.0]);
axel1 = axes('Position',[0.07 0.09 0.91 0.89],'Box','on','XGrid','on','YGrid','on',...
    'FontSize',20,'Parent',figure1);
xlabel(axel1,'Total modelled damage (milj m^3)');
ylabel(axel1,'Total reported damage (milj m^3)');
hold(axel1,'all');

% Plot modelled damage vs reported
plot(axel1,xdata, ydata,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',12);
hold(axel1,'all');

% Show country label for the model vs reported pairs
text(axel1,xdata+1.5,ydata,countrycode, 'Vert','middle', 'Horiz','left', 'FontSize',20);
hold(axel1,'all');

% Plot the equation
plot(axel1,[min(xdata) max(xdata)],[min(xdata)*slope max(xdata)*slope],...
    'Color','k','LineWidth',2,'LineStyle',':','Marker','none');
hold(axel1,'all');

% Show equation and R2 in the upper right corner
legtext = sprintf('%s','y = ',num2str(slope,3),'x, R^2 = ',num2str(R2,2));
text(axel1,axel1.XLim(2)*0.98,axel1.YLim(2)*0.95, legtext,...
    'Vert','bottom', 'Horiz','right', 'FontSize',20);

% Save as emf file
thefile = sprintf('%s',path,'ScatterCalibration');
saveas(figure1,thefile,'emf');
delete(figure1);
%%%


%%% MAKE MAPS COVERING EUROPE
figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 11.0 9.8]);
worldmap([35 70],[-10 35]);

% Uncallibrated modelled damage in severe storm (kg C m-2)
theshow = surfm(lats,lons,stormmap(:,:,3));
setm(gca, 'FontSize', 20);
setm(gca,'GLineWidth',1.0);
borders('countries');
colorbar('FontSize',20);
thefile = sprintf('%s',path,'DamWood_uncall');
saveas(figure1,thefile,'emf');

% Callibrated modelled damage in severe storm (kg C m-2)
set(theshow,'CData',stormmap_call(:,:,3));
thefile = sprintf('%s',path,'DamWood_call');
saveas(figure1,thefile,'emf');

% Wind damage probability map (yr-1)
set(theshow,'CData',wdpmap);
thefile = sprintf('%s',path,'wdp');
saveas(figure1,thefile,'emf');

% Uncallibrated expected modelled damage (kg C m-2 yr-1)
set(theshow,'CData',damagemap_uncall);
thefile = sprintf('%s',path,'Damage_uncall');
saveas(figure1,thefile,'emf');

% Callibrated expected modelled damage (kg C m-2 yr-1)
set(theshow,'CData',damagemap_call);
thefile = sprintf('%s',path,'Damage_call');
saveas(figure1,thefile,'emf');
%%%
