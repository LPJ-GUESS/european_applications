% Make plots of new allometries and how they compare to the previous LPJ-GUESS allometries
%
% T. Pugh
% 30.07.24

pftname={'Abi alb','Bet pen','Bet pub','Car bet','Cor ave','Fag syl','Fra exc','Jun oxy',...
    'Lar dec','Pic abi','Pic sit','Pin syl','Pin hal','Pop tre','Que coc','Que ile',...
    'Que pub','Que rob','Til cor','Ulm gla'};
pftclass={'needle_cold','broad_cold','broad_cold','broad_cold','broad_cold','broad_cold','broad_warm','needle_warm',...
    'needle_cold','needle_cold','needle_cold','needle_cold','needle_warm','broad_warm','broad_warm','broad_warm',...
    'broad_warm','broad_cold','broad_cold','broad_cold'};
tol={'tol','intol','intol','inter','inter','tol','inter','intol',...
    'inter','tol','inter','inter','intol','intol','inter','inter',...
    'inter','inter','inter','inter'};
sla=[14.85,32.6,32.6,39.9,34.2,57.0,29.4,13.0,...
    13.0,12.6,12.6,10.9,11.2,29.1,25.5,13.3,...
    24.9,22.8,45.4,31.3];
k_allom2=[54.2,58.2,58.2,34.3,29.3,34.2,58.7,29.3,...
    51.0,61.2,40.1,30.7,20.6,49.2,36.1,9.6,...
    15.4,24.6,44.2,41.1];
k_allom3=[0.796,0.675,0.675,0.476,0.461,0.445,0.684,0.461,...
    0.5,0.818,0.743,0.582,0.650,0.623,0.548,0.333,...
    0.387,0.503,0.575,0.643];
k_allom1=[61.0,193.4,193.4,213.8,109.1,157.4,104.8,104.8,...
    61.0,59.4,130.3,118.6,130.8,106.6,191.8,196.0,...
    144.3,176.4,30.3,207.6];
k_rp=[1.109,1.430,1.430,1.277,1.0,1.242,1.358,1.358,...
    1.109,1.103,1.300,1.552,1.407,1.235,1.054,1.617,...
    1.480,1.466,1.0,1.253];
wooddens=[252.5,285.0,285.0,340.0,292.5,337.5,332.5,287.5,...
    262.5,225.0,182.5,282.5,307.5,238.8,257.5,420.0,...
    377.5,307.5,227.5,312.5];

npft=length(pftname);

% Standard allometries
k_allom2_stan=40;
k_allom3_stan=0.67;
k_rp_stan=1.6;
k_allom1_stan=NaN(size(k_allom1));
for nn=1:npft
    if strcmp(pftclass(nn),'needle_cold') || strcmp(pftclass(nn),'needle_warm')
        k_allom1_stan(nn)=150;
    elseif strcmp(pftclass(nn),'broad_cold') || strcmp(pftclass(nn),'broad_warm')
        k_allom1_stan(nn)=250;
    end
end
clear nn

% Calculate dimension ranges for each PFT
dbh_range=5:1:100; % in cm
dbh_range_m=dbh_range/100; % in m

pft_height=NaN(npft,length(dbh_range));
pft_height_stan=NaN(npft,length(dbh_range));
pft_CA=NaN(npft,length(dbh_range));
pft_CA_stan=NaN(npft,length(dbh_range));
for nn=1:npft
    pft_height(nn,:)=k_allom2(nn)*(dbh_range_m.^k_allom3(nn));
    pft_height_stan(nn,:)=k_allom2_stan*(dbh_range_m.^k_allom3_stan);
    pft_CA(nn,:)=k_allom1(nn)*(dbh_range_m.^k_rp(nn));
    pft_CA_stan(nn,:)=k_allom1_stan(nn)*(dbh_range_m.^k_rp_stan);
end
clear nn

% Make figures
figure
for nn=1:npft
    subplot(5,4,nn)
    hold on
    plot(dbh_range,pft_height_stan(nn,:))
    plot(dbh_range,pft_height(nn,:))
    title(pftname{nn})
end
legend('v4.1','Updated')

figure
for nn=1:npft
    subplot(5,4,nn)
    hold on
    plot(dbh_range,pft_CA_stan(nn,:))
    plot(dbh_range,pft_CA(nn,:))
    title(pftname{nn})
end
legend('v4.1','Updated')