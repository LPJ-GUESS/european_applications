
/// Version that uses files without added PFT column, sorted by plot id and names in all header columns (line_no) (fileName_sort); new file is created with PNV column (_addFT). A species_PFT map file (fileName_species_map) is required.

//************************************************************ HAOMING MODIFY *****************************************
// Define this first time you use a plot file, then this can be commented out for faster run.
#define USE_FILE_WITHOUT_PFT
//*********************************************************************************************************************

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>

#define MAXLINE 20000
#define MAXRECORDS 100
#define MAXNAMESIZE 20
#define NSPECIESTOT 8950
#define NSPECIESTOTACTUAL 2670
#define MAXNSPECIESINDATA 500

#define FIRSTPRINTYEAR 1500		// 1268			; Ages will be capped at FIRSTPRINTYEAR+1
#define LASTPRINTYEAR 2008

#define NPLOTSMAX 8800
#define NCENSUSMAX 15800

double const PI = 4 * atan(1.0);

double basal_area(double diam) {

	return PI * pow(0.5 * diam, 2);
}

int getspeciespos(int specieslist[], int speciesid) {

	int pos = -1;

	for(unsigned int i=0; i< MAXNSPECIESINDATA; i++) {

		if(specieslist[i] == speciesid) {
			pos = i;
			break;
		}
	}

	return pos;
}

struct PlotData {

	char plotid[MAXNAMESIZE];
	double latitude;
	double longitude;
	double area;
	int elevation;

	PlotData() {
		for(int i=0; i<MAXNAMESIZE; i++) {
			plotid[i] = 0;
		}
		elevation = 0;
		latitude = 0;
		longitude = 0;
		area = 0;
	}
};

int getcoordpos_from_struct(PlotData plotlist[], char* plotidX) {

	int pos = -1;

	for(unsigned int i=0; i< NPLOTSMAX; i++) {

		if(!strcmp(plotlist[i].plotid, plotidX) && strcmp(plotlist[i].plotid, "")) {
			pos = i;
			break;
		}
	}
	return pos;
}

struct CensusData {

	char plotid[MAXNAMESIZE];
	int census_year;
	int census_no;
	int stand_age;

	CensusData() {
		for(int i=0; i<MAXNAMESIZE; i++) {
			plotid[i] = 0;
		}
		census_year = -1;
		stand_age = -1;
	}
};

int getagepos_from_struct(CensusData censuslist[], char* plotidX, int year) {

	int pos = -1;

	for(unsigned int i=0; i< NCENSUSMAX; i++) {

		if(!strcmp(censuslist[i].plotid, plotidX) && censuslist[i].census_year == year) {
			pos = i;
			break;
		}
	}
	return pos;
}

struct tree_struct {

	char pft_name[30];
	int census_no;
	int census_year;
	int age;
	int dead;
	double diam;
	double height;
	double agb;
	double ba;
	double dens;
//	double cmass_leaf_est;
//	double cmass_root_est;
};


const int MAX_NTREES_PER_PLOT = 500;
const int MAX_NESTABLISH_INDIV = 100;

struct est_struct {

	char pft_name[30];
	int census_no;
	int census_year;
	int age;
	double diam;
	double height;
	double agb;
	double ba;
	double dens;
//	double cmass_leaf_est;
//	double cmass_root_est;

};

class Establishment_this_year {

public:
	tree_struct tree_array[MAX_NTREES_PER_PLOT];
	tree_struct est_array[MAX_NESTABLISH_INDIV];

	int ntrees;

	Establishment_this_year() {
		zero();
	}

	void zero() {

		ntrees = 0;

		for(int i=0;i<MAX_NTREES_PER_PLOT;i++) {
			tree_array[i].pft_name[0] = '\0';
			tree_array[i].census_no = 0;
			tree_array[i].census_year = 0;
			tree_array[i].age = 0;
			tree_array[i].diam = 0.0;
			tree_array[i].agb = 0.0;
			tree_array[i].ba = 0.0;
			tree_array[i].height = 0.0;
			tree_array[i].dens = 0.0;
			tree_array[i].dead = 0;
		}
		for(int i=0;i<MAX_NESTABLISH_INDIV;i++) {
			est_array[i].pft_name[0] = '\0';
			est_array[i].census_no = 0;
			est_array[i].census_year = 0;
			est_array[i].age = 0;
			est_array[i].diam = 0.0;
			est_array[i].agb = 0.0;
			est_array[i].ba = 0.0;
			est_array[i].height = 0;
			est_array[i].dens = 0.0;
		}
	}

	void add_tree(int tree_no, char* pft_name, int census_no, int census_year, double diameter, double height, double agb, double ba, double dens) {

		// add data to tree list
		strcpy(tree_array[ntrees].pft_name, pft_name);
		tree_array[ntrees].census_no = census_no;
		tree_array[ntrees].census_year = census_year;
		tree_array[ntrees].diam = diameter;
		tree_array[ntrees].height = height;
		tree_array[ntrees].agb = agb;
		tree_array[ntrees].ba = ba;
		tree_array[ntrees].dens = dens;

		ntrees++;
	}
};

int main(int argc,char* argv[]) {

//	printf("HHH\n");
	FILE *ifp=NULL, *ifp_cens=NULL, *ofp=NULL, *ofp_sp=NULL, *ofp_lu=NULL, *ofp_st=NULL, *ofp_firstmanageyear=NULL, *ofp_init=NULL, *ofp_log=NULL;
	FILE *ifp_sort = NULL, *ifp_sp_map = NULL, *ofp_addPFT = NULL;
	int size=0, count1=0, nnames=0;
	bool species_used[NSPECIESTOT] = {0};
	int species_ids[MAXNSPECIESINDATA] = {0};
	double ba[MAXNSPECIESINDATA] = {0};
	char line[MAXLINE]={0}, line2[MAXLINE]={0},ext[]=".txt", tag[]="_calc", tag_init[]="_init", s1[MAXRECORDS][MAXNAMESIZE]={'\0'};
	char *outFileNameAddPFTnames=NULL,*outFileName=NULL,*chp=NULL,*outFileNameInit=NULL,*outFileNameLog=NULL, *p=NULL;
	bool firstline=true;

//	typedef enum {BNE, BINE, BNS, TeBS, IBS, TeBE, NPFTS} pfttype;
//	const char* PFTNAMES[] = {"BNE", "BINE", "BNS", "TeBS", "IBS", "TeBE"};
	typedef enum {Abi_alb, Bet_pen, Bet_pub, Car_bet, Cor_ave, Fag_syl, Fra_exc, Lar_dec, Pic_abi, Pin_syl, Pin_hal, Pop_tre, Que_coc, Que_ile, Que_pub, Que_rob, Til_cor, Ulm_gla, NPFTS} pfttype;
	const char* PFTNAMES[] = {"Abi_alb", "Bet_pen", "Bet_pub", "Car_bet", "Cor_ave", "Fag_syl", "Fra_exc", "Lar_dec", "Pic_abi", "Pin_syl", "Pin_hal", "Pop_tre", "Que_coc", "Que_ile", "Que_pub", "Que_rob", "Til_cor", "Ulm_gla"};

	PlotData plot_data[NPLOTSMAX];
	CensusData census_data[1];

///***///	File names specified here

	char fileName_species[] = "C:\\pugh_entire_US.csv\\US_REF_SPECIES.csv";

#ifndef USE_FILE_WITHOUT_PFT
// This option is to use input file with PFT already in one column, this file is created the first run without this option, it's called, xxx_addPFT.txt
//************************************************************ HAOMING MODIFY *****************************************
//	char fileName[] = "C:\\Forest_init_C\\Windows text format\\NFI_SwedenGermany3_Biomass_addPFT.txt";
	char fileName[] = "C:\\Forest_init_C\\\\NFI_SwedenGermany3_Biomass_CRLF_noNGE_addPFT.txt";	// Input files sorted by plot name
//*********************************************************************************************************************
#else
// This option needed is needed the first run, to create input file with PFT in one column, called, xxx_addPFT.txt. After this, it's faster to use the other option
//************************************************************ HAOMING MODIFY *****************************************
	char fileName_species_map[] = "E:\\Forest_init\\Haoming\\Species Mapping_Final Result_NoShrub_fixed.txt";	// Input file with species and PFT columns
//	char fileName_sort[] = "C:\\Forest_init_C\\Windows text format\\NFI_SwedenGermany3_Biomass.txt";	// Input files sorted by plot name
	char fileName_sort[] = "E:\\Forest_init\\NFI_SwedenGermany3_Biomass_CRLF_noNGE.txt\\NFI_SwedenGermany3_Biomass_CRLF_noNGE.txt";	// Input files sorted by plot name
//*********************************************************************************************************************	
	char *fileName = NULL;	// Temporary input file with PFT column added
	char tag_addPFT[] = "_addPFT";

// Add PFT column to tree input file
	ifp_sort = fopen(fileName_sort, "r");
	if(!ifp_sort) {
		printf("File %s could not be opened for input !\n\n", fileName_sort);
		return 1;
	}
	ifp_sp_map = fopen(fileName_species_map, "r");
	if(!ifp_sp_map) {
		printf("File %s could not be opened for input !\n\n", fileName_species_map);
		return 1;
	}

	if(chp = strrchr(fileName_sort, '.')) {
		*chp = '\0';
	}

	size=strlen(fileName_sort) + strlen(tag_addPFT) + strlen(ext)+1;
	outFileNameAddPFTnames = new char[size];
	memset(outFileNameAddPFTnames, 0, size); 
	strncpy(outFileNameAddPFTnames, fileName_sort, strlen(fileName_sort));
	strncat(outFileNameAddPFTnames, tag_addPFT, strlen(tag_addPFT));;
	strncat(outFileNameAddPFTnames, ext, strlen(ext));
	ofp_addPFT = fopen(outFileNameAddPFTnames, "w");
	if(!ofp_addPFT) {
		printf("File %s could not be opened for output !\n\n", outFileNameAddPFTnames);
		if(outFileNameAddPFTnames)
			delete[] outFileNameAddPFTnames;
		return 1;
	}


	int sp_map_species_col, sp_map_PFT_col;

	if(fgets(line,sizeof(line),ifp_sp_map)) {

		p=strtok(line, "\t\n,\"");	
		if(p) {
			strncpy(s1[count1], p, MAXNAMESIZE-1);
		}
		count1++;
		do {
			p=strtok(NULL, "\t\n,\"");
			if(p) {
				strncpy(s1[count1], p, MAXNAMESIZE-1);
				count1++;
			}
		}
		while(p);

		p=NULL;

		for(int i=0;i<count1;i++) {

			if(!strcmp(s1[i], "species.cor"))
				sp_map_species_col = i;
			else if(!strcmp(s1[i], "PFT"))
				sp_map_PFT_col = i;
		}
	}

	rewind(ifp_sp_map);


	int sort_species_col;
	firstline = true;
	int line_no_sort = 0;
	char raw_line[20000] = {0};

	while(!feof(ifp_sort)) {

		count1=0;	

		if(fgets(raw_line,sizeof(raw_line),ifp_sort)) {

			char species_name[200] = {0};
			char PFT_name[200] = {0};

			strcpy(line, raw_line);	
			p=strtok(line, "\t\n,\"");					// Names contain spaces

			if(p) {
				strncpy(s1[count1], p, MAXNAMESIZE-1);
			}
			count1++;
			do {
				p=strtok(NULL, "\t\n,\"");
				if(p) {
					strncpy(s1[count1], p, MAXNAMESIZE-1);
					count1++;
				}
			}
			while(p);

			p=NULL;

			//Header row:
			if(firstline) {

				for(int i=0;i<count1;i++) {

					if(!strcmp(s1[i], "species.cor"))
						sort_species_col = i;
				}
				firstline=false;

				if(chp = strrchr(raw_line, '\n'))
					*chp = '\0';
				strcat(raw_line, "\t");
				strcat(raw_line, "PFT\n");
				fprintf(ofp_addPFT, raw_line);
			}
			else {

				if(sort_species_col != -1)
					strcpy(species_name, s1[sort_species_col]);

				char sp_map_line[20000] = {0}, s1_sp_map[MAXRECORDS][MAXNAMESIZE]={'\0'};;
				bool firstline_sp_map = true;

				while(!feof(ifp_sp_map)) {

					count1=0;	

					if(fgets(sp_map_line,sizeof(sp_map_line), ifp_sp_map)) {

						p=strtok(sp_map_line, "\t\n,\"");					// Names contain spaces

						if(p) {
							strncpy(s1_sp_map[count1], p, MAXNAMESIZE-1);
						}
						count1++;
						do {
							p = strtok(NULL, "\t\n,\"");
							if(p) {
								strncpy(s1_sp_map[count1], p, MAXNAMESIZE-1);
								count1++;
							}
						}
						while(p);

						p = NULL;

						//Header row:
						if(firstline_sp_map) {
							firstline_sp_map = false;
						}
						else {

							if(!strcmp(s1_sp_map[sp_map_species_col], species_name)) {

								if(sp_map_PFT_col != -1)
									strcpy(PFT_name, s1_sp_map[sp_map_PFT_col]);
								rewind(ifp_sp_map);
								break;
							}
						}
					}
				}

				if(chp = strrchr(raw_line, '\n'))
					*chp = '\0';
				strcat(raw_line, "\t");
				strcat(raw_line, PFT_name);
				strcat(raw_line, "\n");
				fprintf(ofp_addPFT, raw_line);
				if(!strcmp(PFT_name, "")) {
					printf("Tree species name %s not found in species PFT map file\n", species_name);
					return 0;
				}

				line_no_sort++;
			}
		}
	}

	if(ifp_sort)
		fclose(ifp_sort);

	if(ofp_addPFT)
		fclose(ofp_addPFT);

	size = strlen(outFileNameAddPFTnames) + 1;
	fileName = new char[size];
	memset(fileName, 0, size); 
	strncpy(fileName, outFileNameAddPFTnames, strlen(outFileNameAddPFTnames));

////////////////////////////////////////////
#endif

	ifp=fopen(fileName, "r");
	if(!ifp) {
		printf("File %s could not be opened for input !\n\n", fileName);
		return 1;
	}

	if(chp=strrchr(fileName, '.')) {
		*chp='\0';
	}

	size=strlen(fileName)+strlen(tag)+strlen(ext)+1;
	outFileName=new char[size];
	memset(outFileName, 0, size); 
	strncpy(outFileName, fileName, strlen(fileName));
	strncat(outFileName, tag, strlen(tag));
	strncat(outFileName, ext, strlen(ext));
	ofp=fopen(outFileName, "w");
	if(!ofp) {
		printf("File %s could not be opened for output !\n\n", outFileName);
		if(outFileName)
			delete[] outFileName;
		return 1;
	}

	fprintf(ofp, "%s\t%s\t%s\t%s\t%s\t%s", "Lon", "Lat", "No", "Year", "Age", "Elevation");
	for(int pft=0;pft<NPFTS;pft++) {
		fprintf(ofp, "\t%s", PFTNAMES[pft]);
	}
	fprintf(ofp, "\n");

//////////////////
	size=strlen(fileName)+strlen(tag_init)+strlen(ext)+1;
	outFileNameInit=new char[size];
	memset(outFileNameInit, 0, size); 
	strncpy(outFileNameInit, fileName, strlen(fileName));
	strncat(outFileNameInit, tag_init, strlen(tag_init));
	strncat(outFileNameInit, ext, strlen(ext));
	ofp_init=fopen(outFileNameInit, "w");
	if(!ofp_init) {
		printf("File %s could not be opened for output !\n\n", outFileNameInit);

		if(outFileNameInit)
			delete[] outFileNameInit;
		return 1;
	}

	size=strlen(fileName)+strlen(tag_init)+strlen(ext)+1;
	outFileNameLog=new char[size];
	memset(outFileNameLog, 0, size); 
	strncpy(outFileNameLog, fileName, strlen(fileName));
	strcat(outFileNameLog, "_log");
	strncat(outFileNameLog, ext, strlen(ext));
	ofp_log=fopen(outFileNameLog, "w");
	if(!ofp_log) {
		printf("File %s could not be opened for output !\n\n", outFileNameLog);

		if(outFileNameLog)
			delete[] outFileNameLog;
		return 1;
	}

	char last_stateid[100] ={0};
	char last_countyid[100] ={0};
	char last_plotid[100] ={0};
	double last_lon = 200;
	double last_lat = 200;

	int last_year = 0;

///////////////////////////////
/////// MAIN TREE INPUT FILE///
///////////////////////////////
	int pftid_col = -1, treeid_col = -1, diam_col = -1, ba_col = -1, agb_col = -1, plotid_col = -1, age_col = -1, year_col = -1, census_no_col = -1, treedens_col = -1;
	int height_col = -1, bole_height_col = -1, crownarea_col = -1, lai_col = -1, tree_age_col = -1, dead_col = -1;
	int longitude_col = -1, latitude_col = -1, area_inv_col = -1;

///***///	// Read tree file

	double ba_plot = 0.0;
	double ba_pft[NPFTS] = {0};
	double agb_plot = 0.0;
	double agb_pft[NPFTS] = {0};

//************************************************************ HAOMING POSSIBLY MODIFY *****************************************
	const int ndiam_bins = 15;
	const int diam_bin_min = 0;	// Minimum size 0 cm
	const int diam_bin_size = 100;	// 10-cm bins (in mm !)

	const bool use_last_census = true;
	const bool allow_dead_census = true;					// Use first or last census, depending on use_last_census, even if the trees are all dead (no trees established) DEFAULT true
	const bool use_census_closest_to_default_year = true;	// Overrides use_last_census value
	const bool replace_empty_census_year_with_default_year = false;	//  DEFAULT false
	const int default_year = 2010;					// Set this to 0 if trees (plots) with(out?) census year is to be ignored (replace_empty_census_year_with_default_year used instead).
//******************************************************************************************************************************
	if(use_census_closest_to_default_year && !default_year) {
		printf("\ndefault_year must be set with use_census_closest_to_default_year option\n");
		return 1;
	}

	Establishment_this_year establish_list;

	int max_ntrees_per_plot_census = 0;		// Could possibly change between laps if trees are disallowed by some rule in lap 2 (lap 2 values are printed)
	int max_ntrees_per_plot_total = 0;		// Could possibly change between laps if trees are disallowed by some rule in lap 2 (lap 2 values are printed)
	int max_ncohorts = 0;					// Could possibly change between laps if trees are disallowed by some rule in lap 2 (lap 2 values are printed)
	int max_ncohorts_save = 0;				// This is the value from lap 2, which is used for sizing the umber of output lines
	int max_ncensuses = 0;					// Could possibly change between laps if trees are disallowed by some rule in lap 2 (lap 2 values are printed)

	const int census_no_max = 25;
	int census_nos[census_no_max] = {0};

	int line_no = 0;
	int lastline_no = 0;

	int nplots = 0;


	const int nlaps = 2;
	for (int lap_no=0;lap_no<nlaps;lap_no++) {

		lastline_no = line_no - 1;	// Used in second lap
		max_ncohorts_save = max_ncohorts;

		max_ntrees_per_plot_census = 0;
		max_ntrees_per_plot_total = 0;
		max_ncohorts = 0;
		max_ncensuses = 0;	// Only used in lap 0, so not really needed to xero here

		firstline = true;
		line_no = 0;
		int tree_no = 0;
		bool second_census_last_time = false;
		int use_census = 0;
		int use_census_year = 0;
		int ncensuses2 = 0;

		double smallest_dif = 1000;

		while(!feof(ifp)) {

			count1=0;

			if(fgets(line,sizeof(line),ifp)) {

/////////////// Tokenise string (put /0 between "words", copy to string array s1)
				p=strtok(line, "\t\n,\"");		// Names contain spaces

				if(p) {
					strncpy(s1[count1], p, MAXNAMESIZE-1);
				}
				count1++;
				do {
					p=strtok(NULL, "\t\n,\"");		// Names contain spaces
					if(p) {
						strncpy(s1[count1], p, MAXNAMESIZE-1);
						count1++;
					}
				}
				while(p);

				p=NULL;
/////////////// 

	//Header row:
				if(firstline) {

					for(int i=0;i<count1;i++) {

						if(!strcmp(s1[i], "PFT"))
							pftid_col = i;
						else if(!strcmp(s1[i], "d"))
							diam_col = i;
						else if(!strcmp(s1[i], "tree.id"))
							treeid_col = i;
						else if(!strcmp(s1[i], "ba"))	// mm2
							ba_col = i;
//						else if(!strcmp(s1[i], "agb"))	// ton/tree
						else if(!strcmp(s1[i], "biomass"))	// Updated 220324, again 240912
							agb_col = i;
						else if(!strcmp(s1[i], "height"))
							height_col = i;
						else if(!strcmp(s1[i], "bole_height"))	// NY
							bole_height_col = i;
						else if(!strcmp(s1[i], "crownarea"))	// NY
							crownarea_col = i;
						else if(!strcmp(s1[i], "lai"))			// NY
							lai_col = i;
						else if(!strcmp(s1[i], "tree_age"))		// NY
							tree_age_col = i;
//						else if(!strcmp(s1[i], "plot.id"))
						else if(!strcmp(s1[i], "tmt.plot.id"))	// NB. WE NOW USE THE UNIQUE tmt.plot.id, NOT plot.id !!!!!
							plotid_col = i;
						else if(!strcmp(s1[i], "STDAGE"))
							age_col = i;
						else if(!strcmp(s1[i], "census.date"))
							year_col = i;
						else if(!strcmp(s1[i], "census.n"))
							census_no_col = i;
						else if(!strcmp(s1[i], "n.ha"))
							treedens_col = i;
						else if(!strcmp(s1[i], "tree.status"))
							dead_col = i;
						else if(!strcmp(s1[i], "plot.lon"))
							longitude_col = i;
						else if(!strcmp(s1[i], "plot.lat"))
							latitude_col = i;
						else if(!strcmp(s1[i], "n.ha"))			// Finns redan som dens !
							area_inv_col = i;
					}
					firstline=false;

					if(lap_no) {
							fprintf(ofp_init, "%s\t%s\t%s\t%s\t", "Lon", "Lat", "tree_no", "str_pft");
							if(diam_col != -1)
								fprintf(ofp_init, "%s\t", "dbh");
							if(height_col != -1)
								fprintf(ofp_init, "%s\t", "height");
//							if(treedens_col != -1)
							fprintf(ofp_init, "%s\t", "dens");
							if(agb_col != -1)
//								fprintf(ofp_init, "%s\t", "agb");
								fprintf(ofp_init, "%s\t", "biomass");	// 240912
							if(ba_col != -1)
								fprintf(ofp_init, "%s\t", "ba");
							if(bole_height_col != -1)
								fprintf(ofp_init, "%s\t", "bole_height");
							if(crownarea_col != -1)
								fprintf(ofp_init, "%s\t", "crownarea");
							if(lai_col != -1)
								fprintf(ofp_init, "%s\t", "lai");
							if(tree_age_col != -1)
								fprintf(ofp_init, "%s\t", "tree_age");
							fprintf(ofp_init, "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "firstyear", "plotid", "plot_area", "census_id", "census_year", "init_year", "ntrees");
					}
				}
				else {	// if(firstline)

					int census_no = -1;
					int census_year = -1;
					int dead = -1;
					double diameter = -1;
					double height = -1;
					double agb = -1;
					double ba = -1;
					double dens = -1;

					if(chp = strrchr(s1[year_col], '.'))		// remove trailing decimals from census year (not needed for strtol() to work.
						*chp='\0';

					if(census_no_col != -1)
						census_no = strtol(s1[census_no_col], NULL, 0);
					if(year_col != -1)
						census_year = strtol(s1[year_col], NULL, 0);	// returns 0 if "NA"
					if(!census_year && replace_empty_census_year_with_default_year && !use_census_closest_to_default_year)
						census_year = default_year;						// Set default year for trees with census year = NA.(don't use it when using closest to census year)
																		// Set this to 0 if trees (plots) with(out?) census year are to be ignored.
					if(dead_col != -1)
						dead = strtol(s1[dead_col], NULL, 0);	
					if(diam_col != -1)
						diameter = strtod(s1[diam_col], NULL);
					if(height_col != -1)
						height = strtod(s1[height_col], NULL);		// returns 0 if "NA"
					if(agb_col != -1)
						agb = strtod(s1[agb_col], NULL);
					if(ba_col != -1)
						ba = strtod(s1[ba_col], NULL);
					if(treedens_col != -1)
						dens = strtod(s1[treedens_col], NULL);	

					if(line_no != 1) {
						tree_no++;
					}
					if(tree_no > MAX_NTREES_PER_PLOT) {	// More trees (all censuses) than room in tree object
						printf("Too many trees in plot %s: Increase MAX_NTREES_PER_PLOT !\n", s1[plotid_col]);
						fprintf(ofp_log, "Too many trees in plot %s: Increase MAX_NTREES_PER_PLOT !\n",s1[plotid_col]);
					}

					if(!lap_no) {	// This lap is just to obtain max_ncohorts_save and lastline, to minimise output file size
								// (prints max_ncohorts, max_ncensuses, max_ntrees_per_plot_census and max_ntrees_per_plot_total)

						// Same plot as last line
						if(!strcmp(s1[plotid_col], last_plotid)) {

//							if(census_year)	{	// Ignore lines with no census year !
							if(census_year && (!dead || allow_dead_census))	{	// Ignore dead trees !
								census_nos[census_no] = 1;
								if(census_no != use_census) {

									if(use_census_closest_to_default_year) {

										double dif_to_default_year = fabs(double(census_year - default_year));
										if(dif_to_default_year < smallest_dif) {
											smallest_dif = dif_to_default_year;
											use_census = census_no;
											use_census_year = census_year;
										}
									}
									else if(!use_census || (census_no > use_census && use_last_census) || (census_no < use_census && !use_last_census)) {
										use_census = census_no;
									}
								}
								// add data to tree list
								establish_list.add_tree(tree_no, s1[pftid_col], census_no, census_year, diameter, height, agb, ba, dens);
							}
						}
						else {	// New plot line
							if(line_no != 1) {

								int ncohorts = 0;

								for(int pft=0;pft<NPFTS;pft++) {
//									for(int db=0; db<ndiam_bins; db++) {
									for(int db=0; db<ndiam_bins + 1; db++) {

										int ntrees_cohort = 0;
										for(int i=0; i<MAX_NTREES_PER_PLOT;i++) {
											if(establish_list.tree_array[i].census_no == use_census) {
												if(!strcmp(establish_list.tree_array[i].pft_name, PFTNAMES[pft])) {
													double diam_indiv = establish_list.tree_array[i].diam;
													if(diam_indiv != 0.0 && diam_indiv >= diam_bin_min + (diam_bin_size * db)
//														&& diam_indiv < diam_bin_min + (diam_bin_size * (db +1))) {
														&& (diam_indiv < diam_bin_min + (diam_bin_size * (db + 1)) || (db == ndiam_bins))) {

														ntrees_cohort++;
													}
												}
											}
										}
										if(ntrees_cohort) {

											if((ncohorts + 1) > MAX_NESTABLISH_INDIV) {
												printf("Too many cohorts in plot %s: Increase MAX_NESTABLISH_INDIV !\n", last_plotid);
												fprintf(ofp_log, "Too many cohorts in plot %s: Increase MAX_NESTABLISH_INDIV !\n", last_plotid);
											}
											ncohorts++;

											if(ntrees_cohort > max_ntrees_per_plot_census)
												max_ntrees_per_plot_census = ntrees_cohort;		// Only used for printed info
										}
									}
								}

								if(tree_no > max_ntrees_per_plot_total)
									max_ntrees_per_plot_total = tree_no;						// Only used for printed info
								if(ncohorts > max_ncohorts)
									max_ncohorts = ncohorts;									// Used to size output line number

								for(int c=0;c<census_no_max;c++) {
									if(census_nos[c])
										ncensuses2++;
								}
								if(ncensuses2 > max_ncensuses)
									max_ncensuses = ncensuses2;									// Only used for printed info

								establish_list.zero();
								tree_no = 0;
								use_census = 0;
								ncensuses2 = 0;
								for(int c=0;c<census_no_max;c++) {
									census_nos[c] = 0;
//									trees_per_census[c] = 0;
								}
								smallest_dif = 1000;
							}

//							if(census_year)	{	// Ignore lines with no census year !
							if(census_year && (!dead || allow_dead_census))	{	// Ignore dead trees !
								census_nos[census_no] = 1;
								if(census_no != use_census) {
									if(use_census_closest_to_default_year) {

										double dif_to_default_year = fabs(double(census_year - default_year));
										if(dif_to_default_year < smallest_dif) {
											smallest_dif = dif_to_default_year;
											use_census = census_no;
											use_census_year = census_year;
										}
									}
									else if(!use_census || (census_no > use_census && use_last_census) || (census_no < use_census && !use_last_census)) {
										use_census = census_no;
									}
								}

								// add data to tree list for first tree in plot
								establish_list.add_tree(tree_no, s1[pftid_col], census_no, census_year, diameter, height, agb, ba, dens);
							}
							nplots++;
						}
					}	// if(!r)
					else {			// Write lap

						// Same plot as last line; add data to tree list
						if(!strcmp(s1[plotid_col], last_plotid)) {
//							if(census_year)	{	// Ignore lines with no census year !
							if(census_year && (!dead || allow_dead_census))	{	// Ignore dead trees !
								census_nos[census_no] = 1;
								if(census_no != use_census) {
									if(use_census_closest_to_default_year) {

										double dif_to_default_year = fabs(double(census_year - default_year));
										if(dif_to_default_year < smallest_dif) {
											smallest_dif = dif_to_default_year;
											use_census = census_no;
											use_census_year = census_year;
										}
									}
									else if(!use_census || (census_no > use_census && use_last_census) || (census_no < use_census && !use_last_census)) {
										use_census = census_no;
										use_census_year = census_year;
									}
								}

								establish_list.add_tree(tree_no, s1[pftid_col], census_no, census_year, diameter, height, agb, ba, dens);
							}
						}
						// New plot
						if(strcmp(s1[plotid_col], last_plotid) || line_no == lastline_no) {
		// Last plot will not be printed if only one census !

		///***///	// DO things that need to be done before moving on to new plot (transfer mean values to dimaeter cohorts, print to file)

							// Add data to establishment list, then zero tree list and add new line to tree list
							if(line_no == 1) {
								establish_list.zero();
							}
							else {

								double plot_area;

								if(dens)
									plot_area = 1.0 / dens;

								// Determine which census to use

								int use_census2 = 0;;
								int use_census_year2 = 0;
								int census_nos2[census_no_max] = {0};

								double smallest_dif = 1000;

								for(int i=0; i<establish_list.ntrees;i++) {

									int census_no_local = establish_list.tree_array[i].census_no;
									int census_year_local = establish_list.tree_array[i].census_year;	
									int dead_local = establish_list.tree_array[i].dead;

									if(census_year_local && (!dead_local || allow_dead_census))	{	// Ignore dead trees !
										census_nos2[census_no] = 1;

										if(use_census_closest_to_default_year) {

											double dif_to_default_year = fabs(double(establish_list.tree_array[i].census_year - default_year));
											if(dif_to_default_year < smallest_dif) {
												smallest_dif = dif_to_default_year;
												use_census2 = census_no_local;
												use_census_year2 = census_year_local;
											}

										}
										else if(!use_census2 || (census_no_local > use_census2 && use_last_census) || (census_no_local < use_census2 && !use_last_census)) {
											use_census2 = census_no_local;
											use_census_year2 = census_year_local;
										}
									}
								}

								use_census = use_census2;
								use_census_year = use_census_year2;

								int ncohorts = 0;
								for(int pft=0;pft<NPFTS;pft++) {
//									for(int db=0; db<ndiam_bins; db++) {
									for(int db=0; db<ndiam_bins + 1; db++) {

										// NB: trees with diameter 0 are ignored here !
										double diam_cohort = 0.0;
										double height_cohort = 0.0;
										double agb_cohort = 0.0;
										double ba_cohort = 0.0;
										double dens_cohort = 0.0;
										int ntrees_cohort = 0;
										for(int i=0; i<MAX_NTREES_PER_PLOT;i++) {
											if(establish_list.tree_array[i].census_no == use_census) {
												if(!strcmp(establish_list.tree_array[i].pft_name, PFTNAMES[pft])) {	// Går inte in här !!!
													double diam_indiv = establish_list.tree_array[i].diam;
													if(diam_indiv != 0.0 && diam_indiv >= diam_bin_min + (diam_bin_size * db)
//														&& diam_indiv < diam_bin_min + (diam_bin_size * (db +1))) {
														&& (diam_indiv < diam_bin_min + (diam_bin_size * (db + 1)) || (db == ndiam_bins))) {

														diam_cohort += diam_indiv;
														height_cohort += establish_list.tree_array[i].height;
														agb_cohort += establish_list.tree_array[i].agb;
														ba_cohort += establish_list.tree_array[i].ba;
														dens_cohort += establish_list.tree_array[i].dens;
														ntrees_cohort++;
													}
												}
											}
										}

										if(ntrees_cohort) {

											if((ncohorts + 1) > MAX_NESTABLISH_INDIV) {
												printf("Too many cohorts in plot %s: Increase MAX_NESTABLISH_INDIV !\n", last_plotid);
												fprintf(ofp_log, "Too many cohorts in plot %s: Increase MAX_NESTABLISH_INDIV !\n", last_plotid);
											}

											// Print to file !

											diam_cohort /= ntrees_cohort;
											height_cohort /= ntrees_cohort;
											agb_cohort /= ntrees_cohort;	// 240912

											fprintf(ofp_init, "%.9f\t%.9f\t", last_lon, last_lat);
											fprintf(ofp_init, "%d\t", ncohorts);
											fprintf(ofp_init, "%s\t", PFTNAMES[pft]);
		//									fprintf(ofp_init, "%f\t", establish_list.est_array[ncohorts].diam);
											fprintf(ofp_init, "%f\t", diam_cohort / 10.0);		// mm to cm
											fprintf(ofp_init, "%f\t", height_cohort);		// m ?
											fprintf(ofp_init, "%f\t", dens_cohort);		// tree/ha
											// 240612
//											fprintf(ofp_init, "%f\t", agb_cohort / plot_area);						// sum of agb of trees in chohort	UNIT ??	; per ha
//											fprintf(ofp_init, "%f\t", ba_cohort / plot_area / 1000000.0);			// sum of ba of trees in chohort; mm2 to m2	; m2 per ha
											fprintf(ofp_init, "%f\t", agb_cohort);						// sum of agb of trees in chohort	tree cohort mean
											fprintf(ofp_init, "%d\t", -1);			// ba: unclear which unit, don't use for now
											fprintf(ofp_init, "%d\t", 0);			// firstyear
											fprintf(ofp_init, "%s\t", last_plotid);
//											fprintf(ofp_init, "%f\t", plot_area);	// ha
											fprintf(ofp_init, "%d\t", -1);	// ha (not known in this plot input version)
											fprintf(ofp_init, "%d\t", use_census);
											fprintf(ofp_init, "%d\t", use_census_year);
											if(use_census_closest_to_default_year)
												fprintf(ofp_init, "%d\t", default_year);
											else
												fprintf(ofp_init, "%d\t", use_census_year);
											fprintf(ofp_init, "%d\t", ntrees_cohort);
											fprintf(ofp_init, "\n");

											ncohorts++;

											if(ntrees_cohort > max_ntrees_per_plot_census)
												max_ntrees_per_plot_census = ntrees_cohort;
										}
									}
								}

								if(tree_no > max_ntrees_per_plot_total)
									max_ntrees_per_plot_total = tree_no;
								if(ncohorts > max_ncohorts)
									max_ncohorts = ncohorts;

								for(int c=0;c<census_no_max;c++) {
									if(census_nos[c])
										ncensuses2++;
								}
								if(ncensuses2 > max_ncensuses)
									max_ncensuses = ncensuses2;

								establish_list.zero();
								tree_no = 0;
								ncensuses2 = 0;
								for(int c=0;c<census_no_max;c++) {
									census_nos[c] = 0;
								}
								smallest_dif = 1000;

								for(int i=ncohorts;i<max_ncohorts_save;i++) {

									fprintf(ofp_init, "%.9f\t%.9f\t", last_lon, last_lat);
									fprintf(ofp_init, "%d\t", i);
									fprintf(ofp_init, "%s\t", "-1");
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%s\t", last_plotid);
									fprintf(ofp_init, "%d\t", -1);
									if(use_census && !i)
										fprintf(ofp_init, "%d\t", use_census);		// print selected census_no that has no living trees
									else
										fprintf(ofp_init, "%d\t", -1);
									if(use_census_year && !i)
										fprintf(ofp_init, "%d\t", use_census_year);
									else
										fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "%d\t", -1);
									fprintf(ofp_init, "\n");
								}

								use_census = 0;
								use_census_year = 0;

								if(!ncohorts) {
									printf("No trees printed for plot %s\n", last_plotid);
									fprintf(ofp_log, "No trees printed for plot %s\n", last_plotid);
								}
							}

//							if(census_year)	{	// Ignore lines with no census year !
							if(census_year && (!dead || allow_dead_census))	{	// Ignore dead trees !
								census_nos[census_no] = 1;
								if(census_no != use_census) {
									if(use_census_closest_to_default_year) {

										double dif_to_default_year = fabs(double(census_year - default_year));
										if(dif_to_default_year < smallest_dif) {
											smallest_dif = dif_to_default_year;
											use_census = census_no;
											use_census_year = census_year;
										}
									}
									else if(!use_census || (census_no > use_census && use_last_census) || (census_no < use_census && !use_last_census)) {
										use_census = census_no;
										use_census_year = census_year;
									}
								}

								if(line_no != lastline_no)
									establish_list.add_tree(tree_no, s1[pftid_col], census_no, census_year, diameter, height, agb, ba, dens);
							}
						}
					}	// if(r == 1)
						strcpy(last_plotid, s1[plotid_col]);
						last_lon = strtod(s1[longitude_col], NULL);
						last_lat = strtod(s1[latitude_col], NULL);
				}	// if(!firstline)
			
				line_no++;
			}	// if(fgets(line,sizeof(line),ifp))
		}	// while(!feof(ifp))

		if(!lap_no) {
			rewind(ifp);
			if(max_ncohorts > MAX_NESTABLISH_INDIV) {	// This should never happen anymore !
				printf("Too many cohorts: Increase MAX_NESTABLISH_INDIV to %d!\n", max_ncohorts);
				fprintf(ofp_log, "Too many cohorts: Increase MAX_NESTABLISH_INDIV to %d!\n", max_ncohorts);
				if(outFileName)
					delete[] outFileName;
				exit(1);
			}
		}
	}	// for (int r=0;r<nlaps;r++)

	printf("\n");
	printf("nplots = %d\n", nplots);
	printf("max_ntrees_per_plot_census = %d\n", max_ntrees_per_plot_census);
	printf("max_ntrees_per_plot_total = %d\n", max_ntrees_per_plot_total);
	printf("max_ncohorts = %d\n", max_ncohorts);
	printf("max_ncensuses = %d\n", max_ncensuses);

	fprintf(ofp_log, "\n");
	fprintf(ofp_log, "nplots = %d\n", nplots);
	fprintf(ofp_log, "max_ntrees_per_plot_census = %d\n", max_ntrees_per_plot_census);
	fprintf(ofp_log, "max_ntrees_per_plot_total = %d\n", max_ntrees_per_plot_total);
	fprintf(ofp_log, "max_ncohorts = %d\n", max_ncohorts);
	fprintf(ofp_log, "max_ncensuses = %d\n", max_ncensuses);

	if(outFileName)
		delete[] outFileName;
	if(outFileNameInit)
		delete[] outFileNameInit;

	if(ifp)
		fclose(ifp);
	if(ifp_sp_map)
		fclose(ifp_sp_map);
	if(ofp)
		fclose(ofp);
	if(ofp_sp)
		fclose(ofp_sp);
	if(ofp_init)
		fclose(ofp_init);

	return 0;
}