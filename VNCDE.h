/* VNCDE */

#ifndef VNCDE_H
#define VNCDE_H

CEC2013 *pFunc;		

#define MAXDIM  20
#define MAXFES  400000
#define POPSIZE 500

/*parameters of the algorithm*/
int popsize = 100;  
double F = 0.5;
double CR = 0.9;
int nn_list[POPSIZE];           /* Voronoi neighbor list */
int nn_mat[POPSIZE][POPSIZE];   /* Voronoi neighbor matrix */
int nn_size[POPSIZE];
int nn_better[POPSIZE][POPSIZE];
int nn_bsize[POPSIZE];
int type[POPSIZE];   
double pi = 3.1415926535897932;

struct Individual
{
	double x[MAXDIM];
	double u[MAXDIM];
	double fit;
	double ufit;
	double theta;
	int    limit;
};
Individual pop[POPSIZE];
Individual archive[POPSIZE];
bool replace[POPSIZE];
int archive_size;

double gbestval;            //use the found global best fitness to estimate the global optima
int global_update;

#define betterthan >
#define worsethan  <
const double accuracy_level[] = {0.1, 0.01, 0.001, 0.0001, 0.00001};//five levels of accuracy 
int dim;
int maxfes,fes;
double xmax[MAXDIM];
double xmin[MAXDIM];

vector<vector<double>> solu; 
vector<vector<double>> seed; //solu & seed are used in the final output of the algorithm

/*performance measure*/
double peak_ratio[5];
double succ_rate[5];
double run_time;

inline double rand_uni(double low,double high);
inline double rand_gau(double mu, double thegma);
double rand_cau(double alpha, double beta);
inline void boundsctl(double & x_id, int d);
inline int best_in_niche(int i);
inline int worst_in_niche(int i);
inline double eu_dis(double p1[], double p2[], int dd);
void cal_niche_dis();
int find_arc_worst();
inline void reinitialize_ind(int i);
void reinitialize_niche(int i) ;
void outputsolu();
void add_to_archive(int i);
void add_to_archive(Individual& p);
double caleval(double pos[]);
int cmp ( const void *a , const void *b );
double cos_theta(Individual& p1, Individual& p2);

#endif
