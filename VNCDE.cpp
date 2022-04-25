/*	VNCDE */

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<ctime>
#include<math.h>
using std::vector;

#include"cec2013.h"
#include"VNCDE.h"

void initialization()
{
	int i, j;
	fes = 0;
	gbestval = -1e300;

	for(i=0; i<popsize; i++){
		for(j=0; j<dim; j++)
			pop[i].x[j] = rand_uni(xmin[j], xmax[j]);
		pop[i].fit = caleval(pop[i].x);
		if(pop[i].fit > gbestval)
			gbestval = pop[i].fit;
		pop[i].limit = 0;
	}
	global_update = 0;

	archive_size = 100;
	for(i=0; i<archive_size; i++){
		for(j=0; j<dim; j++)
			archive[i].x[j] =  rand_uni(xmin[j], xmax[j]);
		archive[i].fit = caleval(archive[i].x);
	}
	memset(replace, 0, sizeof(replace));
}

int find_nearest_v(Individual& ind, bool used[])
{
	int nei;
	double mindis = 1e300, dis;
	for(int j=0; j<popsize; j++){
		if(used[j]==1)continue;
		dis = eu_dis(ind.x, pop[j].x, dim);
		if(dis < mindis){
			mindis = dis;
			nei = j;
		}
	}
	return nei;
}

int find_nearest_v_ij(Individual& ind, bool used[], int p1, int p2)
{
	int nei;
	double mindis = 1e300, dis;
	for(int j=0; j<popsize; j++){
		if(used[j]==1)continue;
		dis = eu_dis(ind.x, pop[j].x, dim);
		if(dis < mindis){
			mindis = dis;
			nei = j;
		}
		else if(dis-mindis==0 && (j==p1 || j==p2)){
			nei = j;
		}
	}
	return nei;
}

int find_neighbors_v(int i)
{
	int nsize = 0;
	int nei;
	bool used[POPSIZE] = {0};
	Individual mid;
	for(int j=0; j<popsize; j++){
		if(j==i)continue;
		for(int k=0; k<dim; k++){
			mid.x[k] = (pop[i].x[k] + pop[j].x[k])/2.0;
		}
		//nei = find_nearest_v(mid, used);
		nei = find_nearest_v_ij(mid, used, i, j);
		if(nei==i || nei==j){
			nn_list[nsize] = j;
			nsize++;
		}
	}
	return nsize;
}

int compare_fit(const int& a, const int& b)
{
	return pop[a].fit > pop[b].fit;
}

int compare_type(const int& a, const int& b)
{
	if(type[a]==0 && type[b]!=0)
		return 1;
	else if(type[b]==0 && type[a]!=0)
		return -1;
	else 
		return 0;
}

double cos_theta(Individual& p1, Individual& p2)
{
	double pro = 0;
	double norm1 = 0;
	double norm2 = 0;
	for(int j=0; j<dim; j++){
		pro = p1.x[j]*p2.x[j];
		norm1 += p1.x[j]*p1.x[j];
		norm2 += p2.x[j]*p2.x[j];
	}
	norm1 = sqrt(norm1);
	norm2 = sqrt(norm2);
	return pro/(norm1*norm2);
}

int define_type(int i, int neighbors[], int nsize)
{
	Individual vec[300];
	for(int j=0; j<nsize; j++){
		for(int k=0; k<dim; k++){
			vec[j].x[k] = pop[neighbors[j]].x[k]-pop[i].x[k];
		}
	}
	double avg = 0;
	int cnt = 0;
	double min = pi;
	double max = 0;
	for(int j=0; j<nsize; j++){
		for(int k=j+1; k<nsize; k++){
			double cs = cos_theta(vec[j], vec[k]);
			cs = acos(cs);
			if(min>cs) min = cs;
			if(cs>max) max = cs;
			avg += cs;
			cnt++;
		}
	}
	if(cnt==0) avg = 0; else avg /= cnt;
	pop[i].theta = avg;

	if(pop[i].theta<pi/2.0)
		return 1;
	else
		return 2;
}

void find_neighbors_all()
{
	for(int i=0; i<popsize; i++){
		nn_size[i] = find_neighbors_v(i);
		std::sort(nn_list, nn_list+nn_size[i], compare_fit);
		for(int j=0; j<nn_size[i]; j++){
			nn_mat[i][j] = nn_list[j];
		}
	}

	// sort the indivdiuals according to their fitness values
	int ind[POPSIZE];
	for(int i=0; i<popsize; i++)
		ind[i] = i;
	std::sort(ind, ind+popsize, compare_fit);
	
	// define type
	for(int i=0; i<popsize; i++){
		int pt = ind[i];
		int nbest = nn_mat[pt][0];
		int nbetter[300];
		int nsize = 0;
		
		for(int k=0; k<nn_size[pt]; k++){
			if(pop[nn_mat[pt][k]].fit > pop[pt].fit){
				nbetter[nsize] = nn_mat[pt][k];
				nsize++;
			}
		}
		if(nsize==0){
			type[pt] = 0;
			pop[pt].theta = rand_uni(0, pi);

			nn_bsize[pt] = 0;
		}
		else
		{
			type[pt] = define_type(pt, nbetter, nsize);
			nn_bsize[pt] = nsize;
			for(int k=0; k<nsize; k++)
				nn_better[pt][k] = nbetter[k];
		}

	}
}

double FindBestfit()
{
	double bestfit = pop[0].fit;
	for(int i=1; i<popsize; i++){
		if(pop[i].fit > bestfit)
			bestfit = pop[i].fit;
	}
	return bestfit;
}

double FindWorstfit()
{
	double worstfit = pop[0].fit;
	for(int i=1; i<popsize; i++){
		if(pop[i].fit < worstfit)
			worstfit = pop[i].fit;
	}
	return worstfit;
}

void crossover_mutation(int i)
{
	int r1, r2, r3, jrand;
	int nsize = nn_size[i];
	jrand=rand()%dim;

	if(type[i]==0){ // dominators  
		double sigma;
		sigma = pow(0.1, 5.0+rand()%5-3);
		double rg;
		for(int j=0; j<dim; j++){
			if(rand_uni(0,1)<=CR || j==jrand)
			{
				do{rg=rand_gau(0, sigma);}while(rg>1||rg<-1);
				pop[i].u[j]= pop[i].x[j] + rg; // local search
				boundsctl(pop[i].u[j],j);
			}
			else
				pop[i].u[j]=pop[i].x[j];
		}
	}
	else if(type[i]==1 || nsize<3){ // challengers
		if(nsize>=3)
		{
			int ni = 0;
			double sum = 0;
			double prob[POPSIZE];
			double bestfit = FindBestfit();
			double worstfit = FindWorstfit();
			for(int k=0; k<nn_bsize[i]; k++){
				sum += (pop[nn_better[i][k]].fit-worstfit)/(bestfit-worstfit+1E-10);
			}
			for(int k=0; k<nn_bsize[i]; k++){
				prob[k] = (pop[nn_better[i][k]].fit-worstfit)/(bestfit-worstfit+1E-10)/sum;
			}
			double partial_sum = prob[0];
			double r = rand_uni(0, 1.0);
			while(partial_sum<r){
				if(ni==nn_bsize[i]-1)break;
				ni++;
				partial_sum += prob[ni];
			}
			ni = nn_better[i][ni];
		
			r1=rand()%nsize;
			r2=rand()%nsize;
			r1=nn_mat[i][r1];
			r2=nn_mat[i][r2];

			double rg, rg2;
			for(int j=0; j<dim; j++){
				if(rand_uni(0,1)<=CR || j==jrand)
				{
					do{rg=rand_gau(0, 0.5);}while(rg>1||rg<-1);
					do{rg2=rand_gau(0, 1E-05);}while(rg2>1||rg2<-1);
					pop[i].u[j]= pop[i].x[j] + (pop[ni].x[j]-pop[i].x[j])*rg + (pop[r1].x[j]-pop[r2].x[j])*rg2; 
					boundsctl(pop[i].u[j],j);
				}
				else 
					pop[i].u[j]=pop[i].x[j];
			}	
		}
		else
		{
			int ni = 0;
			double sum = 0;
			double prob[POPSIZE];
			double bestfit = FindBestfit();
			double worstfit = FindWorstfit();
			for(int k=0; k<nsize; k++){
				sum += (pop[nn_mat[i][k]].fit-worstfit)/(bestfit-worstfit+1E-10);
			}
			for(int k=0; k<nsize; k++){
				prob[k] = (pop[nn_mat[i][k]].fit-worstfit)/(bestfit-worstfit+1E-10)/sum;
			}
			double partial_sum = prob[0];
			double r = rand_uni(0, 1.0);
			while(partial_sum<r){
				if(ni==nsize-1)break;
				ni++;
				partial_sum += prob[ni];
			}
			ni = nn_mat[i][ni];

			double rg;
			for(int j=0; j<dim; j++){
				if(rand_uni(0,1)<=CR || j==jrand)
				{
					do{rg=rand_gau(0, 0.5);}while(rg>1||rg<-1);
					pop[i].u[j]= pop[i].x[j] + (pop[ni].x[j]-pop[i].x[j])*rg; 
					boundsctl(pop[i].u[j],j);
				}
				else 
					pop[i].u[j]=pop[i].x[j];
			}	
		}
	}
	else if(type[i]==2){  // explorers
		r1=rand()%nsize;
		do{r2=rand()%nsize;}while(r2==r1);
		do{r3=rand()%nsize;}while(r3==r2 || r3==r1 );
		r1 = nn_mat[i][r1];
		r2 = nn_mat[i][r2];
		r3 = nn_mat[i][r3];

		for(int j=0; j<dim; j++){
			if(rand_uni(0,1)<=CR || j==jrand){
				pop[i].u[j]= pop[r1].x[j]+F*(pop[r2].x[j]-pop[r3].x[j]); 
				boundsctl(pop[i].u[j],j);
			}
			else 
				pop[i].u[j]=pop[i].x[j];
		}	
	}
}

int find_nearest(int i)
{
	int nei = 0;
	double mindis = 1e300, dis;
	for(int j=0; j<popsize; j++){
		dis = eu_dis(pop[i].u, pop[j].x, dim);
		if(dis < mindis){
			mindis = dis;
			nei = j;
		}
	}
	return nei;
}

int find_nearest_in_archive(Individual& p)
{
	int nei;
	double mindis = 1e300, dis;
	for(int j=0; j<archive_size; j++){
		dis = eu_dis(p.x, archive[j].x, dim);
		if(dis < mindis){
			mindis = dis;
			nei = j;
		}
	}
	return nei;
}

void add_to_archive(Individual& p)
{
	int nn = find_nearest_in_archive(p);
	if(archive[nn].fit < p.fit){
		archive[nn] = p;
		replace[nn] = 1;
	}
}

void selection()
{
	int i, j;
	find_neighbors_all();

	// sort the indivdiuals according to their types
	int ind[POPSIZE];
	for(i=0; i<popsize; i++)
		ind[i] = i;
	std::sort(ind, ind+popsize, compare_type);

	for(int rank=0; rank<popsize; rank++)
	{
		i = ind[rank];
		crossover_mutation(i);
		for(j=0; j<dim; j++)
			boundsctl(pop[i].u[j], j);

		pop[i].ufit = caleval(pop[i].u);
		int nn = find_nearest(i);
		if(pop[i].ufit > pop[nn].fit){
			add_to_archive(pop[nn]);
			for(j=0; j<dim; j++)
				pop[nn].x[j] = pop[i].u[j];
			pop[nn].fit = pop[i].ufit;
			if(type[i]==0)
				pop[i].limit=0;
		}
		else{
			if(type[i]==0)
				pop[i].limit++;
		}
	}

	// re-initialization, a method to avoid stagnation
	int sid = 0;
	for(i=0; i<popsize; i++){
		if(type[i]==0 && ((pop[i].limit > 0.1*popsize && fabs(pop[i].fit-gbestval)<1e-07)||pop[i].limit>2*popsize )){
			pop[i].limit = 0;
			add_to_archive(pop[i]);
			while(sid<archive_size && replace[sid]==1)sid++;
			if(sid>=archive_size){
				for(j=0; j<dim; j++)
					pop[i].x[j] = rand_uni(xmin[j], xmax[j]);
				pop[i].fit = caleval(pop[i].x);
			}
			else{
				for(j=0; j<dim; j++)
					pop[i].x[j] = archive[sid].x[j];
				pop[i].fit = archive[sid].fit;
				replace[sid] = 1;
			}
			int nsize = nn_size[i];
			for(j=0; j<nsize; j++){
				int id = nn_mat[i][j];
				add_to_archive(pop[id]);
				while(sid<archive_size && replace[sid]==1)sid++;
				if(sid>=archive_size){
					for(int k=0; k<dim; k++)
						pop[id].x[k] = rand_uni(xmin[k], xmax[k]);
					pop[id].fit = caleval(pop[id].x);
				}
				else{
					for(int k=0; k<dim; k++)
						pop[id].x[k] = archive[sid].x[j];
					pop[id].fit = archive[sid].fit;
					replace[sid] = 1;
				}
			}
		}
	}
}

void process()
{
	initialization();
	while(fes < maxfes)
	{	
		selection();	
	}
	printf("num global update: %d\n", global_update);
}

void main()
{	
	int start, end, runtime;
	FILE* input_param = fopen("parameters.txt", "r");
	if(input_param!=NULL)
	{
		fscanf(input_param, "%d", &start);
		fscanf(input_param, "%d", &end);
		fscanf(input_param, "%d", &runtime);
		fclose(input_param);
	}
	else{
		start = 1;
		end = 20;
		runtime = 50;
	}

	int func_id, d;
	int nkp;    //the number of known global optima
	int nfp[5]; //the number of found global optima (with different accuracy levels)
	double start_time, end_time;

	srand((unsigned)time(NULL));
	FILE* total_pr = fopen("total_pr.txt", "a");
	FILE* total_sr = fopen("total_sr.txt", "a");
    fprintf(total_pr,"func\taccruacy level\t1.0E-01\t1.0E-02\t1.0E-03\t1.0E-04\t1.0E-05\n");
	fclose(total_pr);
	fprintf(total_sr,"func\taccruacy level\t1.0E-01\t1.0E-02\t1.0E-03\t1.0E-04\t1.0E-05\n");
	fclose(total_sr);
	
	for(func_id = start; func_id <= end; func_id ++){
		/*initialize benchmark instance*/
		pFunc = new CEC2013(func_id);
		dim = pFunc->get_dimension();	
		maxfes = pFunc->get_maxfes();
		for(d=0;d<dim;++d) {
			xmax[d]=pFunc->get_ubound(d);
			xmin[d]=pFunc->get_lbound(d);
		}
		printf("F%d is running...\n",func_id);
		total_pr = fopen("total_pr.txt", "a");
		fprintf(total_pr, "%d\t", func_id);
		total_sr = fopen("total_sr.txt", "a");
		fprintf(total_sr, "%d\t", func_id);

		FILE* stat_pr = fopen("stat_pr.txt", "a");
		FILE* stat_np = fopen("stat_np.txt", "a");

		/*run algorithm & output result*/
		nkp = pFunc->get_no_goptima();
		printf("%d\n", nkp); 
		memset(peak_ratio, 0, sizeof(peak_ratio));
		memset(succ_rate, 0, sizeof(succ_rate));
		run_time = 0;
		for(int i=0; i<runtime; i++){
			start_time = clock();
			process();
			end_time = clock();
			run_time = (end_time-start_time)/double(CLOCKS_PER_SEC);
			printf("run time: %lf\n", run_time);

			outputsolu();
			for(int j=0; j<5; j++){ //five levels of accuracy
				nfp[j] = how_many_goptima(solu,seed,pFunc,accuracy_level[j],pFunc->get_rho());
				printf("%d\t", nfp[j]);
				fprintf(stat_np, "%d\t", nfp[j]);
				fprintf(stat_pr, "%lf\t", (double)nfp[j]/nkp);
				peak_ratio[j] += nfp[j];
				if(nfp[j] >= nkp)succ_rate[j]++;
			}
			printf("\n"); 
			fprintf(stat_np, "\n");
			fprintf(stat_pr, "\n");
		}
		fprintf(stat_np, "\n");
		fprintf(stat_pr, "\n");

		fprintf(total_pr, "peak_ratio\t");
		for(int j=0; j<5; j++){
			peak_ratio[j] /= (runtime*nkp + 0.0);
			succ_rate[j] /= (runtime + 0.0);
			printf("level %d : %f\t%f\n", j+1, peak_ratio[j],succ_rate[j]);
			fprintf(total_pr, "%f\t", peak_ratio[j]);
		}
		fprintf(total_pr, "\n");
		fprintf(total_sr, "\tsuccess_rate\t");
		for(int j=0; j<5; j++)fprintf(total_sr, "%f\t", succ_rate[j]);
		fprintf(total_sr, "\n");
		
		fclose(total_pr);
		fclose(total_sr);
		fclose(stat_pr);
		fclose(stat_np);

		delete pFunc;
	}
	system("Pause");
}

/*           auxiliary functions          */

inline double rand_uni(double low,double high)//generate uniformly random numbers
{
	return (double(rand())/RAND_MAX)*(high-low)+low;
}

static int phase = 0;
inline double rand_gau(double mu, double thegma)//generate Gaussian distributed random numbers
{
	static double V1, V2, S; 
	double X;    
	if ( phase == 0 ) {
		do {
			double U1 = (double)rand() /(double)RAND_MAX;
			double U2 = (double)rand() /(double)RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	return mu+X*thegma;
}

/* alpha: position parameter, beta: scale parameter */
double rand_cau(double alpha, double beta)
{
	#define PI 3.1415926535897932384626433832795029
	return alpha+beta*tan(PI*(double(rand())/RAND_MAX-0.5));
}

inline void boundsctl(double & x_id, int d)//bound control, to restrict the variable in the range
{
	if(rand()%20<1){
		if(x_id < xmin[d])
			x_id = xmin[d];
		else if(x_id > xmax[d])
			x_id = xmax[d];
	}
	else{
		if(x_id<xmin[d]){
			x_id = 2*xmin[d] - x_id;
			if(x_id>xmax[d])x_id = xmax[d];
		}
		else if(x_id>xmax[d]){
			x_id=2*xmax[d] - x_id;
			if(x_id<xmin[d]) x_id=xmin[d];
		}
	}
}

inline double eu_dis(double p1[], double p2[], int dd)
{
	double dis = 0;
	for(int d=0; d<dd; d++)
		dis += (p1[d]-p2[d])*(p1[d]-p2[d]);
	return sqrt(dis);
}

void outputsolu()
{
	vector<double> s;
	solu.clear();
	for(int i=0; i<popsize; i++){
		s.clear();
		for(int j=0; j< dim; j++) s.push_back(pop[i].x[j]);
		solu.push_back(s);
	}
	for(int i=0; i<archive_size; i++){
		s.clear();
		for(int j=0; j< dim; j++) s.push_back(archive[i].x[j]);
		solu.push_back(s);
	}
}

double caleval(double pos[])
{
	fes++;
	double fit = pFunc->evaluate(pos);
	if(fit > gbestval) global_update++;
	if(fit > gbestval) gbestval = fit;
	return fit;
}

int cmp ( const void *a , const void *b ) 
{ 
	struct Individual * p = (struct Individual *) a;
	struct Individual * q = (struct Individual *) b;
	return (p->fit - q->fit < 0) ? 1 : -1; ; 
} 