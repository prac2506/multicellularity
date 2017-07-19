// compila-se da seguinte forma c++ -O3 random-neutral.cpp -o random-neutral -lm -lgsl -lgslcblas

#include<math.h>
#include<iostream>
#include<string.h>
#include<stdio.h>
#include<stdlib.h> 
#include<iomanip> 				
#include<fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>

using namespace std;

long idum;
gsl_rng *gerador;

struct grupo{
  int *number;
  int **numb_stick;
  double *fitness;
  int kmax;
  double kmed;
};


struct celula{
  double c;
  int p;
  double sigma[2];
};



struct taxa{
  double k_plus;
  double k_minus;
  double k_death;
  double k_reprod;
  double **K_AGG0;
  double **K_AGG1;
  double **K_REP;
  double **K_DEATH;
  double **K_DISSOC;
  double K_total;
  int poptotal;
  double p_single;
};

void initialization(struct taxa *rates, struct celula *cell, struct grupo *group);
void rates_evaluation(struct taxa *rates, struct celula *cell, struct grupo *group);
void reactions(struct taxa *rates, struct celula *cell, struct grupo *group);
void fitness_evaluation(struct taxa *rates, struct celula *cell, struct grupo *group);

double R1(double C, double p);

int i1=0;

int main(int ac, char **av)
{

  FILE *ptt1, *ptt2;

  grupo group;

  taxa rates;

  celula cell;

  int conf, config, t, tmax, soma, i, j, k, sum_groups;

  long seed;  //define a seed

  char arq1[200], arq2[200];

  if (ac!=12)
    {
      cout  <<  "start the program like this:\n" << av[0]
            << " <L> <U> <s> <sd> <alpha> <Nlevel> <tmax> <Tterm> <config> <Nbins> <Nmax> <semente> \n"
            << endl;
      exit (-1);
    }
  
  /**** Leitura de dados por meio de um script ****/

  j = 0;
  cell.c = atof (av[++j]);
  cell.p = atoi (av[++j]);
  rates.k_plus = atof (av[++j]);
  rates.k_minus = atof (av[++j]);
  rates.k_reprod = atof (av[++j]);
  rates.k_death = atof (av[++j]);
  rates.p_single = atof (av[++j]);
  cell.sigma[0] = atof (av[++j]);
  cell.sigma[1] = atof (av[++j]);
  tmax = atoi (av[++j]);
  config = atoi (av[++j]);

  cout  << "#invocation: ";
  for (int i=0; i<ac; i++){
    cout << av[i] << " ";
  }
  cout << endl;

  seed = time (NULL) * getpid();  //set the seed to system time
  gerador=gsl_rng_alloc(gsl_rng_mt19937);  
  gsl_rng_set(gerador, seed); //give the seed to random generator

  group.number = new int[10000];
  group.fitness = new double[10000];

  group.numb_stick = new int*[10000];
  for( k=1; k<10000; k++ )
    group.numb_stick[k] = new int[k+1];

  rates.K_AGG0 = new double*[10000];
  for( k=1; k<10000; k++ )
    rates.K_AGG0[k] = new double[k+1];

  rates.K_AGG1 = new double*[10000];
  for( k=1; k<10000; k++ )
    rates.K_AGG1[k] = new double[k+1];
  
  rates.K_REP = new double*[10000];
  for( k=1; k<10000; k++ )
    rates.K_REP[k] = new double[k+1];

  rates.K_DEATH = new double*[10000];
  for( k=1; k<10000; k++ )
    rates.K_DEATH[k] = new double[k+1];

  rates.K_DISSOC = new double*[10000];
  for( k=1; k<10000; k++ )
    rates.K_DISSOC[k] = new double[k+1];
  
  conf = 0;
  
  while( (++conf)<=config )
    {
      initialization(&rates,&cell,&group);

      fitness_evaluation(&rates,&cell,&group);
      t = 0;
      
      while( (++t<tmax) )
	{
	  rates_evaluation(&rates,&cell,&group);

	  reactions(&rates,&cell,&group);

	  if( (t%100)==0 )
	    {
	      sprintf(arq1,"COMPLEXITY-c%g-k_plus%g-k_minus%g-k_reprod%g-k_death%g.dat",cell.c,rates.k_plus,rates.k_minus,rates.k_reprod,rates.k_death);
	      ptt1 = fopen(arq1,"a");
	      
	      group.kmed = 0;
	      sum_groups = 0;
	      for( k=1; k<group.kmax; k++ )
		for( i=0; i<=k; i++ )
		  {
		    group.kmed += k*group.numb_stick[k][i];
		    sum_groups += group.numb_stick[k][i];
		    //	    cout << k << "\t" << i << "\t" << group.numb_stick[k][i] << endl;
		  }
	      group.kmed /= sum_groups;
	      fprintf(ptt1,"%d    %d    %g \n",t,group.kmax,group.kmed);
	      fclose(ptt1);
	    }  
	}
      
    } 
  
}




void initialization(taxa *rates, celula *cell, grupo *group)
{
  int i, j, k;

  group->number[1] = 10000;
  group->numb_stick[1][0] = 5000;
  group->numb_stick[1][1] = 5000;

  group->kmax = 1;

  for( j=2; j<10000; j++ )
    {
      group->number[j] = 0;
      for( k=0; k<=j; k++ )
	group->numb_stick[j][k] = 0;
    } 

}


void rates_evaluation(taxa *rates, celula *cell, grupo *group)
{
  int i, j, k, kmax;

  /***** Aggregation *****/

  kmax = group->kmax+1;
  rates->K_total = 0;
  rates->poptotal = 0;
  for( k=1; k<=group->kmax; k++ )
    for( i=0; i<=k; i++ )
      rates->poptotal += k*group->numb_stick[k][i];
  
  for( k=2; k<=kmax; k++ )  /***** this loop goes until kmax because at this point the group size can become equal to group->kmax+1 *****/
    {
      for( i=0; i<=k; i++ )
	{
	  rates->K_AGG0[k][i] = 0;
	  rates->K_AGG1[k][i] = 0;
	}

      if( k!=2 )
	{
	  for( i=0; i<=(k-1); i++ )
	    {
	      rates->K_AGG1[k][i+1] += rates->k_plus*cell->sigma[0]*group->numb_stick[k-1][i]*group->numb_stick[1][1]*pow((k-1.),(2./3))*(i*cell->sigma[0]+(k-i-1)*cell->sigma[1])/(k-1.);
	      rates->K_AGG0[k][i] += rates->k_plus*cell->sigma[1]*group->numb_stick[k-1][i]*group->numb_stick[1][0]*pow((k-1.),(2./3))*(i*cell->sigma[0]+(k-i-1)*cell->sigma[1])/(k-1.);
	    }
	}

      if( k==2 )
	{
	  rates->K_AGG0[2][0] = rates->k_plus*cell->sigma[1]*group->numb_stick[1][0]*(group->numb_stick[1][0]-1)*cell->sigma[1];
	  rates->K_AGG1[2][2] = rates->k_plus*cell->sigma[0]*group->numb_stick[1][1]*(group->numb_stick[1][1]-1)*cell->sigma[0];
	  rates->K_AGG1[2][1] = rates->k_plus*cell->sigma[0]*group->numb_stick[1][1]*group->numb_stick[1][0]*cell->sigma[1];
	}

      for( i=0; i<=k; i++ )
	rates->K_total += rates->K_AGG0[k][i] + rates->K_AGG1[k][i];

    }
  /***** Reproduction *****/

  for( k=1; k<kmax; k++ )
    for( i=0; i<=k; i++ )
      {
	rates->K_REP[k][i] = rates->k_reprod*group->fitness[k]*group->numb_stick[k][i]*(1-(i*cell->sigma[0]+(k-i)*cell->sigma[1])/(k*10));
	rates->K_total += rates->K_REP[k][i];       
      }

  /***** Death *****/

  for( k=1; k<kmax; k++ )
    for( i=0; i<=k; i++ )
      {
	rates->K_DEATH[k][i] = rates->k_death*rates->poptotal*group->numb_stick[k][i]*k;
	rates->K_total += rates->K_DEATH[k][i];
      }

  /***** Dissociation *****/

  for( k=2; k<kmax; k++ )
    for( i=0; i<=k; i++ )
      {
	rates->K_DISSOC[k][i] = rates->k_minus*group->numb_stick[k][i]*pow(k,(2./3));
	rates->K_total +=  rates->K_DISSOC[k][i];
      }
  
}

void reactions(taxa *rates, celula *cell, grupo *group)
{
  int i, j, k, kmax, sum_check;
  double sum, r, r1, r2;
   
  kmax = group->kmax + 1;

  r = gsl_ran_flat(gerador, 0., 1.);
  sum = 0;
  for( k=2; k<=kmax; k++ )  /***** Aggregation *****/
    for( i=0; i<=k; i++ )
      {
	sum += rates->K_AGG0[k][i]/rates->K_total; /***** K_AGG0 - o agregado tem sigma[1]  *****/
	if( sum>r )
	  {
	    group->numb_stick[k][i]++;
	    group->numb_stick[k-1][i]--;
	    group->numb_stick[1][0]--;
	    k = kmax + 1;
	    
	    break;
	  }

	sum += rates->K_AGG1[k][i]/rates->K_total;
	if( sum>r )
	  {
	    group->numb_stick[k][i]++;
	    group->numb_stick[k-1][i-1]--;
	    group->numb_stick[1][1]--;
	    k = kmax + 1;

	    break;
	  }

      }

  if( sum<r )  /***** Reproduction *****/
    {
      for( k=1; k<kmax; k++ )
	for( i=0; i<=k; i++ )
	  {
	    sum += rates->K_REP[k][i]/rates->K_total;
	    if( sum>r )
	      {
		r2 = gsl_ran_flat(gerador, 0., 1.);
		if( r2<rates->p_single )
		  {
		    r1 = gsl_ran_flat(gerador, 0., 1.);
		    if( r1<((double)i/k) )
		      group->numb_stick[1][1]++;
		    else
		      group->numb_stick[1][0]++;
		    k = kmax + 1;
		    break;
		  }
		else
		  {
		    group->numb_stick[k][i]--;
		    r1 = gsl_ran_flat(gerador, 0., 1.);
		    if( r1<((double)i/k) )
		      group->numb_stick[k+1][i+1]++;
		    else
		      group->numb_stick[k+1][i]++;
		    k = kmax + 1;
		    break;
		  }
	      }
	  }
    }

  if( sum<r )   /***** Death *****/
    {
      for( k=1; k<kmax; k++ )
	for( i=0; i<=k; i++ )
	  {
	    sum += rates->K_DEATH[k][i]/rates->K_total;
	    if( sum>r )
	      {
		r1 = gsl_ran_flat(gerador, 0., 1.);
		if( r1<((double)i/k) )
		  {
		    group->numb_stick[k][i]--;
		    if( k>1 )
		      group->numb_stick[k-1][i-1]++;
		  }
		else
		  {
		    group->numb_stick[k][i]--;
		    if( k>1 )
		      group->numb_stick[k-1][i]++;
		  }
		k = kmax + 1;
		break;
	      }
	  }
    }

  if( sum<r )   /***** Dissociation *****/
    {
      for( k=2; k<kmax; k++ )
	for( i=0; i<=k; i++ )
	  {
	    sum += rates->K_DISSOC[k][i]/rates->K_total;
	    if( sum>r )
	      {
		r1 = gsl_ran_flat(gerador, 0., 1.);
		if( r1<((double)i/k) )
		  {
		    group->numb_stick[k][i]--;
		    group->numb_stick[k-1][i-1]++;
		    group->numb_stick[1][1]++;
		  }
		else
		  {
		    group->numb_stick[k][i]--;
		    group->numb_stick[k-1][i]++;
		    group->numb_stick[1][0]++;
		  }

		k = kmax+1;
		break;
		
	      }
	    
	  }
    }

  group->kmax = 0;
  for( k=1; k<=(2*kmax); k++ )
    for( i=0; i<=k; i++ )
      if( (group->numb_stick[k][i]>0) && (k>group->kmax) )
	group->kmax = k;
}


void fitness_evaluation(taxa *rates, celula *cell, grupo *group)
{
  int i;
  int p = cell->p;
  double undifferentiated_fitness = R1(cell->c, double(cell->p));

  for( i=1; i<p; i++ )
  {
    group->fitness[i] = i*undifferentiated_fitness;
  }

  
  for( i=p; i<=10000; i++ )
  {
    double k = i%p;
    double m = i/p;
    group->fitness[i] = i*pow(double(p),4.)*pow(m+1.,4.*k/p)*pow(m,4.*(p-k)/p)/(pow(double(i),4.)*432.*cell->c*cell->c);
  }
}

//Newton's method
double R1(double C, double p) {
  double x = 1.;
  for (int i = 0; i < 10; i++)
    x = x - (exp((1.-p)*x*x)*p*p + 4.*C*x*(-3. + 2.*(1.-p)*x*x))/(-2.*((p-1.)*x*exp((1.-p)*x*x)*p*p + 6.*C*(1. + 2.*(p-1.)*x*x)));
  return x*x*(1. - 8.*C/(p*p)*x*exp((p-1.)*x*x));
}


