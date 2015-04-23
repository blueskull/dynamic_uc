#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>

//#define DEBUG

typedef struct
{
	float pmin;
	float pmax;
	float a, b, c;
}GEN;

typedef struct
{
	float pmin;
	float pmax;
	int config;
}COM;

float *load;
GEN *gen;
COM *com;
int *state;
int loads;
int gens;
int coms;

//Read in generator configuration data
int readgen(char *fn)
{
	if(!fn) return -1;
	FILE *fp;
	fp=fopen(fn, "r");
	if(!fp) return -1;
	size_t l=0;
	int i=0;
	char *s;
	int gencnt;
	while(getline(&s, &l, fp)!=-1)
	{
		if(s[0]=='#') continue;
		if(i==0)
		{
			if(sscanf(s, "%d", &gencnt)!=1)
			{
				fclose(fp);
				return -1;
			}
			gen=(GEN*)malloc(gencnt*sizeof(GEN));
			if(!gen || gencnt==0)
			{
				fclose(fp);
				return -1;
			}
#ifdef DEBUG
			printf("Totally %d generators:\n", gencnt);
#endif
			i++;
			continue;
		}
		if(i-1>=gencnt)
		{
			fclose(fp);
			return -1;
		}
		if(sscanf(s, "%f,%f,%f,%f,%f", &gen[i-1].pmin, &gen[i-1].pmax, &gen[i-1].a, &gen[i-1].b, &gen[i-1].c)!=5)
		{
			fclose(fp);
			return -1;
		}
#ifdef DEBUG
		printf("Generator %2d:\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\n", i, gen[i-1].pmin, gen[i-1].pmax, gen[i-1].a, gen[i-1].b, gen[i-1].c);
#endif
		i++;
	}
	fclose(fp);
	if(i-1!=gencnt) return -1;
	gens=gencnt;
	return 0;
}

//Read in load and renewable data and calculate net load
int readload(char *fn)
{
	if(!fn) return -1;
	FILE *fp;
	fp=fopen(fn, "r");
	if(!fp) return -1;
	size_t l=0;
	int i=0, j;
	char *s;
	int loadcnt;
	int renewcnt;
	float tmp[4];
	while(getline(&s, &l, fp)!=-1)
	{
		if(s[0]=='#') continue;
		if(i==0)
		{
			if(sscanf(s, "%d,%d", &loadcnt, &renewcnt)!=2)
			{
				fclose(fp);
				return -1;
			}
			load=(float*)malloc(loadcnt*sizeof(float));
			state=(int*)malloc(loadcnt*sizeof(int));
			if(!load || !state || loadcnt==0)
			{
				fclose(fp);
				return -1;
			}
#ifdef DEBUG
			printf("Totally %d load points and %d renewable sources:\n", loadcnt, renewcnt);
#endif
			i++;
			continue;
		}
		if(i-1>loadcnt)
		{
			fclose(fp);
			return -1;
		}
		if(sscanf(s, "%f,%f,%f,%f,%f", &load[i-1], &tmp[0], &tmp[1], &tmp[2], &tmp[3])!=renewcnt+1)
		{
			fclose(fp);
			return -1;
		}
		for(j=0;j<renewcnt;j++)
			load[i-1]-=tmp[j];
		if(load[i-1]<0) load[i-1]=0;
#ifdef DEBUG
		printf("Time %7d:\t%.2f\n", i, load[i-1]);
#endif
		i++;
	}
	fclose(fp);
	if(i-1!=loadcnt) return -1;
	loads=loadcnt;
	return 0;
}

//Enumerate all generator combinations
int enumcom()
{
	int i, j;
	coms=1;
	for(i=0;i<gens;i++)
		coms*=2;
#ifdef DEBUG
	printf("Totally %d states are possible:\n", coms);
#endif
	com=(COM*)malloc(coms*sizeof(COM));
	for(i=0;i<coms;i++)
	{
		com[i].config=i;
		com[i].pmin=0;
		com[i].pmax=0;
		for(j=0;j<gens;j++)
			if(com[i].config&(1<<j))
			{
				com[i].pmin+=gen[j].pmin;
				com[i].pmax+=gen[j].pmax;
			}
#ifdef DEBUG
		printf("State %6d:\t%.2f\t%.2f\n", i+1, com[i].pmin, com[i].pmax);
#endif
	}
	if(!com) return -1;
	return 0;
}

//Do ecomonic dispatch, return minimum cost to supply the demand
float ecodispatch(int combo, float demand, float *pwr)
{
	if(demand>com[combo].pmax || demand<com[combo].pmin) return FLT_MAX;
	float lambda;
	int i;
	float *pset=(float*)malloc(gens*sizeof(float));
	int *flag=(int*)malloc(gens*sizeof(int));
	if(!pset || !flag) return -1.0;
	int done=0, valid;
	float sum_b, sum_a, sum_p;
	float cost=0;
	for(i=0;i<gens;i++)
		flag[i]=0;
	while(1)
	{
		valid=0;
		sum_b=0;
		sum_a=0;
		sum_p=0;
		for(i=0;i<gens;i++)
		{
			if(flag[i]==1 || (com[combo].config&(1<<i))==0)
				continue;
			sum_b+=gen[i].b/2.0/gen[i].a;
			sum_a+=1.0/2.0/gen[i].a;
		}
		lambda=(demand+sum_b)/sum_a;
		for(i=0;i<gens;i++)
		{
			if(flag[i]==1 || (com[combo].config&(1<<i))==0)
				continue;
			valid=1;
			pset[i]=(lambda-gen[i].b)/2.0/gen[i].a;
			if(pset[i]<gen[i].pmin)
			{
				pset[i]=gen[i].pmin;
				flag[i]=1;
				demand-=pset[i];
				break;
			}
			else if(pset[i]>gen[i].pmax)
			{
				pset[i]=gen[i].pmax;
				flag[i]=1;
				demand-=pset[i];
				break;
			}
			sum_p+=pset[i];
		}
		if(fabs(sum_p-demand)<0.1)
		{
			done=1;
			break;
		}
		if(valid==0) break;
	}
	for(i=0;i<gens;i++)
	{
		if((com[combo].config&(1<<i))==0) continue;
		cost+=gen[i].c+gen[i].b*pset[i]+gen[i].a*pset[i]*pset[i];
	}
	if(done && pwr!=NULL)
		memcpy(pwr, pset, gens*sizeof(float));
	free((void*)pset);
	if(done) return cost;
	return FLT_MAX;
}

//Unit commitment recursive kernel
int uckernel(int time)
{
	if(time==loads) return 1;
	float *cost=(float*)malloc(coms*sizeof(float));
	int i;
	for(i=0;i<coms;i++)
		cost[i]=ecodispatch(i, load[time], NULL);
	float cmin;
	while(1)
	{
		cmin=FLT_MAX;
		for(i=0;i<coms;i++)
			if(cost[i]<cmin)
			{
				cmin=cost[i];
				state[time]=i;
			}
		if(cmin==FLT_MAX)
		{
			free((void*)cost);
			return 0;
		}
		if(uckernel(time+1)==1)
		{
			free((void*)cost);
			return 1;
		}
		cost[state[time]]=FLT_MAX;
	}
}

//Unit commitment
int unitcmt()
{
	return uckernel(0);
}

void clean(int retval)
{
	if(gen) free((void*)gen);
	if(com) free((void*)com);
	if(load) free((void*)load);
	if(state) free((void*)state);
	exit(retval);
}

int main(int argc, char **argv)
{
	if(argc!=5)
	{
		printf("Invalid syntax.\nUsage:\t%s <config_in> <load_in> <uc_out> <ed_out>\n", argv[0]);
		return -1;
	}
	if(readgen(argv[1]))
	{
		printf("Can not read in generator config data.\n");
		clean(-1);
	}
	if(readload(argv[2]))
	{
		printf("Can not read in load and renewable data.\n");
		clean(-1);
	}
	if(enumcom())
	{
		printf("Insufficient memory.\n");
		clean(-1);
	}
	if(unitcmt()!=1)
	{
		printf("Can not find a solution.\n");
		clean(-1);
	}
#ifdef DEBUG
	printf("Unit commitment results:\n");
#endif
	int i, j;
	float sum_c=0;
	for(i=0;i<loads;i++)
	{
#ifdef DEBUG
		printf("Time %7d:\t", i+1);
		for(j=0;j<gens;j++)
		{
			if(com[state[i]].config&(1<<j))
				putchar('x');
			else
				putchar('o');
		}
		putchar('\n');
#endif
		sum_c+=ecodispatch(state[i], load[i], NULL);
	}
#ifdef DEBUG
	printf("Total cost: %.2f\n", sum_c);
#endif
	FILE *fp=fopen(argv[3], "w");
	if(!fp)
	{
		printf("Can not open output file.\n");
		clean(-1);
	}
	fprintf(fp, "%.2f", sum_c);
	for(i=0;i<loads;i++)
	{
		fputc('\n', fp);
		for(j=0;j<gens;j++)
		{
			if(j!=0) fputc(',', fp);
			if(com[state[i]].config&(1<<j))
				fputc('1', fp);
			else
				fputc('0', fp);
		}
	}
	fclose(fp);
	fp=fopen(argv[4], "w");
	if(!fp)
	{
		printf("Can not open output file.\n");
		clean(-1);
	}
	float *pwr=(float*)malloc(gens*sizeof(float));
	if(!pwr)
	{
		printf("Insufficient memory.\n");
		clean(-1);
	}
	for(i=0;i<loads;i++)
	{
		if(i!=0) fputc('\n', fp);
		ecodispatch(state[i], load[i], pwr);
		for(j=0;j<gens;j++)
		{
			if(j!=0) fputc(',', fp);
			fprintf(fp, "%.2f", pwr[j]);
		}
	}
	fclose(fp);
	free((void*)pwr);
	clean(0);
	return 0;
}
