//Dynamic programming UC
//(c) 2015, Bo Gao

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

float dpwr[12];
float gmin[]={25,60};
float gmax[]={250,300};
float gnlc[]={1,1.25};
float ginc[]={81.25,169.75};
//float gavg[]={23.54,20.34,19.74,28.00};
//float gccs[]={350.00,400.00,1100.00,0.00};
//float ghcs[]={150.00,170.00,500.00,0.00};
//int   gctm[]={4,5,5,0};
//int   gton[]={-5,8,8,-6};
//int   gmon[]={4,5,5,1};
//int   gmof[]={2,3,4,1};
int   gpri[]={0,1};
float fmin[]={25,60,85};
float fmax[]={250,300,550};
int   fcfg[]={1,2,3};

int states[12]={0};

/*
int ontime(int time, int gen)
{
	int t, i;
	t=gton[gen];
	t=(t>0)?t:0;
	for(i=0;i<time;i++)
		if((fcfg[states[i]]&(1<<gen))!=0)
			t++;
		else
			t=0;
	return t;
}

int offtime(int time, int gen)
{
	int t, i;
	t=-gton[gen];
	t=(t>0)?t:0;
	for(i=0;i<time;i++)
		if((fcfg[states[i]]&(1<<gen))==0)
			t++;
		else t=0;
	return t;
}
*/

int stage_cost(int time)
{
	int i, j, k;
	int t_on;
	float cost[3];
	int flag;
	if(time==12) return 1;
	for(i=0;i<3;i++)
	{
		if(dpwr[time]>fmax[i] || dpwr[time]<fmin[i]) //Capability check
		{
			cost[i]=FLT_MAX;
			continue;
		}
		/*
		flag=0;
		for(j=0;j<2;j++) //On-off time sanity check
		{
			k=offtime(time, j);
			if((fcfg[i]&(1<<j))!=0 && k!=0 && k<gmof[j]) //Off->On
				flag=1;
			k=ontime(time, j);
			if((fcfg[i]&(1<<j))==0 && k!=0 && k<gmon[j]) //On->Off
				flag=1;
		}
		if(flag)
		{
			cost[i]=FLT_MAX;
			continue;
		}
		*/
		cost[i]=0;
		float pset[2];
		float ptmp=dpwr[time];
		for(j=0;j<2;j++) //Economic dispatch
		{
			if((fcfg[i]&(1<<gpri[j]))!=0)
			{
				if(ptmp>gmax[gpri[j]])
					pset[gpri[j]]=gmax[gpri[j]];
				else
					pset[gpri[j]]=ptmp;
				ptmp-=pset[gpri[j]];
			}
		}
		for(j=0;j<2;j++) //Accumulate cost
			if((fcfg[i]&(1<<j))!=0)
			{
				cost[i]+=gnlc[j]+ginc[j]*pset[j];
				/*
				k=offtime(time, j); //Calculate start cost
				if(k!=0) cost[i]+=ghcs[j];
				if(k!=0 && k>=gctm[j]) cost[i]+=gccs[j]-ghcs[j];
				*/
			}
	}
	float cmin;
	while(1)
	{
		cmin=FLT_MAX;
		for(i=0;i<3;i++) //Find best solution
			if(cost[i]<cmin)
			{
				cmin=cost[i];
				states[time]=i;
			}
		if(cmin==FLT_MAX) return 0; //No feasible choice
		if(stage_cost(time+1)==1) //Search for next stage
			return 1;
		cost[states[time]]=FLT_MAX; //Rule out option
	}
}

int parse(char *s, int i)
{
	float a=0, b=0, c=0;
	sscanf(s, "%f,%f,%f", &a, &b, &c);
	if(a==0 || b==0 || c==0 || a<b+c)
		return -1;
	dpwr[i]=a-b-c;
	printf("%f\n", dpwr[i]);
	return 0;
}

int main()
{
	int i=0;
	FILE *fp;
	size_t l;
	fp=fopen("load.csv", "r");
	if(!fp)
	{
		printf("Error reading input.\n");
		return -1;
	}
	char *line;
	while(getline(&line, &l, fp)!=-1)
	{
		if(line[0]=='#')
			continue;
		if(parse(line, i++))
		{
			printf("Error parsing file.\n");
			return -1;
		}
	}
	if(i!=12)
	{
		printf("Incomplete input file.\n");
		return -1;
	}
	//if(line) free((void*)line);
	fclose(fp);
	if(stage_cost(0)==1)
		for(i=0;i<12;i++)
			printf("%d,%d\n", i+1, fcfg[states[i]]);
	else
		printf("No solution found.\n");
	return 0;
}
