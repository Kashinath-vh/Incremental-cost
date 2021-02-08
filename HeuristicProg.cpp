#include<iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdlib>
using namespace std;
double pavg,socmin[25],socmax[25],socst,price[25];
int no_hrs=24;

double max(double* a, int n)
{
int i;
double x;
x=a[1];
for(i=2;i<=n;i++)
if(a[i]>x)
x=a[i];
return x;
}

double min(double* a, int n)
{
int i;
double x;
x=a[1];
for(i=2;i<=n;i++)
if(a[i]<x)
x=a[i];
return x;
}

int main()
{
int i,j,n,k;
double r,smin,smax,x[26],xg[26],data[1000],delx,ncharge=1.0,ndis=1.0,Pmax,r1,savings,global=0,sum=0,best=0,
Emax,cumax,cumin;
string comment;
fstream fp;
fp.open("price1.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
for(i=1;i<=no_hrs;i++)
{
fp>>r;
price[i]=r;
}
}
else
cout<<"Unable to open file\n";
fp.close();

fp.open("nrlfop.txt",ios::out);
fp<<std::fixed<<std::setprecision(4);
socst=0.5;					//SoC at the start, %
ncharge=sqrt(1.0);				//Charging effeiciency, %
ndis=ncharge;					//Discharging effeiciency, %
Pmax=1;						//Power limit of the battery, kW
Emax=10.0;					//Storage capacity, kWh
smin=0.2;					//Minimum limit of SoC, %
smax=1.0;					//Maximum limit of SoC, %
socst=0.5;
socmin[1]=socst;
socmax[1]=socst;
socmin[25]=socst;
socmax[25]=socst;
for(i=1;i<=(no_hrs/2);i++)
{
socmin[i+1]=max((socmin[i]-Pmax/(Emax*ndis)),smin);
socmin[no_hrs-i+1]=max((socmin[no_hrs-i+2]-Pmax*ncharge/Emax),smin);

socmax[i+1]=min((socmax[i]+Pmax*ncharge/Emax),smax);
socmax[no_hrs-i+1]=min((socmax[no_hrs-i+2]+Pmax/(Emax*ndis)),smax);
}
for(i=1;i<=no_hrs;i++)
sum+=price[i];
pavg=sum/24.0;
cumax=max(price,no_hrs);
cumin=min(price,no_hrs);
srand (time(NULL));				//Reset the randomn number generator
k=1;
fp<<"Ieration\tRevenue\n";
for(k=1;k<=100;k++)
{
x[1]=socst;
global=0;
for(j=1;j<=100;j++)
{
for(i=2;i<=(no_hrs+1);i++)
{
//				Normalizing with (cumax-cumin) gives better results
//r1=((double) rand() / (RAND_MAX)) - 0.5 - (price[i-1]-pavg)/(cumax-cumin);
r1=((double) rand() / (RAND_MAX)) - 0.5 - (price[i-1]-pavg)/pavg;
delx=r1*(Pmax/Emax)*2.0*8.0;
if(delx>(Pmax*ncharge/Emax))
delx=(Pmax*ncharge/Emax);
if(delx<-(Pmax/(Emax*ndis)))
delx=-(Pmax/(Emax*ndis));
x[i]=x[i-1]+delx;
if(x[i]>socmax[i])
{
x[i]=socmax[i];
}
if(x[i]<socmin[i])
{
x[i]=socmin[i];
}
}
savings=0;
for(i=1;i<=(no_hrs);i++)
if(x[i]>x[i+1])
savings+=(x[i]-x[i+1])*Emax*ndis*price[i];
else
savings+=(x[i]-x[i+1])*Emax*price[i]/ncharge;
if(global<savings)
global=savings;
if(best<savings)
{
best=savings;
for(i=1;i<=(no_hrs+1);i++)
xg[i]=x[i];
}
}
data[k]=global;
fp<<k<<"\t\t"<<global<<"\n";
}
fp<<"\nBest solution = "<<best<<"\n";
savings=0;
fp<<"\tSoC\tSavings\n";
for(i=1;i<=(no_hrs+1);i++)
{
fp<<"hour "<<i<<"\t"<<xg[i]<<"\t"<<savings<<"\n";
if(xg[i]>xg[i+1])
savings+=(xg[i]-xg[i+1])*Emax*ndis*price[i];
else
savings+=(xg[i]-xg[i+1])*Emax*price[i]/ncharge;
}

fp.close();
return 0;
}
