#include<iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
using namespace std;

int M;
double del;
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

double i2soc(int i)
{
double x=0;
x=(M +50.0-double(i))*del;
return(x);
}

int soc2i(double x)
{
return(M +50.0 -round(M*x));
}

int main()
{
int i,j,k,kmin,kmax,l,m,n,no_hrs=24,dsocmax,ndel;
double r,price[101],cost[320][103],socmin[101],socmin1[101],socmin2[101],socmax[101],socmax1[101],socmax2[101],socst,smin,smax,temp,soc,
charge[102],ncharge=1.0,ndis=1.0,Pmax=1.0,Emax=10.0,a[500],battpower[101];
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
cout<<std::fixed<<std::setprecision(4);
M=200;						//No. of discrete levels (The example in paper assumed M=10)
del=1.0/double(M);
socst=i2soc(soc2i(0.5));			//SoC at the start, %
Pmax=1;						//Power limit of the battery, kW
Emax=10.0;					//Storage capacity, kWh
dsocmax=2;
ndel=M+1;
smin=0.2;					//Minimum limit of SoC, %
smax=1.0;					//Maximum limit of SoC, %
ncharge=sqrt(1.0);				//Charging effeiciency, %
ndis=ncharge;					//Discharging effeiciency, %
kmin=(Pmax*M)/(Emax*ndis);
kmax=(Pmax*M*ncharge)/(Emax);
fp<<"kmin="<<kmin<<"\tkmax="<<kmax<<"\n";
fp<<"socst="<<socst;
socmin[1]=socst;
socmin1[1]=socst;
socmin2[1]=socst;
socmax[1]=socst;
socmax1[1]=socst;
socmax2[1]=socst;

socmin[int(no_hrs+1)]=socst;
socmin1[int(no_hrs+1)]=socst;
socmin2[int(no_hrs+1)]=socst;
socmax[int(no_hrs+1)]=socst;
socmax1[int(no_hrs+1)]=socst;
socmax2[int(no_hrs+1)]=socst;
for(i=0;i<=(M +50.0);i++)
for(j=0;j<=(no_hrs+3);j++)
cost[i][j]=-999999999;
cost[soc2i(socst)][1]=0;
for(i=1;i<=(no_hrs);i++)
{
socmin1[i+1]=max((socmin1[i]-kmin*del),smin);
socmin2[int(no_hrs-i+1)]=max((socmin2[int(no_hrs-i+2)]-kmax*del),smin);

socmax1[i+1]=min((socmax1[i]+kmax*del),smax);
socmax2[int(no_hrs-i+1)]=min((socmax2[int(no_hrs-i+2)]+kmin*del),smax);
}
for(i=1;i<=(no_hrs);i++)
{
socmin[i+1]=max(socmin1[i+1],socmin2[i+1]);
socmax[i+1]=min(socmax1[i+1],socmax2[i+1]);
}
for(i=1;i<=(no_hrs+1);i++)
cout<<socmin[i]<<"\t"<<socmax[i]<<"\n";

n=soc2i(0.0);
for(j=2;j<=(no_hrs+1);j++)
for(i=soc2i(1.0);i<=n;i++)
{
if((i2soc(i)-socmax[j]>0.000001)||(socmin[j]-i2soc(i)>0.000001))
{
cost[i][j]=-999999999;
}
else
{
l=1;
for(k=-kmin;k<=kmax;k++)
{
if((i+k)<0)
{
a[l]=-999999999;
l++;
continue;
}
if((i+k)>319)
{
a[l]=-999999999;
l++;
continue;
}
if(k<0)
a[l]=cost[i+k][j-1]-del*k*Emax*ndis*price[j-1];
else
a[l]=cost[i+k][j-1]-del*k*Emax*price[j-1]/ncharge;
l++;
}
cost[i][j]=max(a,(kmin+kmax+1));
}}
temp=0;
soc=socst;
charge[int(no_hrs+1)]=0;
fp<<"\tSavings="<<cost[soc2i(socst)][int(no_hrs+1)]<<"\tn="<<(ndis*ncharge)<<"\tM="<<M<<"\n";
cout<<"Savings="<<cost[soc2i(socst)][int(no_hrs+1)]<<"\n";

fp<<"\nsoc\t";
for(j=1;j<=(no_hrs+1);j++)
fp<<j<<"\t";
fp<<"\n";
for(i=soc2i(1.0);i<=n;i++)
{
fp<<i2soc(i)<<"\t";
for(j=1;j<=(no_hrs+1);j++)
if(cost[i][j]==-999999999)
fp<<"x\t";
else
fp<<cost[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"Hour"<<"\t"<<"Revenue"<<"\t\t"<<"charge"<<"\t\t"<<"price"<<"\t\t"<<"SoC"<<"\n";
temp=0;
soc=socst;
i=soc2i(socst);
j=(no_hrs+1);
for(j=(no_hrs+1);j>=1;j--)
{
l=1;
for(k=-kmin;k<=kmax;k++)
{
if((i+k)<0||(i+k)>319)
continue;
if(k<0)
{
if(cost[i][j]==(cost[i+k][j-1] -del*k*Emax*ndis*price[j-1]))
{
charge[j]=k;
break;
}
}
else
if(cost[i][j]==(cost[i+k][j-1] -del*k*Emax*price[j-1]/ncharge))
{
charge[j]=k;
break;
}
l++;
}
i+=charge[j];
fp<<j<<"\t"<<(cost[soc2i(socst)][int(no_hrs+1)]-temp)<<"\t\t"<<charge[j]<<"\t\t"<<price[j-1]<<"\t\t"<<soc<<"\n";
if(charge[j]>0)
temp-=charge[j]*price[j-1]*del*Emax/ncharge;
else
temp-=charge[j]*price[j-1]*del*Emax*ndis;
soc-=charge[j]*del;
}
for(m=1;m<=(no_hrs);m++)
{
if(charge[m+1]>0)
battpower[m]=(charge[m+1]*del*Emax/(ncharge));
else
battpower[m]=(ndis*charge[m+1]*del*Emax);
}

fp.close();
return 0;
}
