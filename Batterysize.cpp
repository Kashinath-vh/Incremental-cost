#include<iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdlib>
using namespace std;

class gen{
public:
int isDG,qlim,type;
double mp,nq,V0,w0,Qmax,a,b,c,Pmax,Smax,Pg,Qg;
gen();
};

gen::gen ()
{
isDG=0;
qlim=0;
V0=1.0;
w0=1.0;
Qmax=10000000;
Smax=10000000;
type=1;
}

class load{
public:
int type;
double Pl0,Ql0,alpha,beta,kpf,kqf;
load();
};

load::load()
{
type=0;
Pl0=0;
Ql0=0;
alpha=0;
kpf=0;
kqf=0;
beta=0;
}

class line{
public:
int isON,from,to;
double r,x;
line();
};

line::line()
{
isON=0;
}

double w0,PI=3.141592654,w=1.0,V0=1.0,Ploss,tolerance,base_V,base_VA,loaddata[101],J[250][250],J1[250][250],battpower[101],
ncharge=1.0,ndis=1.0,del;
int no_lines,no_tie,no_iterations,no_buses,k1,k2,no_gens,no_DG,m,n;

gen DG[120];
load LD[120];
line tran[140];

complex<double> Scal[120],V[120],Z_base,SLtotal,SL[120];

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
x=(1.0/del +50.0-double(i))/(1.0/del);
return(x);
}

int soc2i(double x)
{
return(1.0/del +50.0 -round((1.0/del)*x));
}

complex<double> sum(complex<double> * x,int n)
{
int i;
complex<double> y;
for(i=1;i<=n;i++)
y=y + x[i];
return y;
}

double absmax(double* x,int n)
{int i;
double y;
y=abs(x[1]);
for(i=2;i<=n;i++)
if(abs(x[i])>y)
y=abs(x[i]);
return y;
}

int invshipley(double a[250][250],int n )
{
int i,j,k,precision_no=6;
for(k=1;k<=n;k++)
{
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
{
if(i==k || j==k)
continue;
a[i][j]=a[i][j] - (a[i][k]*a[k][j])/a[k][k];
}
for(i=1;i<=n;i++)
{if (i==k)
continue;
a[i][k]=(-1*a[i][k])/a[k][k];
}

for(j=1;j<=n;j++)
{if (j==k)
continue;
a[k][j]=(-1.0*a[k][j])/a[k][k];
}

a[k][k]=(-1.0/a[k][k]);
}
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
a[i][j]=-1.0*a[i][j];
return 0;
}

int islandnrlf(double lpercent, double Pdg1=0)
{
int i,j,k,l,from_bus,to_bus;
double delP[250],delV[250],Vang,Vmag,r,x;
complex<double> I[120],SG[120],z,Y[120][120];
w=1.0;
V[1]=1.0 + 0i;
for(i=2;i<=no_buses;i++)
V[i]= 1. + 0i;
for(i=1;i<=no_buses;i++)
SG[i]=complex<double>(0,0);

for(k=1;k<=no_iterations;k++)
{
//Ybus
for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
Y[i][j]=complex<double>(0,0);

for(i=1;i<=(no_lines+no_tie);i++)
if(tran[i].isON==1)
{
from_bus=tran[i].from;
to_bus=tran[i].to;
r=tran[i].r;
x=tran[i].x;
z=complex<double>(r,w*x);
Y[from_bus][to_bus]=Y[from_bus][to_bus] - Z_base/z;
Y[from_bus][from_bus]=Y[from_bus][from_bus] + Z_base/z;
Y[to_bus][to_bus]=Y[to_bus][to_bus] + Z_base/z;
Y[to_bus][from_bus]=Y[from_bus][to_bus];
}

for(i=1;i<=no_buses;i++)
if(DG[i].isDG)
{
if(DG[i].type==1)
{
r=(1/DG[i].mp)*(DG[i].w0-w);
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==2)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==3)
{
double a,b;
a=(1/DG[i].mp)*(DG[i].w0-w);
b=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
r=(a+b)/2.0;
x=(b-a)/2.0;
}
if(DG[i].type==4)
{
r=DG[i].Pg;
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==5)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=DG[i].Pg;
}
DG[i].Qmax=sqrt(DG[i].Smax*DG[i].Smax - r*r);
DG[i].qlim=0;
if(x>DG[i].Qmax)
{
x=DG[i].Qmax;
DG[i].qlim=1;
}

if(x<-DG[i].Qmax)
{
x=-DG[i].Qmax;
DG[i].qlim=1;
}
DG[i].Qg=x;
SG[i]=complex<double>(r,x);
}


for(i=1;i<=no_buses;i++)
{
r=LD[i].Pl0*(pow(abs(V[i]),LD[i].alpha))*(1.0+LD[i].kpf*(w-1.0));
x=LD[i].Ql0*(pow(abs(V[i]),LD[i].beta))*(1.0+LD[i].kqf*(w-1.0));
SL[i]=complex<double>(r,x)*lpercent;
}
if(Pdg1>0)
SL[k1]+=complex<double>(Pdg1/ncharge,0);
else
SL[k1]+=complex<double>(Pdg1*ndis,0);
for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
J[i][j]=0;

for(i=1;i<=no_buses;i++)
{complex<double> temp=0;
for(j=1;j<=no_buses;j++)
temp+=Y[i][j]*V[j];
I[i]=temp;
}

for(i=1;i<=no_buses;i++)
Scal[i]=(V[i]*conj(I[i]));

//    J11
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
{
if(i==j)
J[i-1][j-1]=-(imag(Scal[i]) + abs(V[i])*abs(V[i])*imag(Y[i][i]));
else
J[i-1][j-1]=-(abs(V[i])*abs(V[j])*abs(Y[i][j])*sin( arg(Y[i][j]) + arg(V[j]) - arg(V[i])));
}
//J21
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
J[i-2+no_buses][j-1]= real(Scal[i]) - abs(V[i])*abs(V[i])*real(Y[i][i]);
else
J[i-2+no_buses][j-1]= -abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J12
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
{
J[i-1][j-2+no_buses]= real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]) + real(SL[i])*LD[i].alpha;
if(DG[i].isDG)
{
if(DG[i].type==2)
J[i-1][j-2+no_buses]+=abs(V[i])/DG[i].nq;
if(DG[i].type==3)
J[i-1][j-2+no_buses]+=0.5*abs(V[i])/DG[i].nq;
}
}
else
J[i-1][j-2+no_buses]=abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J22
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
{
J[i-2+no_buses][j-2+no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]) + LD[i].beta*imag(SL[i]);
if(DG[i].isDG&&(DG[i].qlim==0))
{
if(DG[i].type==1)
J[i-2+no_buses][j-2+no_buses]+=(1/DG[i].nq)*abs(V[i]);
if(DG[i].type==3)
J[i-2+no_buses][j-2+no_buses]+=(0.5/DG[i].nq)*abs(V[i]);
if(DG[i].type==4)
J[i-2+no_buses][j-2+no_buses]+=(1/DG[i].nq)*abs(V[i]);
}

}
else
J[i-2+no_buses][j-2+no_buses]=-abs(V[i])*abs(V[j])*abs(Y[i][j])*sin(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J13
j=2*no_buses - 1;
for(i=2;i<=no_buses;i++)
{
complex<double> temp=0+0i;
for(l=1;l<=no_buses;l++)
{
if(abs(Y[i][l]))
{
r=real(complex<double>(-1.0,0)/Y[i][l]);
temp+=(V[i]-V[l])*(Y[i][l]/w)*(complex<double>(1,0) + Y[i][l]*complex<double>(r,0))*conj(V[i]);
}
}
complex<double> Yii=0+0i;
int z;
for(z=1;z<=no_buses;z++)
Yii+=Y[i][z];
if(abs(Yii))
{
r=real(complex<double>(1.0,0)/Yii);
temp-=V[i]*(Yii/w)*(complex<double>(1,0)+complex<double>(r,0)*Yii)*conj(V[i]);
}
J[i-1][j]=real(temp) + LD[i].Pl0*pow(abs(V[i]),LD[i].alpha)*LD[i].kpf*w;
J[i-2+no_buses][j]=-imag(temp) + LD[i].Ql0*pow(abs(V[i]),LD[i].beta)*LD[i].kqf*w;
if(DG[i].isDG)
{
if(DG[i].type==1)
J[i-1][j]+=(1/DG[i].mp);
if((DG[i].type==2)&&(DG[i].qlim==0))
J[i-2+no_buses][j]-=(1/DG[i].mp);
if(DG[i].type==3)
{
J[i-1][j]+=(0.5/DG[i].mp);
J[i-2+no_buses][j]-=(0.5/DG[i].mp);
}
if((DG[i].type==5)&&(DG[i].qlim==0))
J[i-2+no_buses][j]-=(1/DG[i].mp);
}
}

//J14
j=2*no_buses;
for(i=2;i<=no_buses;i++)
{
J[i-1][j]=abs(V[i])*abs(V[1])*abs(Y[i][1])*cos(arg(Y[i][1]) + arg(V[1]) - arg(V[i]));
}

//J24
j=2*no_buses;
for(i=2;i<=no_buses;i++)
{
J[i-2+no_buses][j]=-abs(V[i])*abs(V[1])*abs(Y[i][1])*sin(arg(Y[i][1]) + arg(V[1]) - arg(V[i]));
}

//J31
for(j=2;j<=no_buses;j++)
J[2*no_buses - 1][j-1]=-(abs(V[1])*abs(V[j])*abs(Y[1][j])*sin( arg(Y[1][j]) + arg(V[j]) - arg(V[1])));

//J32
for(j=2;j<=no_buses;j++)
{
J[2*no_buses - 1][j-2+no_buses]=abs(V[1])*abs(V[j])*abs(Y[1][j])*cos(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J41
i=2*no_buses;
for(j=2;j<=no_buses;j++)
{
J[i][j-1]=-abs(V[1])*abs(V[j])*abs(Y[1][j])*cos(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J42
i=2*no_buses;
for(j=2;j<=no_buses;j++)
{
J[i][j-2+no_buses]=-abs(V[1])*abs(V[j])*abs(Y[1][j])*sin(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J33 & J43
i=1;
{
complex<double> temp=0+0i;
for(l=1;l<=no_buses;l++)
{
if(abs(Y[i][l]))
{
r=real(complex<double>(-1.0,0)/Y[i][l]);
temp+=(V[i]-V[l])*(Y[i][l]/w)*(complex<double>(1,0) + Y[i][l]*complex<double>(r,0))*conj(V[i]);
}
}
complex<double> Yii=0+0i;
int z;
for(z=1;z<=no_buses;z++)
Yii+=Y[i][z];
if(abs(Yii))
{
r=real(complex<double>(1.0,0)/Yii);
temp-=V[i]*(Yii/w)*(complex<double>(1,0)+complex<double>(r,0)*Yii)*conj(V[i]);
}
J[2*no_buses - 1][2*no_buses - 1]=real(temp) + LD[i].Pl0*pow(abs(V[i]),LD[i].alpha)*LD[i].kpf*w;
J[2*no_buses][2*no_buses - 1]=-imag(temp) + LD[i].Ql0*pow(abs(V[i]),LD[i].beta)*LD[i].kqf*w;
if(DG[i].isDG)
{
if(DG[i].type==1)
J[2*no_buses - 1][2*no_buses - 1]+=(1/DG[i].mp);
if((DG[i].type==2)&&(DG[i].qlim==0))
J[2*no_buses][2*no_buses - 1]-=(1/DG[i].mp);
if(DG[i].type==3)
{
J[2*no_buses - 1][2*no_buses - 1]+=(0.5/DG[i].mp);
J[2*no_buses][2*no_buses - 1]-=(0.5/DG[i].mp);
}
if((DG[i].type==5)&&(DG[i].qlim==0))
J[2*no_buses][2*no_buses - 1]-=(1/DG[i].mp);
}
}

//J34
i=1;
J[2*no_buses - 1][2*no_buses]=real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]) + real(SL[i])*LD[i].alpha;
if(DG[i].isDG)
{
if(DG[i].type==2)
J[2*no_buses - 1][2*no_buses]+=abs(V[i])/DG[i].nq;
if(DG[i].type==3)
J[2*no_buses - 1][2*no_buses]+=0.5*abs(V[i])/DG[i].nq;
}

//J44
i=1;
j=1;
J[2*no_buses][2*no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]);
if(DG[i].isDG&&(DG[i].qlim==0))
{
if(DG[i].type==1)
J[2*no_buses][2*no_buses]+=(1/DG[i].nq)*abs(V[i]);
if(DG[i].type==3)
J[2*no_buses][2*no_buses]+=(0.5/DG[i].nq)*abs(V[i]);
if(DG[i].type==4)
J[2*no_buses][2*no_buses]+=(1/DG[i].nq)*abs(V[i]);
}

for(i=1;i<=2*no_buses;i++)
for(j=1;j<=2*no_buses;j++)
J1[i][j]=J[i][j];

invshipley(J,2*no_buses);

for(i=2;i<=no_buses;i++)
delP[i-1]=-(real(Scal[i]) + real(SL[i]) -real(SG[i]));
for(i=no_buses;i<=2*(no_buses - 1);i++)
delP[i]=-(imag(Scal[i-no_buses+2]) + imag(SL[i-no_buses+2]) - imag(SG[i-no_buses+2]));


i=1;
delP[2*no_buses - 1]=-(real(Scal[i]) + real(SL[i]) -real(SG[i]));
delP[2*no_buses]=-(imag(Scal[i]) + imag(SL[i]) - imag(SG[i]));

for(i=1;i<=2*(no_buses);i++)
{
double temp=0;
for(j=1;j<=2*(no_buses);j++)
temp+=J[i][j]*delP[j];
delV[i]=temp;
}

for(i=1;i<=(no_buses - 1);i++)
{
Vang=arg(V[i+1])+delV[i];
Vmag=abs(V[i+1])*(1.0 + delV[i+ no_buses - 1]);
V[i+1]=polar(Vmag,Vang);
}

w+=delV[2*no_buses - 1];
Vmag=abs(V[1])*(1.0 + delV[2*no_buses]);
V[1]=polar(Vmag,0.0);

if(tolerance>absmax(delP,(2*no_buses)))  
break;
}


if(k<=no_iterations)
{
SLtotal=sum(Scal,no_buses);
Ploss=real(SLtotal)*base_VA/1000;
for(i=1;i<=no_buses;i++)
if(DG[i].isDG)
{
if(DG[i].type==1)
{
r=(1/DG[i].mp)*(DG[i].w0-w);
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==2)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==3)
{
double a,b;
a=(1/DG[i].mp)*(DG[i].w0-w);
b=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
r=(a+b)/2.0;
x=(b-a)/2.0;
}
if(DG[i].type==4)
{
r=DG[i].Pg;
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==5)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=DG[i].Pg;
}
DG[i].qlim=0;
if(x>DG[i].Qmax)
{
x=DG[i].Qmax;
DG[i].qlim=1;
}

if(x<-DG[i].Qmax)
{
x=-DG[i].Qmax;
DG[i].qlim=1;
}
DG[i].Pg=r;
DG[i].Qg=x;
}
return 1;
}
else
cout<<"\nLoad-flow solution did not converge\n";
return 0;
}

int main()
{
int i,j,k,kmin,kmax,l,n1,count,no_hrs=24,no_days=4,dsocmax,ndel,temp,deli,nbatt=0,LBN,e,flag;
double r,x,price[101],cost[320][103],socmin[101],socmin1[101],socmin2[101],socmax[101],socmax1[101],socmax2[101],socst,smin,smax,temp1,soc,
socT[102],charge[102],Ps[101][7],pvsize,Ppv,Pd,fgold,lambda,Pdemand,Wmin,Ws,mpmax,Pmax,Emax,Cbatt,cost1,Jlbn[250],alpha[250],alpha1,penalty,
Plold,Pcontrol,tcostold,delbat,a[200],tstep,ntstep; 
complex<double> z,line_z[120];
string comment;
fstream fp;

fp.open("islandip33.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
fp>>no_buses;
fp>>no_lines;
fp>>tolerance;
fp>>no_iterations;
fp>>base_V;
fp>>base_VA;
fp>>no_DG;
fp>>no_tie;
getline(fp,comment);
getline(fp,comment);
Z_base=(base_V*base_V)/base_VA;
for(i=1;i<=no_lines;i++)
{
fp>>temp;
tran[i].isON=1;
fp>>tran[i].from;
fp>>tran[i].to;
fp>>tran[i].r;
fp>>tran[i].x;
z=complex<double>(tran[i].r,tran[i].x);
line_z[i]=z;
}
getline(fp,comment);
getline(fp,comment);

for(i=1;i<=no_buses;i++)
{
fp>>temp;
fp>>r;
fp>>x;
fp>>temp;
if(temp==1)
{
LD[i].alpha=1.51;
LD[i].beta=3.4;
}
if(temp==2)
{
LD[i].alpha=0.18;
LD[i].beta=6.0;
}
if(temp==3)
{
LD[i].alpha=0.92;
LD[i].beta=4.04;
}
LD[i].Pl0=r*1000/base_VA;
LD[i].Ql0=x*1000/base_VA;
}

getline(fp,comment);
getline(fp,comment);
for(i=1;i<=no_DG;i++)
{
fp>>temp;
fp>>r;
DG[temp].mp=r;
fp>>x;
DG[temp].nq=x;
fp>>x;
DG[temp].Qmax=x;
DG[temp].isDG=1;
fp>>r;
DG[temp].V0=r;
fp>>DG[temp].type;
}
getline(fp,comment);
getline(fp,comment);
for(i=no_lines+1;i<=(no_lines+no_tie);i++)
{
fp>>tran[i].from;
fp>>tran[i].to;
fp>>tran[i].r;
fp>>tran[i].x;
}
}
else
{cout<<"Unable to open file\n";
return 0;
}
fp.close();

tstep=1.0/1.0;
ntstep=1.0/tstep;
fp.open("loaddata.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
l=1;
for(i=1;i<=no_hrs;i++)
{
fp>>r;
for(j=1;j<=ntstep;j++)
{
loaddata[l]=r/100.0;
l++;
}
}
loaddata[int(no_hrs*ntstep+1)]=loaddata[1];
}

else
cout<<"Unable to open file\n";

fp.close();
fp.open("nrlfop.txt",ios::out);
fp<<std::fixed<<std::setprecision(4);

k1=12;
k2=1;
Pdemand=0;
for(i=1;i<=no_buses;i++)
Pdemand+=LD[i].Pl0;

DG[6].a=0.3312;
DG[6].b=0.0156;
DG[6].c=0.0002484;
DG[6].Pmax=(1500.0*1000.0)/base_VA;
DG[6].Smax=(1500.0*1000.0)/base_VA;

DG[k2].a=0.4969*(2500.0/1790.0);
DG[k2].b=0.0116;
DG[k2].c=0.0001987*(1790.0/2500.0);
DG[k2].Pmax=(2500.0*1000.0)/base_VA;

ncharge=sqrt(0.9);
ndis=ncharge;
del=0.005;
delbat=100;
Pmax=0.2042;
Emax=2.042;
Cbatt=0.3248;
socst=i2soc(soc2i(0.7));
dsocmax=2;
ndel=int(1.0/del)+1;
kmin=(Pmax*tstep)/(Emax*del*ndis);
kmax=(Pmax*tstep*ncharge)/(Emax*del);
cout<<"kmin="<<kmin<<"\tkmax="<<kmax<<"\n";
smin=0.4;
smax=1.0;
fp<<"socst="<<socst<<"\tM="<<(1.0/del)<<"\n";
socmin[1]=socst;
socmin1[1]=socst;
socmin2[1]=socst;
socmax[1]=socst;
socmax1[1]=socst;
socmax2[1]=socst;

socmin[int(no_hrs*ntstep+1)]=socst;
socmin1[int(no_hrs*ntstep+1)]=socst;
socmin2[int(no_hrs*ntstep+1)]=socst;
socmax[int(no_hrs*ntstep+1)]=socst;
socmax1[int(no_hrs*ntstep+1)]=socst;
socmax2[int(no_hrs*ntstep+1)]=socst;
for(i=0;i<=(1.0/del +50.0);i++)
for(j=0;j<=(no_hrs*ntstep+3);j++)
cost[i][j]=-999999999;
cost[soc2i(socst)][1]=0;
for(i=1;i<=(no_hrs*ntstep);i++)
{
socmin1[i+1]=max((socmin1[i]-kmin*del),smin);
socmin2[int(no_hrs*ntstep-i+1)]=max((socmin2[int(no_hrs*ntstep-i+2)]-kmax*del),smin);

socmax1[i+1]=min((socmax1[i]+kmax*del),smax);
socmax2[int(no_hrs*ntstep-i+1)]=min((socmax2[int(no_hrs*ntstep-i+2)]+kmin*del),smax);
}

for(i=1;i<=(no_hrs*ntstep);i++)
{
socmin[i+1]=max(socmin1[i+1],socmin2[i+1]);
socmax[i+1]=min(socmax1[i+1],socmax2[i+1]);
}
for(i=1;i<=(no_hrs*ntstep+1);i++)
cout<<socmin[i]<<"\t"<<socmax[i]<<"\n";
tcostold=9999999999999;
if(socst<smin)
{
cout<<"Incorrect data\n";
return 0;
}
n1=soc2i(0.0);
k=1;
no_DG=2;
DG[k2].isDG=1;
DG[6].type=5;
DG[k2].type=5;
for(k=1;k<=5000;k++)
{
if(nbatt==3600)
delbat=1;
cout<<nbatt<<"\n";
cost1=0;
m=10;
flag=0;
if(k==1)
{
fp<<"OPF without battery\n";
fp<<"lambda"<<"\t"<<"DG[6].Pg"<<"\t"<<"DG[6].Qg"<<"\t"<<"DG[6].Sg"<<"\t"<<"DG[1].Pg"<<"\t"<<"DG[1].Qg"<<"\t"<<"DG[1].Sg"<<"\t"<<"penalty"<<"\t"<<"Ploss"<<"\n";
}
for(m=1;m<=no_hrs*ntstep;m++)
{
n=2;
Pd=Pdemand*loaddata[m];
Pd+=battpower[m];
lambda=0.05;
penalty=1.0;
Plold=0;
e=0;
DG[k2].nq=DG[6].nq;
DG[k2].mp=DG[6].mp;
do
{
DG[6].Pg=(lambda - DG[6].b)/(2.0*DG[6].c);
DG[k2].Pg=(lambda - DG[k2].b)/(2.0*DG[k2].c);
lambda+=((Pd*base_VA/1000.0) - (DG[6].Pg + DG[k2].Pg))/(1.0/(2.0*DG[k2].c) + 1.0/(2.0*DG[6].c));
}while(abs((Pd*base_VA/1000.0) - (DG[6].Pg + DG[k2].Pg))>(Pd*base_VA/1000000.0));
LBN=6;
DG[LBN].type=2;
if(LBN==1)
LBN=2*no_buses;
do
{
DG[6].Pg=(lambda/penalty - DG[6].b)/(2.0*DG[6].c);
DG[k2].Pg=(lambda - DG[k2].b)/(2.0*DG[k2].c);
DG[6].Pg*=(1000.0/base_VA);
DG[k2].Pg*=(1000.0/base_VA);
if(islandnrlf(loaddata[m],battpower[m]))
flag=1;
else
{
flag=0;
}
lambda+=0.0002*(-(lambda/penalty - DG[6].b)/(2.0*DG[6].c)-(lambda - DG[k2].b)/(2.0*DG[k2].c)+Ploss+Pd*base_VA/1000.0);
for(i=1;i<=(2*no_buses);i++)
Jlbn[i]=J1[LBN-1][i];
for(i=LBN;i<=(2*no_buses);i++)
for(j=1;j<=(2*no_buses);j++)
J1[i-1][j]=J1[i][j];
invshipley(J1,(2*no_buses-1));
for(i=1;i<=(2*no_buses-1);i++)
{
temp1=0;
for(j=1;j<=(2*no_buses-1);j++)
temp1+=Jlbn[j]*J1[j][i];
alpha[i]=temp1;
}
penalty=(-1.0/alpha[5]);
}while(abs((Pd*base_VA/1000.0) + Ploss - ((lambda/penalty - DG[6].b)/(2.0*DG[6].c) + (lambda - DG[k2].b)/(2.0*DG[k2].c)))>(Pd*base_VA/100000000.0)&&(flag==1));
if(k==1)
fp<<lambda<<"\t"<<DG[6].Pg*base_VA/1000.0<<"\t"<<DG[6].Qg*base_VA/1000.0<<"\t"<<sqrt(DG[6].Qg*DG[6].Qg+DG[6].Pg*DG[6].Pg)*base_VA/1000.0<<"\t"<<DG[k2].Pg*base_VA/1000.0<<"\t"<<DG[k2].Qg*base_VA/1000.0<<"\t"<<sqrt(DG[k2].Qg*DG[k2].Qg+DG[k2].Pg*DG[k2].Pg)*base_VA/1000.0<<"\t"<<penalty<<"\t"<<Ploss<<"\n";
price[m]=lambda;
cost1+=(DG[6].a+DG[6].b*DG[6].Pg*base_VA/1000.0 + DG[6].c*(DG[6].Pg*base_VA/1000.0)*(DG[6].Pg*base_VA/1000.0))*tstep;
cost1+=(DG[k2].a+DG[k2].b*DG[k2].Pg*base_VA/1000.0 + DG[k2].c*(DG[k2].Pg*base_VA/1000.0)*(DG[k2].Pg*base_VA/1000.0))*tstep;
}
if(k==1)
{
fp<<"Total DG cost="<<cost1<<"$/day\n";
fp<<"Total battery cost="<<(nbatt*Cbatt)<<"$/day\n";
fp<<"Total cost="<<(cost1 +nbatt*Cbatt)<<"$/day\n";
}
for(j=2;j<=(no_hrs*ntstep+1);j++)
for(i=soc2i(1.0);i<=n1;i++)
{
if((i2soc(i)-socmax[j]>0.000001)||(socmin[j]-i2soc(i)>0.000001))
{
cost[i][j]=-999999999;
}
else
{
l=1;
for(count=-kmin;count<=kmax;count++)
{
if(count<0)
a[l]=cost[i+count][j-1]-del*count*Emax*ndis*price[j-1];
else
a[l]=cost[i+count][j-1]-del*count*Emax*price[j-1]/ncharge;
l++;
}
cost[i][j]=max(a,(kmin+kmax+1));
//cout<<cost[i][j]<<"\n";
}}
i=soc2i(socst);
deli=0;
temp1=0;
soc=socst;
charge[int(no_hrs*ntstep+1)]=0;
//fp<<"Savings="<<cost[soc2i(socst)][int(no_hrs*ntstep+1)]<<"\n\n";
for(j=(no_hrs*ntstep+1);j>=1;j--)
{
l=1;
for(count=-kmin;count<=kmax;count++)
{
if((i+count)<0)
{
a[l]=-999999999;
l++;
continue;
}
if(count<0)
{
if(cost[i][j]==(cost[i+count][j-1] -del*count*Emax*ndis*price[j-1]))
{
charge[j]=count;
break;
}
}
else
if(cost[i][j]==(cost[i+count][j-1] -del*count*Emax*price[j-1]/ncharge))
{
charge[j]=count;
break;
}
l++;
}
i+=charge[j];
//fp<<j<<"\t"<<(cost[soc2i(socst)][int(no_hrs*ntstep+1)]-temp1)<<"\t\t"<<charge[j]<<"\t\t"<<price[j-1]<<"\t\t"<<soc<<"\n";
if(charge[j]>0)
temp1-=charge[j]*price[j-1]*del*Emax/ncharge;
else
temp1-=charge[j]*price[j-1]*del*Emax*ndis;
soc-=charge[j]*del;
}
if(cost[soc2i(socst)][int(no_hrs*ntstep+1)]>(Cbatt))
{
tcostold=(cost1 +nbatt*Cbatt);
nbatt+=delbat;
for(m=1;m<=(no_hrs*ntstep);m++)
{
if(charge[m+1]>0)
battpower[m]+=(delbat*charge[m+1]*del*Emax*1000.0/(base_VA*ncharge))/tstep;
else
battpower[m]+=(ndis*delbat*charge[m+1]*del*Emax*1000.0/(base_VA))/tstep;
}
}
else
{
cout<<"stopped after "<<k<<"\n";
break;
}
}
socT[1]=socst;
for(m=1;m<=(no_hrs*ntstep);m++)
{
if(battpower[m]>0)
socT[m+1]=socT[m]+ncharge*battpower[m]*tstep*base_VA/(1000.0*nbatt*Emax);
else
{
socT[m+1]=socT[m]+battpower[m]*tstep*base_VA/(1000.0*nbatt*Emax*ndis);
}
}
fp<<"\n\nOptimal # of batteries="<<nbatt<<"\n"<<"nsys="<<(ncharge*ndis)<<"\tsocmin="<<smin<<"\tCbatt="<<Cbatt<<"\n";
fp<<"Total DG cost="<<cost1<<"$/day\n";
fp<<"Total battery cost="<<(nbatt*Cbatt)<<"$/day\n";
fp<<"Total cost="<<(cost1 +nbatt*Cbatt)<<"$/day\n";
fp<<"\nOptimal battery power (Pmax="<<Pmax<<" kW Emax="<<Emax<<" kWh)\n";
fp<<"hour\tpower(kW)\tpower per battery\t\tin pu\t\tSOC\n";
for(m=1;m<=(no_hrs*ntstep+1);m++)
{
fp<<m<<"\t"<<battpower[m]*base_VA/1000.0<<"\t\t"<<battpower[m]*base_VA/(1000.0*nbatt)<<"\t\t\t"<<battpower[m]*base_VA/(1000.0*nbatt*Pmax)<<"\t\t"<<socT[m]<<"\n";
}
fp.close();
return 0;
}
