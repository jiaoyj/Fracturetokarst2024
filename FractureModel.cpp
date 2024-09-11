// FractureModel.cpp: implementation of the FractureModel class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"

#include "FractureModel.h"
#include "RandomGenerator.h"
#include <math.h>
#include <stdexcept>
#include <armadillo>
using namespace arma;
//using namespace std;

#include <iostream>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FractureModel::FractureModel()
{

}

FractureModel::~FractureModel()
{

}

//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::inputdata()
{
	int i,n,n1;

	FILE *fp1_input;
	if((fp1_input=fopen("In_DFractureParabenchmark.dat","r"))==NULL)
	{
		std::cout << "Fail to open In_DFracturePara.dat\n";
		return;
	}

	FILE *fp2_input;
	if((fp2_input=fopen("In_RFracturePara.dat","r"))==NULL)
	{
		std::cout << "Fail to open In_RFracturePara.dat\n";
		return;
	}

	FILE *fp3_input;
	if((fp3_input=fopen("In_InfilPara.dat","r"))==NULL)
	{
		std::cout << "Fail to open In_InfilPara.dat\n";
		return;
	}

	// input boundary conditions from fp1_input
	fscanf(fp1_input,"%d%d%d",&nboun1,&nboun2,&nbounf); //nboundf
	nbound = nboun1+nboun2+nbounf;
	//
	for(i=1;i<=nbound;i++) fscanf(fp1_input,"%lf%lf%lf%lf",&xy12[i][1],&xy12[i][2],&xy12[i][3],&xy12[i][4]);
	for(i=1;i<=nbound;i++) fscanf(fp1_input,"%lf%lf",&xbp[i],&ybp[i]);
	for(i=1;i<=nboun1;i++) fscanf(fp1_input,"%lf",&hb1[i]);
	for(i=1;i<=nboun2;i++) fscanf(fp1_input,"%lf",&qb2[i]);
	n = nbound;

	// input determined fractures from fp1_input
	fscanf(fp1_input,"%d",&n1);  //"n1" is the amount of determined fracture
	if(n1 >= 1)
	{
		for(i=1;i<=n1;i++)
		{
			n = n +1;//
			fscanf(fp1_input,"%lf%lf%lf%lf%lf",&xy12[n][1],&xy12[n][2],&xy12[n][3],&xy12[n][4],&b[n]);
		}
	}
	napt = n; //total amount of initial fractures

	fscanf(fp3_input,"%lf%lf",&Qrain,&con_recharg);//
	Qsurfcorro =  Qrain * con_recharg ;

	// input parameters for generating random fractures from fp2_input
	fscanf(fp2_input,"%lf%lf%lf%lf",&minx,&maxx,&miny,&maxy);
	fscanf(fp2_input,"%d%d",&ndx,&ndy);//X and Y
	fscanf(fp2_input,"%d",&ngroup);//
	if(ngroup<1) return; //ngroup<1
	for(i=1;i<=ngroup;i++)
	{
		fscanf(fp2_input,"%d%lf%lf%lf%lf",&ndisl[i],&meanl[i],&sigmal[i],&minl[i],&maxl[i]);
		fscanf(fp2_input,"%d%lf%lf%lf%lf",&ndisb[i],&meanb[i],&sigmab[i],&minb[i],&maxb[i]);
		fscanf(fp2_input,"%d%lf%lf%lf%lf",&ndisd[i],&meand[i],&sigmad[i],&mind[i],&maxd[i]);
		fscanf(fp2_input,"%lf",&meanc[i]); 
	}									

	fclose(fp1_input);
	fclose(fp2_input);
	fclose(fp3_input);


	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
double FractureModel::amin1(double x1,double x2)
{
	if(x1 <= x2) return x1;
	return x2;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
double FractureModel::amax1(double x1,double x2)
{
	if(x1 >= x2) return x1;
	return x2;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::generate_fracture()
{
	int i,j,k,n,kuai,nfikuai; 
	int nkuai;
	double dx,dy,adxdy,d,xo,yo,d1,d2;
	RandomGenerator random_generator;
 
	
	FILE *fp1_output;
	if((fp1_output=fopen("Out_InitFractureNetwork.dat","w"))==NULL)
	{
		std::cout << "fail to open Out_InitFractureNetwork.dat\n";
		return;
	}
	n = napt;//

	if(ngroup<1) goto L1111; //
	
	// generate the random fracture
	dx = (maxx - minx)/double(ndx);
	dy = (maxy - miny)/double(ndy);
	adxdy = dx*dy;
	nkuai = 0;
	for(i=1;i<=ndx;i++)
	{
		for(j=1;j<=ndy;j++)
		{
			nkuai = nkuai+1;
			xynp[nkuai][1] = double(i-1)*dx + minx;
			xynp[nkuai][2] = double(j-1)*dy + miny;
			                                       
		}
	}
	random_generator.init1(0.999);
	for(i=1; i<=ngroup;i++)
	{           //725
		for(kuai=1;kuai<=nkuai;kuai++)
		{ //720
			nfikuai = int(adxdy/meanl[i]/meanc[i]);
			for(k=1;k<=nfikuai;k++)
			{//718
				n = n+1;
				//c      
				//c       width of frature
				//c
				if(ndisb[i] == 1) b[n] = random_generator.unifrm(minb[i],maxb[i]);
				if(ndisb[i] == 2) b[n] = random_generator.gauss(meanb[i],sigmab[i],minb[i],maxb[i]);
				if(ndisb[i] == 3) b[n] = random_generator.lognor(meanb[i],sigmab[i],minb[i],maxb[i]);
				if(ndisb[i] == 4) b[n] = random_generator.erlang(meanb[i],minb[i],maxb[i]);
				//c
				//c       length of fracture
				//c
				if(ndisl[i] == 1) len = random_generator.unifrm(minl[i],maxl[i]);
				if(ndisl[i] == 2) len = random_generator.gauss(meanl[i],sigmal[i],minl[i],maxl[i]);
				if(ndisl[i] == 3) len = random_generator.lognor(meanl[i],sigmal[i],minl[i],maxl[i]);
				if(ndisl[i] == 4) len = random_generator.erlang (meanl[i],minl[i],maxl[i]);
				//c
				//c       diretion of fractures
				//c
				if(ndisd[i] == 1) d = random_generator.unifrm(mind[i],maxd[i]);
				if(ndisd[i] == 2) d = random_generator.gauss(meand[i],sigmad[i],mind[i],maxd[i]);
				if(ndisd[i] == 3) d = random_generator.lognor(meand[i],sigmad[i],mind[i],maxd[i]);
				if(ndisd[i] == 4) d = random_generator.erlang(meand[i],mind[i],maxd[i]);
				//c
				//c       posision of fracture
				//c
				xo = random_generator.unifrm(xynp[kuai][1],xynp[kuai][1]+dx);
				yo = random_generator.unifrm(xynp[kuai][2],xynp[kuai][2]+dy);

				if(d >= 180.0)  d = d - 180.0;
				d2 = cos(d/pii);
				d1 = sin(d/pii);
				xy12[n][3] = xo + 0.5*len*d1;
				xy12[n][1] = xo - 0.5*len*d1;
				xy12[n][4] = yo + 0.5*len*d2;
				xy12[n][2] = yo - 0.5*len*d2;
			}//718     continue
		}//720     continue
	}//725     continue
 
	//finish generating random fracture

L1111:;
	
	// output initial Fracture network(random & determined)-------------------------
	//Out_InitFractureNetwork.dat，
	fprintf(fp1_output,"%d %d %d %d\n",n,nboun1,nboun2,nbound);
	for(i=1;i<=n;i++) fprintf(fp1_output,"%10.4lf%10.4lf%10.4lf%10.4lf\n",xy12[i][1],xy12[i][2],xy12[i][3],xy12[i][4]);
 
	for(i=1;i<=npmax;i++)//#define	npmax 4999 
		{ //
		for(j=1;j<=2;j++)
		{
			xynp[i][j] = 0.0;  
		}
	}

	napt = n; //total amount of initial fractures

	fclose(fp1_output);
	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::poin()
{
	int i,j,ncount,nw;
	double dta,dta1,t1,t2;
	double x,y,x1,x2,y1,y2,a1,a2,b1,b2,f1,f2,f3,f4;

	dta  = 0.002001;
	dta1 = 0.00005;

	for(i=1;i<=napt;i++)
	{// 20 
		ncount=0;
		for(j=1;j<=napt;j++)
		{//10
			if(i <= nbound && j <= nbound) goto L10;
			if(i == j) goto L10;
			x1 = amin1(xy12[i][1],xy12[i][3]) - dta;
			y1 = amin1(xy12[i][2],xy12[i][4]) - dta;
			x2 = amax1(xy12[i][1],xy12[i][3]) + dta;
			y2 = amax1(xy12[i][2],xy12[i][4]) + dta;
			a1 = amin1(xy12[j][1],xy12[j][3]) - dta;
			b1 = amin1(xy12[j][2],xy12[j][4]) - dta;
			a2 = amax1(xy12[j][1],xy12[j][3]) + dta;
			b2 = amax1(xy12[j][2],xy12[j][4]) + dta;
			if(a1 > x2 || a2 < x1 || b2 < y1 || b1 > y2) goto L10;
			f1 = xy12[i][3] - xy12[i][1];
			f2 = xy12[i][4] - xy12[i][2];
			f3 = xy12[j][3] - xy12[j][1];
			f4 = xy12[j][4] - xy12[j][2];
			if(fabs(f1) <= 1.e-10 && fabs(f3) <= 1.e-10) goto L10;
			
			if(fabs(f3) <= 1.e-10)
			{// then
				x = xy12[j][1];
				y = (x - xy12[i][1])*f2/f1 + xy12[i][2];
				goto L7;
			}//end if
			
			if(fabs(f1) <= 1.e-10)
			{// then
				x = xy12[i][1];
				y = (x - xy12[j][1])*f4/f3 + xy12[j][2];
				goto L7;
			}//end if
			
			t2 = f4/f3;
			t1 = f2/f1;
			if(fabs(t1 - t2) <= 1.e-10) goto L10;
			x = (xy12[j][2] + t1*xy12[i][1] - xy12[i][2] - t2*xy12[j][1])/(t1 - t2);
			y=t2*(x-xy12[j][1])+xy12[j][2];

			L7:;//       continue
				
			   if(x >= x1 && x <= x2)
			   {// then
					if(x >= a1 && x <= a2)
					{//then
						if(y >= y1 && y <= y2) 
						{//then
							if(y >= b1 && y <= b2) 
							{//then
								nw =1;
								if(i > nbound && j > nbound) nw = outin(x,y);
								if(nw > 0) 
								{//then
									ncount = ncount + 1;
									xjoin[i][ncount] = x;//
									yjoin[i][ncount] = y;
									ijpoin[i][ncount] = j;//
								}//end if
							}//end if
						}//end if
					}//end if
			   }//end if
			L10:;//  continue
		}//10
		icount[i] = ncount;
	}//20      continue

	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
int FractureModel::outin(double x0,double y0)
{
	int i,i1,npoint,nw;
	double alpha,alpha1,beta,cobeta,alpha2;
	double x1,x2,dx1,dx2,y1,y2,t;

	npoint = nbound;
	alpha = 0.0;
	for(i=1;i<=npoint;i++)
	{//10
		i1 = i+1;
		if(i == npoint) i1=1;
		x1 = xbp[i];
		x2 = xbp[i1];
		y1 = ybp[i];
		y2 = ybp[i1];
		cobeta = (x1-x0)*(x2-x0) + (y2-y0)*(y1-y0);
		dx1 = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
		dx2 = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
		if(dx1 <= 1.e-10 || dx2 <= 1.e-10)
		{// then
			nw = 1;
			return nw;
		}//	end if
		cobeta = cobeta/dx1/dx2;
		if(fabs(cobeta) > 1.0)
		{// then
			if(cobeta > 0.0) cobeta = 0.9999995;
			if(cobeta < 0.0) cobeta = -0.9999995;
		}//end if

		beta = acos(cobeta);
		t = (x2-x0)*(y1-y0) - (y2-y0)*(x1-x0);
		if(t == 0.0)
		{// then
			nw = 1;
			return nw;
		}//	end if

		if(t < 0.0) beta = -beta;
		alpha = alpha + beta;
	}//10      continue

	alpha1 = 1.5*3.14159;
	alpha2 = 0.1;
	if(fabs(alpha) > alpha1) nw = 1;
	if(fabs(alpha) <= alpha2) nw=0;

	return nw;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::connectedFracture()
{
	int i,j,l,ij,j1,napt1,napt2,napt3;

	FILE *fp2_output;
	if((fp2_output=fopen("Out_ConnectedFracNetwork.dat","w"))==NULL)
	{
		std::cout << "fail to open Out_ConnectedNetwork.dat\n";
		return;
	}

	for(i=0;i<=napt;i++)
	{
		icount[i] = 0;
		for(j=0;j<=macoun;j++)
		{
			ijpoin[i][j] = 0;
			xjoin[i][j] = 0.0;
			yjoin[i][j] = 0.0;
		}
	}//This section is added when writting to C++   
	poin();//====================

	nconec = nbound;
	for(i=1;i<=napt;i++) iiw[i] = 2;

	L20:;//         continue

	napt1 = 0;
	for(i=nbound+1;i<=napt;i++) if(iiw[i] == 2) napt1 = napt1+1;//79

	for(i=nbound+1;i<=napt;i++)
	{//89
		j1 = icount[i];
		if(j1 >= macoun) 
		{
			std::cout << "icount is larger than maxcount"; nconec = 0;
			return;
		}
		if(j1 <= 1) iiw[i] = 0;
	}//89         continue

	for(i=nbound+1;i<=napt;i++)
	{//91
		napt3 = 0;
		for(j=nbound+1; j<=macoun;j++)
		{//90
			ij = ijpoin[i][j];
			if(iiw[ij] == 0 && ij >= 1){napt3 = napt3+1; ijpoin[i][j] = 0;}
		}//90  continue
		icount[i] = icount[i] - napt3;
		if(icount[i] <= 1) iiw[i] = 0;
	} //91 continue

	napt2 = 0;
	for(i=nbound+1;i<=napt;i++) if(iiw[i] == 2) napt2 = napt2+1; //80

	if(napt1 == napt2) goto L100;
		else goto L20;
	L100:;//        continue

	l = 0;
	for(i=1;i<=napt;i++)
	{
		if(iiw[i] == 2)
		{
			l = l+1;
			b[l] = b[i];  // can it cover array data? NO!
			for(j=1;j<=4;j++) xy12[l][j] = xy12[i][j];
		}
	}//300     continue

	for(i=1;i<=l;i++)
	{//805
		icount[i] = 0;
		for(j=1;j<=macoun;j++)
		{
			ijpoin[i][j] = 0;
			xjoin[i][j] = 0.0;
			yjoin[i][j] = 0.0;
		}
	}//805     continue

	napt = l; //total amount of connected fractures
	poin();//+++++++++++++++
	queue();//++++++++++++++

	for(i=nbound+1;i<=l;i++)
	{//302
		ij = icount[i];
		xy12[i][1] = xjoin[i][1];//ijpoin[i][j]
		xy12[i][2] = yjoin[i][1];
		xy12[i][3] = xjoin[i][ij];
		xy12[i][4] = yjoin[i][ij];
	}//302     continue

	fprintf(fp2_output,"%d\n",l);
	for(i=1;i<=l;i++) fprintf(fp2_output,"%10.4lf%10.4lf%10.4lf%10.4lf\n",xy12[i][1],xy12[i][2],xy12[i][3],xy12[i][4]);

	fclose(fp2_output);

	this->cell_parameter_calculate();
      
	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::queue() //acording to length of the connected fracture,regulate ijpoin(i,j)
{
int i,ii,ij,j,jj,ll,c[macoun];
double a[macoun],b[macoun],d[macoun],e[macoun];
double xymax;

xymax=100000.0;
for(i=1;i<=napt;i++){//20

ij = icount[i];
for(j=1;j<=ij;j++) a[j] = sqrt((xjoin[i][j] - xy12[i][1])*(xjoin[i][j] - xy12[i][1]) + (yjoin[i][j] - xy12[i][2])*(yjoin[i][j] - xy12[i][2]));


for(jj=1;jj<=ij;jj++){//18
b[jj] = a[1];
for(ii=1;ii<=ij;ii++){//17
if(a[ii] <= b[jj]){b[jj] = a[ii]; ll=ii;}
}//17 continue

c[jj] = ijpoin[i][ll];
d[jj] = xjoin[i][ll];
e[jj] = yjoin[i][ll];
a[ll] = xymax;
}//18      continue

for(j=1;j<=ij;j++){//19
ijpoin[i][j] = c[j];
xjoin[i][j]  = d[j];
yjoin[i][j]  = e[j];
}//19      continue

}//20      continue

return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
//give the number to the nodes,if the nodes is numbered,or not
int FractureModel::given_no(int ni,int nj)
{
	int i,i1,j1,np1;

	np1 = 0;
	for(i=1;i<=nodes;i++)
	{//10
		i1 = nupoin[i][1];
		j1 = nupoin[i][2];
		if((i1==ni&&j1==nj)||(i1==nj&&j1==ni)) np1 = i;
	}//10  continue

	return np1;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::plane_waterhead()// numbering nodes
{
	int ib2,inode,inode1,inode2,nb2,fenn,qj1,qj2,qj3,qj4;
	int i,j,k,kk,iterma,maj(4),nnode,NRe,m,z,nnumber1,nnumber2,nnumber3,nnumber4;
	double tolmax,b,xl,qq,oldval,err,qq1,qq2,xxl,xxk,sumbl,sumhbl,c2,c;
	double J,kbl,allkb,Vl,reynold,b1,b2,b3,b4,xl1,xl2,xl3,xl4,H1,H2,H3,H4,fenx1,fenx2,fenx0,aa1,bb1,cc1;
	int wshi[5];
	int cshi[5];
	  vec newtonh0(5000);    vec newtonf(5000);     mat newtondif(5000, 5000);
	            
	     
//+++++++++++++++++++++++++confined water flow ++++++++++++++++++++++++++++++



	c=98000.0/12.0 ;   // c is pg/12u and the transferring of the unit
	 c2=464.60533;
	tolmax=0.000001;
       numit=0;

do
{	//----------------------------------------------------------
	numit = numit + 1;// 	

    errmax=0.0;
		
	//---------------------------the sencond boundary condition----------------
 
	ib2 = 0;
	inode  = npb21 - 1;
	inode1 = inode;

	for(nb2=nboun1+1; nb2<=nbound;nb2++)
	{//620
		ib2=ib2+1;
		allkb=0.0;
		for(j=1;j<=icount[nb2];j++)
		{//605
			inode1=inode1+1;
			inode2 = int(nijkl[inode1][maj]);
			b3 = nijkl[inode1][3*maj];
			
			kbl=b3*b3*b3;//kbl=980.0*b3*b3*b3/12.0/0.01;0.01cm**2/s

			allkb = allkb + kbl;
		}//605	continue
		for(k=1;k<=icount[nb2];k++)
		{//610
			inode=inode+1;
			inode2 = int(nijkl[inode][maj]);
			b3 = nijkl[inode][3*maj];
			xl = nijkl[inode][2*maj];
			qq = qb2[ib2]*b3*b3*b3/allkb; 

			if (qq==0)
			  {
				h[inode] = h[inode2];
				newtonf(inode)=h[inode] - h[inode2];
				for(kk=1;kk<=np1n2;kk++)
				{
				 if (kk==inode)      
				 {newtondif(inode,kk)= 1;} 
				 else if (kk==inode2) 
				 {newtondif(inode,kk)=-1;} 
				 else 
				 {newtondif(inode,kk)=0;}
				 
				 newtonh0(inode)=h[inode];
				}
			  }
			else
			  {	Vl=fabs(qq/b3);
				reynold =b3*Vl*2.0/0.01;
				if(reynold <2300)
					{J=qq*12*0.01/b3/b3/b3/980.0;
					h[inode] = h[inode2] + J*xl;}
				else {qq1=fabs(qq)/4.70/b3;
					J=pow(qq1,7.0)*0.01/b3/b3/b3/b3/b3/980.0/980.0/980.0/980.0;
					J=pow(J,0.25)*qq/fabs(qq);
					h[inode] = h[inode2] + J*xl;}
	
	
					newtonf(inode)=h[inode] - h[inode2] - J*xl;
					for(kk=1;kk<=np1n2;kk++)
					{
					if (kk==inode)      
					{newtondif(inode,kk)= 1;} 
					else if (kk==inode2) 
					{newtondif(inode,kk)=-1;} 
					else 
					{newtondif(inode,kk)=0;}
					
					newtonh0(inode)=h[inode];
					}
			  }
				 
				 
//------------------------------------------------the sencond boundary condition---------
				 
		}//610 continue					
	
	}//620	continue


	
	//------------------------------construction of inner nodes matrix begin---------------------------------------------
	
	for(i= np1n1;i<=np1n2;i++)
	{//55 
		        oldval = h[i];
		        sumbl=0.0 ;
		        sumhbl=0.0 ;
				NRe=0;
				m=0;
				z=0;
				nnumber1=0;
				nnumber2=0;
				nnumber3=0;
				nnumber4=0;
				b1=0;
				b2=0;
				b3=0;
				b4=0;
				xl1=0;
				xl2=0;
				xl3=0;
				xl4=0;
				H1=0;
				H2=0;
				H3=0;
				H4=0;
				qj1=0;
				qj2=0;
				qj3=0;
				qj4=0;
				fenx1=0;
				fenx2=0;
				fenx0=0;
				aa1=0;
				bb1=0;
				cc1=0;
				for(j=0;j<5;j++)
					{
						wshi[j]=0;
						cshi[j]=0;
				    }

		        for(j=1;j<=maj;j++)
		        {//15
		        	inode1=i;
		        	inode2 = int(nijkl[inode1][j]);
		        	xl = nijkl[inode1][j+maj];
		        	b = nijkl[inode1][j+2*maj];
		        	if(b == 0.0) continue;
		        	J=(h[inode1]-h[inode2])/xl;
		        	kbl=980*b*b*b/12/0.01;
		        	Vl=kbl*J/b;
		        	Re[inode1][inode2]=2*b*Vl/0.01;
		        
		        	//m = m + 1; cshi[m] = j;
		        	if(fabs(Re[inode1][inode2])>2300)
		        	{NRe=NRe+1;z=z+1;wshi[z]=j;}
		        	else {m=m+1;cshi[m]=j;}
		        
		        }//15 continue
                    if (NRe<1)
			         {
						 //h[i]=sumhbl/sumbl ;
						 fenn=0;
						qj1=int (cshi[1]);
						if(qj1>0)
						{nnumber1=int(nijkl[i][qj1]);//np0 is connected to np1，np2，np3，np4
						 b1=nijkl[i][qj1+2*maj];
						 xl1=nijkl[i][qj1+maj];
						 H1=h[nnumber1];}
						 else{ H1=h[i];b1=0;xl1=1;}

						qj2=int (cshi[2]);
						if(qj2>0)
						{nnumber2=int(nijkl[i][qj2]);
						 b2=nijkl[i][qj2+2*maj];
						 xl2=nijkl[i][qj2+maj];
						 H2=h[nnumber2];}
						 else{ H2=h[i];b2=0;xl2=1;}

						 qj3=int (cshi[3]);
						 if(qj3>0)
						 {nnumber3=int(nijkl[i][qj3]);
						 b3=nijkl[i][qj3+2*maj];
						 xl3=nijkl[i][qj3+maj];
						 H3=h[nnumber3];}
						 else{ H3=h[i];b3=0;xl3=1;}

						 qj4=int (cshi[4]);
						 if(qj4>0)
						 {nnumber4=int(nijkl[i][qj4]);
						 b4=nijkl[i][qj4+2*maj];
						 xl4=nijkl[i][qj4+maj];
						 H4=h[nnumber4];}
						else{ H4=h[i];b4=0;xl4=1;}
						
						newtonf(i)=c*b1*b1*b1*(h[i]-H1)/xl1+c*b2*b2*b2*(h[i]-H2)/xl2+c*b3*b3*b3*(h[i]-H3)/xl3+c*b4*b4*b4*(h[i]-H4)/xl4;
						for(kk=1;kk<=np1n2;kk++)
						{
						 if (kk==nnumber1)      
						 {newtondif(i,kk)=-c*b1*b1*b1/xl1;} 
						 else if (kk==nnumber2) 
						 {newtondif(i,kk)=-c*b2*b2*b2/xl2;} 
						 else if (kk==nnumber3) 
						 {newtondif(i,kk)=-c*b3*b3*b3/xl3;}
						 else if (kk==nnumber4) 
						 {newtondif(i,kk)=-c*b4*b4*b4/xl4;}
						 else if (kk==i) 
						 {newtondif(i,kk)=c*b1*b1*b1/xl1+c*b2*b2*b2/xl2+c*b3*b3*b3/xl3+c*b4*b4*b4/xl4;}
						 else 
						 {newtondif(i,kk)=0;}
						 
						 newtonh0(i)=h[i];
						}
						
					}

					if (NRe==1)
					{
						 fenn=0;
						 qj1=int (wshi[1]);
						if(qj1>0)
						{nnumber1=int(nijkl[i][qj1]);
						 b1=nijkl[i][qj1+2*maj];
						 xl1=nijkl[i][qj1+maj];
                         H1=h[nnumber1];}
						 else{ H1=h[i]+0.00001;b1=0;xl1=1;}

						 qj2=int (cshi[1]);
						if(qj2>0)
						{nnumber2=int(nijkl[i][qj2]);
						 b2=nijkl[i][qj2+2*maj];
						 xl2=nijkl[i][qj2+maj];
                         H2=h[nnumber2];}
						 else{ H2=h[i];b2=0;xl2=1;}

						 qj3=int (cshi[2]);
						if(qj3>0)
						{nnumber3=int(nijkl[i][qj3]);
						 b3=nijkl[i][qj3+2*maj];
						 xl3=nijkl[i][qj3+maj];
                          H3=h[nnumber3];}
						 else{ H3=h[i];b3=0;xl3=1;}

						  qj4=int (cshi[3]);
						if(qj4>0)
						{nnumber4=int(nijkl[i][qj4]);
						 b4=nijkl[i][qj4+2*maj];
						 xl4=nijkl[i][qj4+maj];
                          H4=h[nnumber4];}
						else{ H4=h[i];b4=0;xl4=1;}
						
						newtonf(i)=c2*pow(b1,1.7143)*(h[i]-H1)/fabs(h[i]-H1)*pow(fabs((h[i]-H1)/xl1),0.5714)+c*b2*b2*b2*(h[i]-H2)/xl2+c*b3*b3*b3*(h[i]-H3)/xl3+c*b4*b4*b4*(h[i]-H4)/xl4;
						for(kk=1;kk<=np1n2;kk++)
						{
						 if (kk==nnumber1)      
						 { newtondif(i,kk)=-c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1;}
						 else if (kk==nnumber2) 
						 {newtondif(i,kk)=-c*b2*b2*b2/xl2;} 
						 else if (kk==nnumber3) 
						 {newtondif(i,kk)=-c*b3*b3*b3/xl3;}
						 else if (kk==nnumber4) 
						 {newtondif(i,kk)=-c*b4*b4*b4/xl4;}
						 else if (kk==i) 
						 {newtondif(i,kk)=c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1+c*b2*b2*b2/xl2+c*b3*b3*b3/xl3+c*b4*b4*b4/xl4;}
						 else
						{newtondif(i,kk)=0;}
						 
						 newtonh0(i)=h[i];
						}
					}

					if (NRe==2)
					 {
						 fenn=0;
						 qj1=int (wshi[1]);
						if(qj1>0)
						{nnumber1=int(nijkl[i][qj1]);
						 b1=nijkl[i][qj1+2*maj];
						 xl1=nijkl[i][qj1+maj];
                         H1=h[nnumber1];}
						 else{ H1=h[i]+0.00001;b1=0;xl1=1;}

						 qj2=int (wshi[2]);
						if(qj2>0)
						{nnumber2=int(nijkl[i][qj2]);
						 b2=nijkl[i][qj2+2*maj];
						 xl2=nijkl[i][qj2+maj];
                         H2=h[nnumber2];}
						 else{ H2=h[i]+0.00001;b2=0;xl2=1;}

						 qj3=int (cshi[1]);
						if(qj3>0)
						{nnumber3=int(nijkl[i][qj3]);
						 b3=nijkl[i][qj3+2*maj];
						 xl3=nijkl[i][qj3+maj];
                         H3=h[nnumber3];}
						 else{ H3=h[i];b3=0;xl3=1;}

						 qj4=int (cshi[2]);
						if(qj4>0)
						{nnumber4=int(nijkl[i][qj4]);
						 b4=nijkl[i][qj4+2*maj];
						 xl4=nijkl[i][qj4+maj];
                         H4=h[nnumber4];}
						else{ H4=h[i];b4=0;xl4=1;}
						
						newtonf(i)=c2*pow(b1,1.7143)*(h[i]-H1)/fabs(h[i]-H1)*pow(fabs((h[i]-H1)/xl1),0.5714)+c2*pow(b2,1.7143)*(h[i]-H2)/fabs(h[i]-H2)*pow(fabs((h[i]-H2)/xl2),0.5714)+c*b3*b3*b3*(h[i]-H3)/xl3+c*b4*b4*b4*(h[i]-H4)/xl4;
						for(kk=1;kk<=np1n2;kk++)
						{
						 if (kk==nnumber1)      
						 {newtondif(i,kk)=-c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1;}
						 else if (kk==nnumber2) 
						 {newtondif(i,kk)=-c2*pow(b2,1.7143)*pow(fabs((h[i]-H2)/xl2),-0.4286)*0.5714/xl2;} 
						 else if (kk==nnumber3) 
						 {newtondif(i,kk)=-c*b3*b3*b3/xl3;}
						 else if (kk==nnumber4) 
						 {newtondif(i,kk)=-c*b4*b4*b4/xl4;}
						 else if (kk==i) 
						 {newtondif(i,kk)=c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1+c2*pow(b2,1.7143)*pow(fabs((h[i]-H2)/xl2),-0.4286)*0.5714/xl2+c*b3*b3*b3/xl3+c*b4*b4*b4/xl4;}
						 else
						 {newtondif(i,kk)=0;}
						 
						 newtonh0(i)=h[i];
						
						}
						
					}

					if (NRe==3)
					 {
						 fenn=0;
						 qj1=int (wshi[1]);
						if(qj1>0)
						{nnumber1=int(nijkl[i][qj1]);
						 b1=nijkl[i][qj1+2*maj];
						 xl1=nijkl[i][qj1+maj];
                         H1=h[nnumber1];}
						 else{ H1=h[i]+0.00001;b1=0;xl1=1;}

						 qj2=int (wshi[2]);
						if(qj2>0)
						{nnumber2=int(nijkl[i][qj2]);
						 b2=nijkl[i][qj2+2*maj];
						 xl2=nijkl[i][qj2+maj];
                         H2=h[nnumber2];}
						 else{ H2=h[i]+0.00001;b2=0;xl2=1;}

						 qj3=int (wshi[3]);
						if(qj3>0)
						{nnumber3=int(nijkl[i][qj3]);
						 b3=nijkl[i][qj3+2*maj];
						 xl3=nijkl[i][qj3+maj];
                          H3=h[nnumber3];}
						 else{ H3=h[i]+0.00001;b3=0;xl3=1;}

						 qj4=int (cshi[1]);
						if(qj4>0)
						{nnumber4=int(nijkl[i][qj4]);
						 b4=nijkl[i][qj4+2*maj];
						 xl4=nijkl[i][qj4+maj];
                         H4=h[nnumber4];}
						else{ H4=h[i];b4=0;xl4=1;}
						
						newtonf(i)=c2*pow(b1,1.7143)*(h[i]-H1)/fabs(h[i]-H1)*pow(fabs((h[i]-H1)/xl1),0.5714)+c2*pow(b2,1.7143)*(h[i]-H2)/fabs(h[i]-H2)*pow(fabs((h[i]-H2)/xl2),0.5714)+c2*pow(b3,1.7143)*(h[i]-H3)/fabs(h[i]-H3)*pow(fabs((h[i]-H3)/xl3),0.5714)+c*b4*b4*b4*(h[i]-H4)/xl4;
						for(kk=1;kk<=np1n2;kk++)
						{
						 if (kk==nnumber1)      
						 {newtondif(i,kk)=-c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1;}
						 else if (kk==nnumber2) 
						 {newtondif(i,kk)=-c2*pow(b2,1.7143)*pow(fabs((h[i]-H2)/xl2),-0.4286)*0.5714/xl2;} 
						 else if (kk==nnumber3) 
						 {newtondif(i,kk)=-c2*pow(b3,1.7143)*pow(fabs((h[i]-H3)/xl3),-0.4286)*0.5714/xl3;} 
						 else if (kk==nnumber4) 
						 {newtondif(i,kk)=-c*b4*b4*b4/xl4;}
						 else if (kk==i) 
						 {newtondif(i,kk)=c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1+c2*pow(b2,1.7143)*pow(fabs((h[i]-H2)/xl2),-0.4286)*0.5714/xl2+c2*pow(b3,1.7143)*pow(fabs((h[i]-H3)/xl3),-0.4286)*0.5714/xl3+c*b4*b4*b4/xl4;}
						 else 
						 {newtondif(i,kk)=0;}
						 
						 newtonh0(i)=h[i];
						
						}
						
					}

					if (NRe==4)
					 {
						 fenn=0;
						 qj1=int (wshi[1]);
						if(qj1>0)
						{nnumber1=int(nijkl[i][qj1]);
						 b1=nijkl[i][qj1+2*maj];
						 xl1=nijkl[i][qj1+maj];
                         H1=h[nnumber1];}
						 else{ H1=h[i]+0.00001;b1=0;xl1=1;}

						 qj2=int (wshi[2]);
						if(qj2>0)
						{nnumber2=int(nijkl[i][qj2]);
						 b2=nijkl[i][qj2+2*maj];
						 xl2=nijkl[i][qj2+maj];
                         H2=h[nnumber2];}
						 else{ H2=h[i]+0.00001;b2=0;xl2=1;}

						 qj3=int (wshi[3]);
						if(qj3>0)
						{nnumber3=int(nijkl[i][qj3]);
						 b3=nijkl[i][qj3+2*maj];
						 xl3=nijkl[i][qj3+maj];
                          H3=h[nnumber3];}
						 else{ H3=h[i]+0.00001;b3=0;xl3=1;}

						 qj4=int (wshi[4]);
						if(qj4>0)
						{nnumber4=int(nijkl[i][qj4]);
						 b4=nijkl[i][qj4+2*maj];
						 xl4=nijkl[i][qj4+maj];
                         H4=h[nnumber4];}
						else{ H4=h[i]+0.00001;b4=0;xl4=1;}
						
						newtonf(i)=c2*pow(b1,1.7143)*(h[i]-H1)/fabs(h[i]-H1)*pow(fabs((h[i]-H1)/xl1),0.5714)+c2*pow(b2,1.7143)*(h[i]-H2)/fabs(h[i]-H2)*pow(fabs((h[i]-H2)/xl2),0.5714)+c2*pow(b3,1.7143)*(h[i]-H3)/fabs(h[i]-H3)*pow(fabs((h[i]-H3)/xl3),0.5714)+c2*pow(b4,1.7143)*(h[i]-H4)/fabs(h[i]-H4)*pow(fabs((h[i]-H4)/xl4),0.5714);
						for(kk=1;kk<=np1n2;kk++)
						{
						 if (kk==nnumber1)      
						 {newtondif(i,kk)=-c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1;}
						 else if (kk==nnumber2) 
						 {newtondif(i,kk)=-c2*pow(b2,1.7143)*pow(fabs((h[i]-H2)/xl2),-0.4286)*0.5714/xl2;} 
						 else if (kk==nnumber3) 
						 {newtondif(i,kk)=-c2*pow(b3,1.7143)*pow(fabs((h[i]-H3)/xl3),-0.4286)*0.5714/xl3;} 
						 else if (kk==nnumber4) 
						 {newtondif(i,kk)=-c2*pow(b4,1.7143)*pow(fabs((h[i]-H4)/xl4),-0.4286)*0.5714/xl4;}
						 else if (kk==i) 
						 {newtondif(i,kk)=c2*pow(b1,1.7143)*pow(fabs((h[i]-H1)/xl1),-0.4286)*0.5714/xl1+c2*pow(b2,1.7143)*pow(fabs((h[i]-H2)/xl2),-0.4286)*0.5714/xl2+c2*pow(b3,1.7143)*pow(fabs((h[i]-H3)/xl3),-0.4286)*0.5714/xl3+c2*pow(b4,1.7143)*pow(fabs((h[i]-H4)/xl4),-0.4286)*0.5714/xl4;}
						 else 
						 {newtondif(i,kk)=0;}
						 
						 newtonh0(i)=h[i];
						
						}
					
					}

	}//55 
	
	//------------------------------construction of inner nodes matrix ---------------------------------------------
	

	         
	           int iii;     iii = np1n2 - npb21 + 1;
	         mat newtondifA(iii, iii);
	         
	         vec newtonh0A(iii);    vec newtonfA(iii);
	         vec newtonh0B(iii);    vec newtontol(iii);

	
	        for(i=0;i<=(np1n2-npb21);i++)
			{
				      int jjj;  jjj =i+npb21 ;
				newtonfA(i) =newtonf(jjj);
				newtonh0A(i) = newtonh0(jjj);
				for(j=0;j<=(np1n2-npb21);j++)
				{
					int kkk;      kkk=j+npb21;
					newtondifA(i,j) =newtondif(jjj,kkk);
				}
			}

			newtontol = solve(newtondifA,newtonfA);
				newtonh0B=newtonh0A-newtontol;
	
				for(i=0;i<=(np1n2-npb21);i++)
				{if (fabs(newtontol(i))>errmax) 	errmax=fabs(newtontol(i));}
			
				for(i=0;i<=(np1n2-npb21);i++)
				{ h[i+ npb21]=newtonh0B(i);}
	

} while (errmax > tolmax && numit<200); 
  
	
//+++++++++++++++++++++newton iteration completed +++++++++++++++++++++++++++++++



//--------------------------water budget and error analysis-----------------------------


	for(inode=1;inode<=nodes;inode++)
		   for(j=1;j<=maj;j++)
			{
				inode1=int(nijkl[inode][j]);
				if(inode1<1)continue;
				b = nijkl[inode][2*maj+j];
				xl = nijkl[inode][maj+j];
				J=(h[inode]-h[inode1])/xl;
				kbl=980.0*b*b*b/12.0/0.01;
				Vl=kbl*J/b;
				Re[inode][inode1]=2*b*Vl/0.01;
				if (fabs(Re[inode][inode1])<2300){Qcell[inode][j]=kbl*J;}
			    else{Qcell[inode][j]=464.60533*(h[inode]-h[inode1])/fabs(h[inode]-h[inode1])*pow(b,1.7143)*pow(fabs(J),0.5714);
					Re[inode][inode1]=2*Qcell[inode][j]/0.01;}//
			}


	qq1 =0.0;
	qq2=0.0;
	qoutmax=0.0;
	for(inode=1;inode<=npb22;inode++)
	{//2000
	  q[inode] = Qcell[inode][maj];
	  if(qoutmax>q[inode])qoutmax=q[inode];
	  qq1 = qq1 + q[inode];
	  qq2 = qq2 + fabs(q[inode]);
	}//2000 continue
	qerror = qq1;//qq2*200.0;
	qoutmax=qoutmax/qq2*200.0;
	   qq2 = 0.5*qq2;
	xxl=sqrt((xy12[1][3] - xy12[1][1])*(xy12[1][3] - xy12[1][1]) + (xy12[1][2] - xy12[1][4])*(xy12[1][2] - xy12[1][4]));
	xxk=qq2/xxl*0.01; 

	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::cell_parameter_calculate()
{
  int i,i1,j,n,n1,nb1,j1,ij,maj,np,np1; 
  double detax,detay;

	FILE *fp20_output;
    if((fp20_output=fopen("Out_ConnectedFracNetwork.dat","a"))==NULL)
	{
		std::cout << "fail to open Out_ConnectedFracNetwork.dat\n";
      return;
	}

    maj = 4; 
    for(i=1;i<=npmax;i++)
	{//508
       for(j=1;j<=2;j++) nupoin[i][j] = 0;
	}//508  

    np = 0;
    for(i=1;i<=napt;i++)
	{//500
       ij = icount[i];
       for(j=1;j<=ij;j++)
	   {//450
          j1 = ijpoin[i][j];
           if(j1 >= i)
		   {
              np = np+1;
              xynp[np][1] = xjoin[i][j];
              xynp[np][2] = yjoin[i][j];
              nupoin[np][1] = i;
              nupoin[np][2] = ijpoin[i][j];
		   }
	   }//450     
	}//500  
	
    nodes = np; //total amount of nodes

    fprintf(fp20_output,"%d\n",nodes);
	for(i=1;i<=nodes;i++)
	{
		fprintf(fp20_output,"%10.4lf %10.4lf\n",xynp[i][1],xynp[i][2]);
	}

	fclose(fp20_output);

//       maj is the total nodes conected with one node
//       nodes is the total nodes
//       npb11-npb12,np1n1-np1n2,npb21-npb22 are the numbers
//       begining-ending number given to the nodes of water table boundery,inerior, water flow boundery, respectively.
//       nijkl(conected number+dx+b)

	npb11 = 1;
	npb12 = 0;
	n = 0;
	for(i=1;i<=nboun1;i++)npb12 = npb12 + icount[i];
	npb21 = npb12 + 1;
	npb22 = npb12;

	nb1 = nboun1 + 1;
	for(i=nb1;i<=nboun1+nboun2;i++) npb22 = npb22 + icount[i];

	npbf1=npb22+1;
	npbf2=npb22;
    for(i=nboun1+nboun2+1;i<=nbound;i++) npbf2 = npbf2 + icount[i];

	np1n1 = npbf2 + 1;
	np1n2 = nodes;
	for(np=1;np<=nodes;np++)
	{//600
		if(np <= npbf2)
		{// then
			nijkl[np][1] = 0.0;
			nijkl[np][1+maj] = 0.0;
			nijkl[np][1+maj*2] = 0.0;
			nijkl[np][2] = 0.0;
			nijkl[np][2+maj] = 0.0;
			nijkl[np][2+maj*2] = 0.0;
			nijkl[np][3] = 0.0;
			nijkl[np][3+maj] = 0.0;
			nijkl[np][3+maj*2] = 0.0;
			i = nupoin[np][1];
			j = nupoin[np][2];
			ij = icount[j];
			if(ijpoin[j][1]  == i)  i1 = ijpoin[j][2];
			if(ijpoin[j][ij] == i)  i1 = ijpoin[j][ij-1];
			np1 = given_no (j,i1);//-----------//what does it mean when np1=0;
			nijkl[np][4] = np1;
			detax = xynp[np][1]-xynp[np1][1];
			detay = xynp[np][2]-xynp[np1][2];
			nijkl[np][4+maj]=sqrt(detax*detax + detay*detay);
			nijkl[np][4+maj*2] = b[j];
			goto L600;
		}//end if

		i  = nupoin[np][1];
		j  = nupoin[np][2];
		ij = icount[i];
		for(n=1;n<=ij;n++)
		{//501
			if(ijpoin[i][n] == j) j1=n;
		}// 501 continue

		if(ijpoin[i][j1-1] == 0)
		{// then
			nijkl[np][1]       = 0.0;
			nijkl[np][1+maj]   = 0.0;
			nijkl[np][1+maj*2] = 0.0;
		}//end if

		if(ijpoin[i][j1-1] != 0)
		{// then
			np1 = given_no(i,ijpoin[i][j1-1]);
			nijkl[np][1] = np1;
			detax = xynp[np][1]-xynp[np1][1];
			detay = xynp[np][2]-xynp[np1][2];
			nijkl[np][1+maj] = sqrt(detax*detax + detay*detay);
			nijkl[np][1+maj*2] = b[i];
		}//end if

		if(ijpoin[i][j1+1] == 0)
		{//then
			nijkl[np][2]       =0.0;
			nijkl[np][2+maj]   =0.0;
			nijkl[np][2+maj*2] =0.0;
		}//end if

		if(ijpoin[i][j1+1] != 0)
		{//then
			np1 = given_no(i,ijpoin[i][j1+1]);
			nijkl[np][2] = np1;
			detax = xynp[np][1]-xynp[np1][1];
			detay = xynp[np][2]-xynp[np1][2];
			nijkl[np][2+maj] = sqrt(detax*detax + detay*detay);
			nijkl[np][2+maj*2]=b[i];
		}//end if

		ij = icount[j];
		for(n=1;n<=ij;n++)
		{//502
			if(ijpoin[j][n] == i) i1 = n;
		}//502     continue

		if(ijpoin[j][i1-1] == 0)
		{//then
			nijkl[np][3]       = 0.0;
			nijkl[np][3+maj]   = 0.0;
			nijkl[np][3+maj*2] = 0.0;
		}//end if

		if(ijpoin[j][i1-1] != 0)
		{// then
			np1 = given_no(j,ijpoin[j][i1-1]);
			nijkl[np][3] = np1;
			detax = xynp[np][1]-xynp[np1][1];
			detay = xynp[np][2]-xynp[np1][2];
			nijkl[np][3+maj] = sqrt(detax*detax + detay*detay);
			nijkl[np][3+maj*2] = b[j];
		}//end if

		if(ijpoin[j][i1+1] == 0)
		{// then
			nijkl[np][4]       = 0.0;
			nijkl[np][4+maj]   = 0.0;
			nijkl[np][4+maj*2] = 0.0;
		}//end if

		if(ijpoin[j][i1+1] != 0)
		{// then
			np1 = given_no(j,ijpoin[j][i1+1]);
			nijkl[np][4] = np1;
			detax = xynp[np][1]-xynp[np1][1];
			detay = xynp[np][2]-xynp[np1][2];
			nijkl[np][4+maj] = sqrt(detax*detax + detay*detay);
			nijkl[np][4+maj*2] = b[j];
		}//end if

	L600:;
	}//600  continue

	//-------------------------- Initiation for Water head calculation --------------------------

	n = 0;
	Remax=0;
	h[0]=0.0; 
	for(i=1;i<=nboun1;i++)
	{//81
		n1 = icount[i];
		for(j=1;j<=n1;j++)
		{//73
			n = n + 1;
			h[n] = hb1[i];
		}//73 
	}//81 

	for(i=1;i<=nodes;i++)
	  for(j=1;j<=4;j++)	
	  {
		// Recell[i][j]=0.0;
		  Qcell[i][j]=0.0;
	  }


	maxx=0.0;
	maxy=0.0;
	for(i=1;i<=nodes;i++)
	{
		height[i]=xynp[i][2];
		if(xynp[i][1]>maxx) maxx=xynp[i][1] ; 
		if(xynp[i][2]>maxy) maxy=xynp[i][2] ; 
	}
	for(i=npb21;i<=nodes;i++)h[i]=0.9999*maxy;
	//-----------------------------------   !   ----------------------------------------------
	
	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::outputdata()
{
	int i,j,k;
	
	FILE *fp3_output;
	if((fp3_output=fopen("Out_FractureReadToKarst.dat","w"))==NULL)
	{
		std::cout << "fail to open Out_FractureReadToKarst.dat.dat\n";
		return;
	}

//"fraturereadtokarst.dat"****************************
	fprintf(fp3_output, "初始水头计算迭代次数= %3d\n", numit);//t=0
	fprintf(fp3_output, "\n 最大水头误差errmax(m) %10.8lf\n", errmax);
	fprintf(fp3_output, "\n 水均衡qerror(cm2/s)= %10.8lf\n", qerror);
//	fprintf(fp3_output, "\n qerror(%%)= %10.5lf\n", qerror);
	fprintf(fp3_output, "\n 最大流量占总流量的百分比Max-discharge(%%)%10.5lf\n\n", qoutmax);
											//water budget
	fprintf(fp3_output, "Water Head:\n");
	for (i = 1; i <= nodes; i++) fprintf(fp3_output, "%2d  %10.4lf\n", i, h[i]);


	fprintf(fp3_output,"%5d %5d %5d\n",nboun1,nboun2,nodes);
    
	for(i=1;i<=nbound;i++) fprintf(fp3_output,"%5d ",icount[i]);
    fprintf(fp3_output,"\n");
	
	fprintf(fp3_output,"%5d %5d %5d %5d %5d %5d\n",npb11,npb12,np1n1,np1n2,npb21,npb22);
	
	for(i=1;i<=nboun1;i++)fprintf(fp3_output,"%10.4lf ",hb1[i]);
    fprintf(fp3_output,"\n");
	
	for(i=1;i<=nboun2;i++)fprintf(fp3_output,"%10.4lf ",qb2[i]);
    fprintf(fp3_output,"\n");
    
	for(i=1;i<=nodes;i++) fprintf(fp3_output,"%10.4lf %10.4lf\n",xynp[i][1],xynp[i][2]);



	for(i=1;i<=nodes;i++)
	{
		for(j=1;j<=12;j++)fprintf(fp3_output,"%10.4lf ",nijkl[i][j]);
		fprintf(fp3_output,"\n");
	}


	for(i=1;i<=nodes;i++)
	{
		for(j=1;j<=4;j++)
		{
			    if(h[i]==0)continue;//if(h[i]<0.001)continue;
			    k=int(nijkl[i][j]);
				if(k<1)continue;
				if(i<np1n1 && k<np1n1)continue;//
			fprintf(fp3_output,"雷诺数 %3d %3d %15.4f",i,k,Re[i][k]);
		    fprintf(fp3_output,"  流量 %3d %3d %15.4lf\n",i,k,Qcell[i][j]);
		    //fprintf(fp3_output,"\n");
		}
	}
	fclose(fp3_output);
	return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void FractureModel::profile_waterhead()
{
	int i,ii,iii,j,jj,jjj,k,kk,iterma,maj(4),inode, inode0, inode1,inode2, inode3,inode4,inode5,ib2,nb2,maxerrnode1,maxerrnode2,m,z,gan,gan2,fenn;
	int nnumber,nnode,available1, available2, ninfil,NRe,nnumber1,nnumber2,nnumber3,nnumber4,qj1,qj2,qj3,qj4;
	double tolmax,errmax1,errmax2,sumbl,sumhbl,sumqbl,err1,err2,reynold;
	double Qspring,b,b1,b2,b3,b4,xl,xl1,xl2,xl3,xl4,H1,H2,H3,H4,q3,q4,c,c2,oldval,alpha,J,kbl,allkb,Vl,kbt,qq,qq1,qq2;
	double tran_coeff,Qsurfcorro1,sumwidth,fenx1,fenx2,fenx0,aa1,bb1,cc1;
	int wshi[5];
	int cshi[5];
	double qbalance[5000], qqbudget; 
	
	  vec newtonh0(5000);    vec newtonf(5000);     mat newtondif(5000, 5000);      mat newtondifA0(np1n2, np1n2);
	      newtonh0.fill(0);      newtonf.fill(0);       newtondif.fill(0);	            newtondifA0.fill(0);
	
	iterma=100000;
	tolmax=0.0001;

	errmax=0;

	//-------------------------- Initiation for original Water head calculation -----------------------------
	//----------------------------------------------   !   --------------------------------------------------
	//----------------------------------------------   !   --------------------------------------------------


	int n,n1;n = 0;n1=0;
	h[0]=0.0;

	for(i=1;i<=nboun1;i++)
	{//81
		n1 = icount[i];
		for(j=1;j<=n1;j++)
		{//73
			n = n + 1;
			h[n] = hb1[i];
		}//73 
	}//81 

	for(i=1;i<=nodes;i++)
	  for(j=1;j<=4;j++)	
	  {		
		  Qcell[i][j]=0.0;
	  }
	//--------------------------Assuming a lower confined water table-----------------------------------------------

	maxx=0.0;
	maxy=0.0;
	for(i=1;i<=nodes;i++)
	{
		height[i]=xynp[i][2];
		if(xynp[i][1]>maxx) maxx=xynp[i][1] ; 
		if(xynp[i][2]>maxy) maxy=xynp[i][2] ; 
	}
	for (i = np1n1; i <= nodes; i++)
	{
		h[i] = 0;
		if (height[i] <= 100.0)
			h[i] = 200.0 + 100 * ( maxx - xynp[i][1]) / maxx;
	}



	//-----------------------------------   !   ----------------------------------------------
	//-------------------------- Inner iteration for the original cofined water flow --------------------------

	numit = 0;
	do   
	{                                 


		errmax1 = 0.0;
		errmax2 = 0.0;
		errmax = 0.0;
		numit = numit + 1;

		int n, n1; n = 0; n1 = 0;
		h[0] = 0.0;

		for (i = 1; i <= nboun1; i++)
		{//81
			n1 = icount[i];
			for (j = 1; j <= n1; j++)
			{//73
				n = n + 1;
				h[n] = hb1[i];
				newtonf(n) = h[n]- hb1[i];
				newtondif(n, n) = 1;
				newtonh0(n) = h[n];
			}//73 
		}//81 



		for (i = np1n1; i <= nodes; i++) 
		{
			if (h[i] < height[i])
				h[i] = 0;
		}



		//-------------------------------------the second boundary condition----------------------------------
		//----------------------------------------------------------------------------------------------------
		ib2 = 0;
		inode = npb21 - 1;
		inode1 = inode;

		for (nb2 = nboun1 + 1; nb2 <= nboun1 + nboun2; nb2++)//
		{//620------------------------
			ib2 = ib2 + 1;
			allkb = 0.0;
			for (j = 1; j <= icount[nb2]; j++)//
			{//605
				inode1 = inode1 + 1;
				inode2 = int(nijkl[inode1][maj]);
				b3 = nijkl[inode1][3 * maj];
				
				kbl = b3 * b3 * b3;//kbl=980.0*b3*b3*b3/12.0/0.01;0.01cm**2/s
				
				allkb = allkb + kbl;
			}//605	continue

			for (k = 1; k <= icount[nb2]; k++)
			{//610
				inode = inode + 1;
				inode2 = int(nijkl[inode][maj]);
				b3 = nijkl[inode][3 * maj];
				xl = nijkl[inode][2 * maj];
				qq = qb2[ib2] * b3 * b3 * b3 / allkb; 

				if (qq == 0)
				{
					h[inode] = h[inode2];
					newtonf(inode) = h[inode] - h[inode2];
					for (kk = 1; kk <= np1n2; kk++)
					{
						if (kk == inode)
						{
							newtondif(inode, kk) = 1;
						}
						else if (kk == inode2)
						{
							newtondif(inode, kk) = -1;
						}
						else
						{
							newtondif(inode, kk) = 0;
						}

						newtonh0(inode) = h[inode];
					}
				}
				else
				{
					Vl = fabs(qq / b3);
					reynold = b3 * Vl * 2.0 / 0.01;
					if (reynold < 2300)
					{
						J = qq * 12 * 0.01 / b3 / b3 / b3 / 980.0;
						h[inode] = h[inode2] + J * xl;
					}
					else {
						qq1 = fabs(qq) / 4.70 / b3;
						J = pow(qq1, 7.0) * 0.01 / b3 / b3 / b3 / b3 / b3 / 980.0 / 980.0 / 980.0 / 980.0;
						J = pow(J, 0.25) * qq / fabs(qq);//           
					}


					newtonf(inode) = h[inode] - h[inode2] - J * xl;
					for (kk = 1; kk <= np1n2; kk++)
					{
						if (kk == inode)
						{
							newtondif(inode, kk) = 1;
						}
						else if (kk == inode2)
						{
							newtondif(inode, kk) = -1;
						}
						else
						{
							newtondif(inode, kk) = 0;
						}

						newtonh0(inode) = h[inode];
					}
				}
			}//610 continue
		}//620	continue
	//------------------------------------------------the second boundary condition------------------------------------------



	//-------------------------------------------------the third boundary condition-------------------------------------------------

		for (i = npbf1; i <= npbf2; i++)
		{
			h[i] = 0.0;//
			newtonf(i) = h[i];
			newtondif(i, i) = 1;
			newtonh0(i) = h[i];
		}


	//the following deals with the spring boundary----------------------------------------------------------------------------

		Qspring = 0.0;
		for (i = np1n1; i <= nodes; i++)//
			for (j = 1; j <= 4; j++)
			{
				nnumber = int(nijkl[i][j]);//
				b = nijkl[i][j + 2 * maj];
				xl = nijkl[i][j + maj];
				if (fabs(xynp[nnumber][1]) > 0.01 && xynp[nnumber][1] < (maxx - 0.01))continue;
				if (nnumber > npb22 && nnumber < np1n1)
				{
					J = (h[i] - height[nnumber]) / xl;
					kbl = 980.0 * b * b * b / 12.0 / 0.01;
					Vl = kbl * J / b;
					Re[i][nnumber] = 2 * b * Vl / 0.01;
					if (Re[i][nnumber] > 2300) { Qspring = 464.60533 * pow(b, 1.7143) * pow(fabs(J), 0.5714); }
					else Qspring = kbl * J;
					// when Qsrping is greater than 0 , it means it can be regared as a spring
					if (h[nnumber] < (height[nnumber]))//20220322
					{
						h[nnumber] = 0;
						newtonf(nnumber) = h[nnumber];
						newtondif(nnumber, nnumber) = 1;
						newtonh0(nnumber) = h[nnumber];
					}

					if (h[nnumber] >= (height[nnumber]))//20220322
					{
						h[nnumber] = height[nnumber];//目的找到内点相邻的渗出边界节点
						newtonf(nnumber) = h[nnumber] - height[nnumber];
						newtondif(nnumber, nnumber) = 1;
						newtonh0(nnumber) = h[nnumber];
					}
				}
			}
		//---------------------------------------------------------------------------------------------------------------------



		//---------------------------------------------------------------------------------------------------------------
		//------------------------------------------rainfall recharge on the water table --------------------------------
		//---------------------------------------------------------------------------------------------------------------


		//----------------------------------find the number that total flow rate is divided-------------------------------
		for (i = 0; i < 500; i++)
			for (j = 0; j < 4; j++)  qab2[i][j] = 0;


		number_flow = 0;
		for (i = np1n1; i <= nodes; i++)
		{			
			nnumber1 = 0; nnumber2 = 0;
			if (height[i] <= h[i]) //wet node
			{
				for (j = 1; j <= 4; j++)
				{
					nnumber = int(nijkl[i][j]);
					if (nnumber == 0)continue;
					if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

					if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
					{
						nnumber1 = nnumber1 + 1;
					}

					if (h[nnumber] >= height[nnumber])
					{
						nnumber2 = nnumber2 + 1;
					}
				}

				if (nnumber1 != 0 && nnumber2 >= 2)
				   for (j = 1; j <= 4; j++)
				   {
					  nnumber = int(nijkl[i][j]);
					  if (nnumber == 0)continue;
					  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

					   if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
                       {  
					   number_flow = number_flow + 1;
					   qab2[number_flow][1] = i;					
					   qab2[number_flow][2] = nnumber;
					   }

				    }

				
			}						   
		}



		//-----------------the following calculates the flow rate of each fracture segment-------------
		tran_coeff = 10.0 / 365.0 / 24.0 / 60.0 / 60.0;  //unit width 1 cm
				 //the unit of Qsurfcorro is mm/a, and the unit of Qsurfcorro1 is ml/s
				 //tran_coeff is the transferring coefficient.
		Qsurfcorro1 = tran_coeff * (Qsurfcorro * (maxx - 0.0) * 1.0);//*(maxx-minx)
		sumwidth = 0;
		ninfil = 0;
		for (i = 1; i <= number_flow; i++)
		{
			inode = int(qab2[i][1]);
			for (k = 1; k <= 4; k++)
			{
				nnumber = int(nijkl[inode][k]);
				
				if (nnumber == int(qab2[i][2]))
				{
					sumwidth = sumwidth + nijkl[inode][k + 2 * maj];
					ninfil++;
				}
			}
		}
		for (i = 1; i <= number_flow; i++)
		{
			inode = int(qab2[i][1]);
			for (k = 1; k <= 4; k++)
			{
				nnumber = int(nijkl[inode][k]);
				if (nnumber == int(qab2[i][2]))//&& h[inode]<maxy
				{
					qab2[i][3] = Qsurfcorro1 / ninfil;//evenly distribute infiltration																	
				}

			}

		}
 //------------------------------------------------boundary conditions completed------------------------------------------
//------------------------------------------------------------------------------------------------------------------------

 


 //----------------------------------------the following begins to calculate the water head of the inner nodes-----------------
 //----------------------------------------------------------------------------------------------------------------------------
 //

		c = 98000.0 / 12.0;   // c is pg/12u and the transferring of the unit
		c2 = 464.60533;


		for (i = np1n1; i <= nodes; i++)
		{ // 255内节点i                 

	  //-----------------------------------------------inner dry nodes and their matrix-------------------------------------------

			if (height[i] > h[i])
			{
				h[i] = 0.0;//eliminate the dry point

				newtonf(i) = h[i];//
							
				newtondif(i, i) = 1;
				newtonh0(i) = h[i];
			}

			//-----------------------------------------inner wet nodes and their matrix-------------------------------

			if (height[i] <= h[i]) 
			{ // available//1951



						// situ1: the inner nodes with some unsatured a neigboring nodes
				available1 = 0; available2 = 0;
				sumbl = 0.0;
				sumhbl = 0.0;
				sumqbl = 0.0;
				//errmax1=0.0 ;
				err1 = 0.0;
				NRe = 0;
				m = 0;
				z = 0;
				gan = 0;
				gan2 = 0;
				nnumber1 = 0;
				nnumber2 = 0;
				nnumber3 = 0;
				nnumber4 = 0;
				b1 = 0;
				b2 = 0;
				b3 = 0;
				b4 = 0;
				xl1 = 0;
				xl2 = 0;
				xl3 = 0;
				xl4 = 0;
				H1 = 0;
				H2 = 0;
				H3 = 0;
				H4 = 0;
				qj1 = 0;
				qj2 = 0;
				qj3 = 0;
				qj4 = 0;
				fenx1 = 0;
				fenx2 = 0;
				fenx0 = 0;
				aa1 = 0;
				bb1 = 0;
				cc1 = 0;
				for (j = 0; j < 5; j++)
				{
					wshi[j] = 0;
					cshi[j] = 0;
				}

				for (j = 1; j <= 4; j++)
				{
					nnumber = int(nijkl[i][j]);
					if (nnumber == 0)continue;
					if (h[nnumber] < height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
					{
						available1 = available1 + 1;//available1 = 1;
					}   //h[nnumber]=0;改这里也可以，从降雨补给的分配节点上判断available有几个
					if (h[nnumber] >= height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
					{
						available2 = available2 + 1;//
					}   

				}

				//--------------------------------------已经判断出来有一个、两个还是0个干邻点-----------------------------------
				//--------------------------------------------------------------------------------------------------------------




				//-----------------------------------------1 wet neighbour nodes------------------------------------------
				//--------------------------------------------------------------------------------------------------------
			  if (available2 == 1)
				{
					newtonh0(i) = h[i];
					for (j = 1; j <= 4; j++)
					{
						nnumber = int(nijkl[i][j]);
						if (nnumber == 0)continue;
						if (h[nnumber] >= height[nnumber])
						{
							newtonf(i) = h[i] - h[nnumber];
							for (kk = 1; kk <= np1n2; kk++)
							{
								if (kk == nnumber)
								{
									newtondif(i, kk) = -1;
								}
								
								else if (kk == i)
								{
									newtondif(i, kk) = 1;
								}
								else
								{
									newtondif(i, kk) = 0;
								}								
							}							
						}
					}
				}
				//------------------------------------2 or more wet neighbour nodes---------------------------------------------
			  if (available2 >= 2)
			  { 


				//-----------------------------------------有1个高位的干邻点--------------------------------------------
				if (available1 == 1)//中心节点是湿的，邻点是干的，有一个干节点就执行这个
				{ // available1
					for (j = 1; j <= 4; j++)
					{ //186
						nnumber = int(nijkl[i][j]);
						if (nnumber == 0)continue;
						if (height[nnumber] <= h[nnumber])//这个是邻节点周围的湿节点,只有找到湿节点才判断层流紊流
						{
							b = nijkl[i][j + 2 * maj];
							xl = nijkl[i][j + maj];
							J = (h[i] - h[nnumber]) / xl;
							kbl = 980.0 * b * b * b / 12.0 / 0.01;
							Vl = kbl * J / b;
							Re[i][nnumber] = 2 * b * Vl / 0.01;
							if (fabs(Re[i][nnumber]) > 2300)
							{
								NRe = NRe + 1; z = z + 1; wshi[z] = j;
							}
							else
							{
								m = m + 1; cshi[m] = j;
							}
							//sumbl=sumbl+b*b*b/xl;//sumhbl=sumhbl+h[nnumber]*b*b*b/xl;		
						}
					} //! 186


					for (j = 1; j <= 4; j++)
					{ //187
						nnumber = int(nijkl[i][j]);
						if (nnumber == 0)continue;
						if (height[nnumber] > h[nnumber])//如果是干节点
						{
							for (k = 1; k <= number_flow; k++)
								if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
								{
									//sumqbl=sumqbl+qab2[k][3]/c/sumbl;
									//Vl=qab2[k][3]/b;
									//Re[i][nnumber]=2*b*Vl/0.01;
									gan = k;
								}

						}
					} // 187


					if (NRe < 1)
					{
						//h[i]=sumhbl/sumbl ;
						fenn = 0;
						qj1 = int(cshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i]; b1 = 0; xl1 = 1; }

						qj2 = int(cshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i]; b2 = 0; xl2 = 1; }

						qj3 = int(cshi[3]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i]; b3 = 0; xl3 = 1; }

						qj4 = int(gan);

						

						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][3];//这点比较特殊，只要把流量在i节点的水均衡计算中表示出来就可以了


						newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];
						}

					}

					if (NRe == 1)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(cshi[1]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i]; b2 = 0; xl2 = 1; }

						qj3 = int(cshi[2]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i]; b3 = 0; xl3 = 1; }

						qj4 = int(gan);
						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][3];

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];
						}

					}
					if (NRe == 2)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(wshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

						qj3 = int(cshi[1]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i]; b3 = 0; xl3 = 1; }

						qj4 = int(gan);
						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][3];

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];

						}

					}

					if (NRe == 3)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(wshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

						qj3 = int(wshi[3]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

						qj4 = int(gan);
						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][3];

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];

						}

					}





				} // if available1>0


		//-----------------------------------------有2个高位的干邻点--------------------------------------------
		//---------------------------------有两个干节点的降雨情况-------------------------------------
		//--------------------------两个干节点的高度如果都低于中心节点高度，不接受降雨补给------------
				if (available1 == 2)
				{ // available1
					for (j = 1; j <= 4; j++)
					{ //186
						nnumber = int(nijkl[i][j]);
						if (nnumber == 0)continue;
						if (height[nnumber] <= h[nnumber])
						{
							b = nijkl[i][j + 2 * maj];
							xl = nijkl[i][j + maj];
							J = (h[i] - h[nnumber]) / xl;
							kbl = 980.0 * b * b * b / 12.0 / 0.01;
							Vl = kbl * J / b;
							Re[i][nnumber] = 2 * b * Vl / 0.01;
							if (fabs(Re[i][nnumber]) > 2300)
							{
								NRe = NRe + 1; z = z + 1; wshi[z] = j;
							}
							else
							{
								m = m + 1; cshi[m] = j;
							}
								
						}
					} // !186

					for (j = 1; j <= 4; j++)
					{ //187
						nnumber = int(nijkl[i][j]);
						if (nnumber == 0)continue;
						if (height[nnumber] > h[nnumber])//干节点
						{
							for (k = 1; k <= number_flow; k++)
								if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
								{
									
									gan = k;
								}
						}
					} //! 187

					for (j = 1; j <= 4; j++)
					{ //188
						nnumber = int(nijkl[i][j]);
						if (nnumber == 0)continue;//
						if (height[nnumber] > h[nnumber])//是干节点
						{
							for (k = 1; k <= number_flow; k++)
								if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
								{
									if (k != gan)
										gan2 = k;
								}
						}
					} //! 188
					


					if (NRe < 1)
					{
						//h[i]=sumhbl/sumbl ;
						fenn = 0;
						qj1 = int(cshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i]; b1 = 0; xl1 = 1; }

						qj2 = int(cshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i]; b2 = 0; xl2 = 1; }

						qj3 = int(gan);
						nnumber3 = int(qab2[qj3][2]);
						q3 = qab2[qj3][3];

						qj4 = int(gan2);
						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][3];

						newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];
						}

					}

					if (NRe == 1)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(cshi[1]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i]; b2 = 0; xl2 = 1; }

						qj3 = int(gan);
						nnumber3 = int(qab2[qj3][2]);
						q3 = qab2[qj3][3];

						qj4 = int(gan2);
						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][2];

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];
						}

					}
					if (NRe == 2)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(wshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

						qj3 = int(gan);
						nnumber3 = int(qab2[qj3][2]);
						q3 = qab2[qj3][3];

						qj4 = int(gan2);
						nnumber4 = int(qab2[qj4][2]);
						q4 = qab2[qj4][2];

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) - q3 - q4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = 0;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + 0 + 0;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];

						}

					}


				} // if available1 = 2,有两个干节点，完成。





		//------------------------------没有干节点，湿节点周围全都是湿节点-------------------

				if (available1 == 0)
				{ // available1
					for (j = 1; j <= maj; j++)
					{//15
						inode1 = i;
						inode2 = int(nijkl[inode1][j]);
						xl = nijkl[inode1][j + maj];
						b = nijkl[inode1][j + 2 * maj];
						if (b == 0.0) continue;
						if (height[inode2] <= h[inode2])//available为0代表的是没有降雨，但仍有可能有低于中心点的干节点，
							//如果不判断，这类干节点（所有干节点都低于中心点没有降雨补给的情况）被当做层流，会以水头0参与运算
						{
							J = (h[inode1] - h[inode2]) / xl;
							kbl = 980 * b * b * b / 12 / 0.01;
							Vl = kbl * J / b;
							Re[inode1][inode2] = 2 * b * Vl / 0.01;

							//m = m + 1; cshi[m] = j;
							if (fabs(Re[inode1][inode2]) > 2300)
							{
								NRe = NRe + 1; z = z + 1; wshi[z] = j;
							}
							else { m = m + 1; cshi[m] = j; }

						}//

					}//15 continue
					if (NRe < 1)
					{
						//h[i]=sumhbl/sumbl ;
						fenn = 0;
						qj1 = int(cshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);//np0 connected with np1，np2，np3，np4
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i]; b1 = 0; xl1 = 1; }

						qj2 = int(cshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i]; b2 = 0; xl2 = 1; }

						qj3 = int(cshi[3]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i]; b3 = 0; xl3 = 1; }

						qj4 = int(cshi[4]);
						if (qj4 > 0)
						{
							nnumber4 = int(nijkl[i][qj4]);
							b4 = nijkl[i][qj4 + 2 * maj];
							xl4 = nijkl[i][qj4 + maj];
							H4 = h[nnumber4];
						}
						else { H4 = h[i]; b4 = 0; xl4 = 1; }

						newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];
						}

					}

					if (NRe == 1)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(cshi[1]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i]; b2 = 0; xl2 = 1; }

						qj3 = int(cshi[2]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i]; b3 = 0; xl3 = 1; }

						qj4 = int(cshi[3]);
						if (qj4 > 0)
						{
							nnumber4 = int(nijkl[i][qj4]);
							b4 = nijkl[i][qj4 + 2 * maj];
							xl4 = nijkl[i][qj4 + maj];
							H4 = h[nnumber4];
						}
						else { H4 = h[i]; b4 = 0; xl4 = 1; }

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];
						}
					}

					if (NRe == 2)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(wshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

						qj3 = int(cshi[1]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i]; b3 = 0; xl3 = 1; }

						qj4 = int(cshi[2]);
						if (qj4 > 0)
						{
							nnumber4 = int(nijkl[i][qj4]);
							b4 = nijkl[i][qj4 + 2 * maj];
							xl4 = nijkl[i][qj4 + maj];
							H4 = h[nnumber4];
						}
						else { H4 = h[i]; b4 = 0; xl4 = 1; }

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];

						}

					}

					if (NRe == 3)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(wshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

						qj3 = int(wshi[3]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

						qj4 = int(cshi[1]);
						if (qj4 > 0)
						{
							nnumber4 = int(nijkl[i][qj4]);
							b4 = nijkl[i][qj4 + 2 * maj];
							xl4 = nijkl[i][qj4 + maj];
							H4 = h[nnumber4];
						}
						else { H4 = h[i]; b4 = 0; xl4 = 1; }

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c * b4 * b4 * b4 / xl4;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];

						}

					}

					if (NRe == 4)
					{
						fenn = 0;
						qj1 = int(wshi[1]);
						if (qj1 > 0)
						{
							nnumber1 = int(nijkl[i][qj1]);
							b1 = nijkl[i][qj1 + 2 * maj];
							xl1 = nijkl[i][qj1 + maj];
							H1 = h[nnumber1];
						}
						else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

						qj2 = int(wshi[2]);
						if (qj2 > 0)
						{
							nnumber2 = int(nijkl[i][qj2]);
							b2 = nijkl[i][qj2 + 2 * maj];
							xl2 = nijkl[i][qj2 + maj];
							H2 = h[nnumber2];
						}
						else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

						qj3 = int(wshi[3]);
						if (qj3 > 0)
						{
							nnumber3 = int(nijkl[i][qj3]);
							b3 = nijkl[i][qj3 + 2 * maj];
							xl3 = nijkl[i][qj3 + maj];
							H3 = h[nnumber3];
						}
						else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

						qj4 = int(wshi[4]);
						if (qj4 > 0)
						{
							nnumber4 = int(nijkl[i][qj4]);
							b4 = nijkl[i][qj4 + 2 * maj];
							xl4 = nijkl[i][qj4 + maj];
							H4 = h[nnumber4];
						}
						else { H4 = h[i] + 0.00001; b4 = 0; xl4 = 1; }

						newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c2 * pow(b4, 1.7143) * (h[i] - H4) / fabs(h[i] - H4) * pow(fabs((h[i] - H4) / xl4), 0.5714);
						for (kk = 1; kk <= np1n2; kk++)
						{
							if (kk == nnumber1)
							{
								newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
							}
							else if (kk == nnumber2)
							{
								newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
							}
							else if (kk == nnumber3)
							{
								newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
							}
							else if (kk == nnumber4)
							{
								newtondif(i, kk) = -c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
							}
							else if (kk == i)
							{
								newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
							}
							else
							{
								newtondif(i, kk) = 0;
							}

							newtonh0(i) = h[i];

						}

					}



				} // ! if available1=0

              }

			} // i node, water budget equation and matrix
			  //1951



		}  // 255inner node i，head vector，water budget vector，and matrix。
	  //-------------------------------------------------construction of matrix completed----------------------------------


	  //------------------------------matrix computing---------------------------------------------

		int iii;     iii = np1n2;
		mat newtondifA(iii, iii), newtondifB(iii, iii), newtondifB0(iii, iii);

		vec newtonh0A(iii);    vec newtonfA(iii);
		vec newtonh0B(iii);    vec newtontol(iii);


		for (i = 0; i <= (np1n2 -1); i++)
		{
			int jjj;  jjj = i + 1;
			newtonfA(i) = newtonf(jjj);
			newtonh0A(i) = newtonh0(jjj);
			for (j = 0; j <= (np1n2 - 1); j++)
			{
				int kkk;      kkk = j + 1;
				newtondifA(i, j) = newtondif(jjj, kkk);
			}
		}
		newtondifB0 = 100.0 * newtondifA;

		if (isnan(det(newtondifB0)))//
		{
			double deti2;
			deti2 = 0;
			deti2 = det(newtondifB0);
			//newtondifA0 = newtondifA;
			//newtontol = newtondifA0 * newtonfA; //solve(newtondifA, newtonfA);
		}
		else
		{
			double deti2;
			deti2 = 0;
			deti2 = det(newtondifB0);
			//if (fabs(deti2）
			newtondifB = inv(newtondifB0);
			newtontol = newtondifB * 100.0 * newtonfA;//solve(newtondifA,newtonfA);

			newtondifA0 = newtondifB * 100.0;
			double newtonfAAB[5000]; double newtonfAAC[5000];


			for (i = 0; i <= (np1n2 - 1); i++)
			{
				for (j = 0; j <= (np1n2 - 1); j++)
				{
					newtondifAA[i][j] = newtondifA0(i, j);
				}


				newtonfAAB[i] = newtonfA(i);

			}
			for (i = 0; i <= (np1n2 - 1); i++)
			{
				newtonfAAC[i] = newtonfAAB[i];
			}




		}

		if (isnan(det(newtondifB0)))//|| det(newtondifB0) == 0 
		{
			double deti2;
			deti2 = 0;
			deti2 = det(newtondifB0);
			//newtondifA0 = newtondifA;
			 //solve(newtondifA, newtonfA);

			double newtonfAA[5000]; double newtonfAB[5000];
			for (i = 0; i <= (np1n2 - 1); i++)
			{
				for (j = 0; j <= (np1n2 - 1); j++)
				{
					newtondifA0(i, j) = newtondifAA[i][j];
				}
				newtonfAA[i] = newtonfA(i);

			}
			for (i = 0; i <= (np1n2 - 1); i++)
			{
				newtonfAB[i] = newtonfAA[i];
			}
			newtontol = newtondifA0 * newtonfA;

		}

		//newtontol = solve(newtondifA,newtonfA);
		newtonh0B = newtonh0A - newtontol;

		for (i = 0; i <= (np1n2 - 1); i++)
		{
			if (fabs(newtontol(i)) > errmax) 	errmax = fabs(newtontol(i));
		}

		for (i = 0; i <= (np1n2 -1); i++)
		{
			h[i + 1] = newtonh0B(i);
		}




		//-----------------------------------------  water budget and error----------------------------
		for ( inode = 1; inode < np1n1; inode++ )
		{
			qbalance[inode] = 0;
			for (j = 1; j <= maj; j++)
			{
				inode1 = int(nijkl[inode][j]);
				Qcell[inode][j] = 0.0; Re[inode][inode1] = 0;
				if (fabs(h[inode]) < height[inode])continue;
				if (inode1 < 1)continue;
				if (inode < np1n1 && inode1 < np1n1)continue;
				if (h[inode1] >= height[inode1])//h[inode1]!= 0
				{
					b = nijkl[inode][2 * maj + j];
					xl = nijkl[inode][maj + j];
					J = (h[inode] - h[inode1]) / xl;
					kbl = 980.0 * b * b * b / 12.0 / 0.01;
					Vl = kbl * J / b;
					Re[inode][inode1] = 2 * b * Vl / 0.01;
					if (fabs(Re[inode][inode1]) < 2300) { Qcell[inode][j] = kbl * J; }
					else Qcell[inode][j] = 464.60533 * (h[inode] - h[inode1]) / fabs(h[inode] - h[inode1]) * pow(b, 1.7143) * pow(fabs(J), 0.5714);
				}
				else
				{
					for (k = 1; k <= number_flow; k++)
						if (inode == int(qab2[k][1]) && inode1 == int(qab2[k][2]))
							Qcell[inode][j] = -qab2[k][3];

				}
				qbalance[inode] = qbalance[inode] + Qcell[inode][j];
			}
		}

		for (inode = np1n1; inode <= nodes; inode++)
		{
			qbalance[inode] = 0;
			for (j = 1; j <= maj; j++)
			{
				inode1 = int(nijkl[inode][j]);
				Qcell[inode][j] = 0.0; Re[inode][inode1] = 0;//先全部置0，后对饱和区计算相应流量值。
				if (fabs(h[inode]) < height[inode])continue;
				if (inode1 < 1)continue;
				if (inode < np1n1 && inode1 < np1n1)continue;//当中心节点和邻点都是边界节点的时候才不执行循环，
				if (h[inode1] >= height[inode1])//h[inode1]!= 0，邻节点也是湿节点
				{
					b = nijkl[inode][2 * maj + j];
					xl = nijkl[inode][maj + j];
					J = (h[inode] - h[inode1]) / xl;
					kbl = 980.0 * b * b * b / 12.0 / 0.01;
					Vl = kbl * J / b;
					Re[inode][inode1] = 2 * b * Vl / 0.01;
					if (fabs(Re[inode][inode1]) < 2300) { Qcell[inode][j] = kbl * J; }// 980.0 * b * b * b / 12.0 / 0.01*(h[inode] - h[inode1]) / xl
					else Qcell[inode][j] = 464.60533 * (h[inode] - h[inode1]) / fabs(h[inode] - h[inode1]) * pow(b, 1.7143) * pow(fabs(J), 0.5714);
				}
				else
				{
					for (k = 1; k <= number_flow; k++)
						if (inode == int(qab2[k][1]) && inode1 == int(qab2[k][2]))
							Qcell[inode][j] = -qab2[k][3];

				}
				qbalance[inode] = qbalance[inode] + Qcell[inode][j];
			}
		}
		qqbudget = 0;

		for (inode = np1n1; inode <= nodes; inode++)
		{
		  qqbudget = qqbudget + qbalance[inode];
		}
		//----------------------------------------------------------------------------------------------------------------------


	} while ((errmax > tolmax && numit < 200) || (fabs(qqbudget) > 0.00001 && numit < 200));

	cout << " numit origin " << numit << " errmax00  " << errmax << " qqbudget00 " << qqbudget;


	double qqboundout0 = 0; double qqboundin0 = 0;
	for (i = 1; i <= npbf2; i++)
	{
		qqboundout0 = qqboundout0 + qbalance[i];//Qcell[i][maj]node
	}
	for (i = 1; i <= number_flow; i++)
	{
		qqboundin0 = qqboundin0 + qab2[i][3];

	}
	cout << " qqboundout0 " << qqboundout0 << " qqboundin0 " << qqboundin0 << " qoutbudget0  " << (qqboundout0 + qqboundin0) << endl;
   




	//---------------------------start----- 3 layers iterations of searching for the free surface----------------------------
	//------------------------------------------------------------   !   ----------------------------------------------------
	//------------------------------------------------------------   !   ----------------------------------------------------


	
	int drytest; 
	 numitouter = 0;
	 

do           
//-------------------------------------------------------- outer layer iteration-----------------------------------------
//------------------------------------------------------------   !   ----------------------------------------------------
{
	          inode0 = 0; inode = 0; inode1 = 0; inode2 = 0; inode3 = 0; inode4 = 0; inode5 = 0;
	          if(numitouter==0) wettest01 = 0;

			  if (numitouter > 0 && numitouter <= 100)

			  {


				  inode2 = 0; 	  int layernodei = 0;	  wettest01 = 0;
				  int layernode[1000];	  for (i = 0; i <= 999; i++)  layernode[i] = 0;




				  for (inode = np1n1; inode <= nodes; inode++)
				  {


					  int wetnode = 0; //int wetnode2 = 0;
					  for (j = 1; j <= maj; j++)/////////&& (fabs(h[inode1]) < 0.001)    (height[inode1]-2.0)
					  {
						  inode1 = int(nijkl[inode][j]);//  np1n1												
						  if ((inode1 >= 1) && (height[inode] > h[inode]) && h[inode1] >= (height[inode] + 0))
						  {
							  wetnode = wetnode + 1;
						  }
					  }

					  if (  wetnode >= 1 )       
					  {
						  layernodei = layernodei + 1;
						  layernode[layernodei] = inode;
					  }


				  }





				  int ii = 0; int iii = 0; int tempii = 0; int iijj = 0;


				  for (i = 1; i <= layernodei; i++)
				  {
					  for (ii = i; ii <= layernodei; ii++)
					  {
						  inode1 = layernode[ii];
						  if ((inode1 >= npbf1) && (inode1 <= npbf2))
						  {
							  iii = ii;
							  tempii = layernode[i];
							  layernode[i] = layernode[iii];
							  layernode[iii] = tempii;
						  }
					  }
				  }

				  //if ( (numitouter %2) ==0 )
				  for (i = 1; i <= layernodei; i++)
				  {
					  if (layernode[i] < np1n1) continue;
					  double dertah = 0;
					  for (ii = i; ii <= layernodei; ii++)
					  {
						  inode = layernode[ii];
						  for (iijj = 1; iijj <= maj; iijj++)/////////&& (fabs(h[inode1]) < 0.001)    (height[inode1]-2.0)
						  {
							  inode1 = int(nijkl[inode][iijj]);// 												
							  if ((inode1 >= 1) && (height[inode] > h[inode]) && h[inode1] >= (height[inode] + 0))
							  								 
								  if (dertah < ( h[inode1]-height[inode]) )
								  {
									  dertah = h[inode1] - height[inode];
									  iii = ii;
								  }
							  
						  }

					  }

					  tempii = layernode[i];
					  layernode[i] = layernode[iii];
					  layernode[iii] = tempii;
				  }

                      cout << " layernodei   " << layernodei << endl;



				  int inodeup; //try uplifting water head node




				  //-------------------------------------------------------- midlle layer iteration----------------------------------------
				  //------------------------------------------------------------   !   ----------------------------------------------------
				  ////////////////////////开始计算一次水面的上升，一次水面的上升包括多个点的多次水头计算，有些点升上去有些点升不上去////////
				  //////////////////////////////按水头差从大到小排序，然后依次升水头////////////////////////////////////////////////////////				  


				  for (inodeup = 1; inodeup <= layernodei; inodeup++)
				  {


					  /////////////////////////////升水头之前，判断出潜水面上的湿节点///////////////////////////////////////////////////////
					  
					  double layerwetnode[5000][2]; int layerdrynode[5000];
					  for (i = 0; i <= 4999; i++)
					  {
						  layerwetnode[i][0] = 0;
						  layerwetnode[i][1] = 0;
						  layerdrynode[i] = 0;
					  }
					  int layerwetnodei = 0; int layerdrynodei = 0;
					  for (inode = 1; inode <= nodes; inode++)//np1n1
					  {

						  if (h[inode] >= height[inode])////////不只是潜水面，原来所有节点都有可能受到影响，潜水面下面也可能
						  {
							  layerwetnodei = layerwetnodei + 1;
							  layerwetnode[layerwetnodei][0] = inode;
							  layerwetnode[layerwetnodei][1] = h[inode];
						  }
						  else
						  {
							  layerdrynodei = layerdrynodei + 1;
							  layerdrynode[layerdrynodei] = inode;

						  }
					  }


					  inode = layernode[inodeup];

					  int layernodeA = 0;

					  int wetnode = 0;
					  for (j = 1; j <= maj; j++)/////////&& (fabs(h[inode1]) < 0.001)    (height[inode1]-2.0)
					  {
						  inode1 = int(nijkl[inode][j]);// 												
						  if ((inode1 >= 1 ) && (height[inode] > h[inode]) && h[inode1] >= ( height[inode] + 0) )//&& height[inode] > height[inode1]
							  wetnode = wetnode + 1;
					  }

					  numit = 0;

					if (wetnode >= 1)
					{
						for (j = 1; j <= maj; j++)
						{
							inode1 = int(nijkl[inode][j]);// 												
							if ((inode1 >= 1 ) && (height[inode] > h[inode]) && h[inode1] >= (height[inode] + 0))

							{							
								if (inode >= np1n1) h[inode] = h[inode1];
								else  h[inode] = height[inode];
							}
							   
						}
						
					  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						  layernodeA = inode;



						  //------------------------------------  inner iteration of confined water table,  first time------------------------------------------
						  //------------------------------------------------------------   inner core   --------------------------------------------------------

						  numit = 0;	      newtonh0.fill(0);      newtonf.fill(0);       newtondif.fill(0);	            newtondifA0.fill(0);
						  do  
						  {                                

							  errmax1 = 0.0;
							  errmax2 = 0.0;
							  errmax = 0.0;
							  numit = numit + 1;

						  //-------------------------- Initiation for Water head calculation --------------------------

							  int n, n1; n = 0; n1 = 0;
							  h[0] = 0.0;

							  for (i = 1; i <= nboun1; i++)
							  {//81
								  n1 = icount[i];
								  for (j = 1; j <= n1; j++)
								  {//73
									  n = n + 1;
									  h[n] = hb1[i];
									  newtonf(n) = h[n] - hb1[i];
									  newtondif(n, n) = 1;
									  newtonh0(n) = h[n];
								  }//73 
							  }//81 



							  //---------------the following deals with the second boundary----------------------------------------
							  ib2 = 0;
							  inode = npb21 - 1;
							  inode1 = inode;

							  for (nb2 = nboun1 + 1; nb2 <= nboun1 + nboun2; nb2++)
							  {//620
								  ib2 = ib2 + 1;
								  allkb = 0.0;
								  for (j = 1; j <= icount[nb2]; j++)
								  {//605
									  inode1 = inode1 + 1;
									  inode2 = int(nijkl[inode1][maj]);
									  b3 = nijkl[inode1][3 * maj];
									  
									  kbl = b3 * b3 * b3;

									  allkb = allkb + kbl;
								  }//605	continue

								  for (k = 1; k <= icount[nb2]; k++)
								  {//610
									  inode = inode + 1;
									  inode2 = int(nijkl[inode][maj]);
									  b3 = nijkl[inode][3 * maj];
									  xl = nijkl[inode][2 * maj];
									  qq = qb2[ib2] * b3 * b3 * b3 / allkb; 

									  if (qq == 0)
									  {
										  h[inode] = h[inode2];
										  newtonf(inode) = h[inode] - h[inode2];
										  for (kk = 1; kk <= np1n2; kk++)
										  {
											  if (kk == inode)
											  {
												  newtondif(inode, kk) = 1;
											  }
											  else if (kk == inode2)
											  {
												  newtondif(inode, kk) = -1;
											  }
											  else
											  {
												  newtondif(inode, kk) = 0;
											  }

											  newtonh0(inode) = h[inode];
										  }
									  }
									  else
									  {
										  Vl = fabs(qq / b3);
										  reynold = b3 * Vl * 2.0 / 0.01;
										  if (reynold < 2300)
										  {
											  J = qq * 12 * 0.01 / b3 / b3 / b3 / 980.0;
											  h[inode] = h[inode2] + J * xl;
										  }
										  else {
											  qq1 = fabs(qq) / 4.70 / b3;
											  J = pow(qq1, 7.0) * 0.01 / b3 / b3 / b3 / b3 / b3 / 980.0 / 980.0 / 980.0 / 980.0;
											  J = pow(J, 0.25) * qq / fabs(qq);
											  h[inode] = h[inode2] + J * xl;
										  }


										  newtonf(inode) = h[inode] - h[inode2] - J * xl;
										  for (kk = 1; kk <= np1n2; kk++)
										  {
											  if (kk == inode)
											  {
												  newtondif(inode, kk) = 1;
											  }
											  else if (kk == inode2)
											  {
												  newtondif(inode, kk) = -1;
											  }
											  else
											  {
												  newtondif(inode, kk) = 0;
											  }

											  newtonh0(inode) = h[inode];
										  }
									  }
								  }//610 continue
							  }//620	continue

						  //--------------------the following deals with the third boundary condition-------------------------------------------------

							  for (i = npbf1; i <= npbf2; i++)
							  {
								  h[i] = 0.0;
								  newtonf(i) = h[i];
								  newtondif(i, i) = 1;
								  newtonh0(i) = h[i];
							  }


						  //--------------------------the following deals with the spring boundary--------------------------------------------------

							  Qspring = 0.0;
							  for (i = np1n1; i <= nodes; i++)
								  for (j = 1; j <= 4; j++)
								  {
									  nnumber = int(nijkl[i][j]);
									  b = nijkl[i][j + 2 * maj];
									  xl = nijkl[i][j + maj];
									  if (fabs(xynp[nnumber][1]) > 0.01 && xynp[nnumber][1] < (maxx - 0.01))continue;
									 
									  if (nnumber > npb22 && nnumber < np1n1)
									  {
										  J = (h[i] - height[nnumber]) / xl;
										  kbl = 980.0 * b * b * b / 12.0 / 0.01;
										  Vl = kbl * J / b;
										  Re[i][nnumber] = 2 * b * Vl / 0.01;
										  if (Re[i][nnumber] > 2300) { Qspring = 464.60533 * pow(b, 1.7143) * pow(fabs(J), 0.5714); }
										  else Qspring = kbl * J;
										  // when Qsrping is greater than 0 , it means it can be regared as a spring
										  if (h[i] < (height[nnumber]))
										  {
											  h[nnumber] = 0;
											  newtonf(nnumber) = h[nnumber];
											  newtondif(nnumber, nnumber) = 1;
											  newtonh0(nnumber) = h[nnumber];
										  }

										  if (h[i] >= (height[nnumber]))
										  {
											  h[nnumber] = height[nnumber];
											  newtonf(nnumber) = h[nnumber] - height[nnumber];
											  newtondif(nnumber, nnumber) = 1;
											  newtonh0(nnumber) = h[nnumber];
										  }
									  }
								  }



							 
							  //----------------------the following deals with the rainfall recharge on water table-------------------------------

							  for (i = 0; i < 500; i++)
								  for (j = 0; j < 4; j++)  qab2[i][j] = 0;


							  number_flow = 0;
							  for (i = np1n1; i <= nodes; i++)
							  {
								  nnumber1 = 0; nnumber2 = 0;
								  if (height[i] <= h[i]) 
								  {
									  for (j = 1; j <= 4; j++)
									  {
										  nnumber = int(nijkl[i][j]);
										  if (nnumber == 0)continue;
										  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

										  if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
										  {
											  nnumber1 = nnumber1 + 1;
										  }

										  if (h[nnumber] >= height[nnumber])
										  {
											  nnumber2 = nnumber2 + 1;
										  }
									  }

									  if (nnumber1 != 0 && nnumber2 >= 2)
										  for (j = 1; j <= 4; j++)
										  {
											  nnumber = int(nijkl[i][j]);
											  if (nnumber == 0)continue;
											  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

											  if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
											  {
												  number_flow = number_flow + 1;
												  qab2[number_flow][1] = i;
												  qab2[number_flow][2] = nnumber;
											  }

										  }

									  
								  }	
							  }



							  //-----------------the following calculates the flow rate of each fracture segment-------------
							  tran_coeff = 10.0 / 365.0 / 24.0 / 60.0 / 60.0;  
									   //the unit of Qsurfcorro is mm/a, and the unit of Qsurfcorro1 is ml/s
									   //tran_coeff is the transferring coefficient.
							  Qsurfcorro1 = tran_coeff * (Qsurfcorro * (maxx - 0.0) * 1.0);//*(maxx-minx)

							  sumwidth = 0;
							  ninfil = 0;
							  for (i = 1; i <= number_flow; i++)
							  {
								  inode = int(qab2[i][1]);
								  for (k = 1; k <= 4; k++)
								  {
									  nnumber = int(nijkl[inode][k]);
									 
									  if (nnumber == int(qab2[i][2]))
									  {
										  sumwidth = sumwidth + nijkl[inode][k + 2 * maj];
										  ninfil++;
									  }
								  }
							  }
							  for (i = 1; i <= number_flow; i++)
							  {
								  inode = int(qab2[i][1]);
								  for (k = 1; k <= 4; k++)
								  {
									  nnumber = int(nijkl[inode][k]);
									  if (nnumber == int(qab2[i][2]))//&& h[inode]<maxy
									  {
										  
										  qab2[i][3] = Qsurfcorro1 / ninfil;//evenly distribute infiltration																	
									  }
									  
								  }

							  }

							 



					   //-------------------------------------------------inner node matrix------------------------------------
					   //------------------------------------------------------------------------------------------------------
					   //the following begins to calculate the water head of the inner nodes


							  c = 98000.0 / 12.0;   // c is pg/12u and the transferring of the unit
							  c2 = 464.60533;

							  for (i = np1n1; i <= nodes; i++)
							  { // 255内节点i                 

							//-----------------------------------------------inner dry nodes and their matrix-------------------------------------------

								  if (height[i] > h[i])
								  {
									  h[i] = 0.0;//eliminate the dry point

									  newtonf(i) = h[i];//

									  newtondif(i, i) = 1;
									  newtonh0(i) = h[i];
								  }

								  //-----------------------------------------inner wet nodes and their matrix-------------------------------

								  if (height[i] <= h[i])
								  { // available//1951



											  // situ1: the inner nodes with some unsatured a neigboring nodes
									  available1 = 0; available2 = 0;
									  sumbl = 0.0;
									  sumhbl = 0.0;
									  sumqbl = 0.0;
									  //errmax1=0.0 ;
									  err1 = 0.0;
									  NRe = 0;
									  m = 0;
									  z = 0;
									  gan = 0;
									  gan2 = 0;
									  nnumber1 = 0;
									  nnumber2 = 0;
									  nnumber3 = 0;
									  nnumber4 = 0;
									  b1 = 0;
									  b2 = 0;
									  b3 = 0;
									  b4 = 0;
									  xl1 = 0;
									  xl2 = 0;
									  xl3 = 0;
									  xl4 = 0;
									  H1 = 0;
									  H2 = 0;
									  H3 = 0;
									  H4 = 0;
									  qj1 = 0;
									  qj2 = 0;
									  qj3 = 0;
									  qj4 = 0;
									  fenx1 = 0;
									  fenx2 = 0;
									  fenx0 = 0;
									  aa1 = 0;
									  bb1 = 0;
									  cc1 = 0;
									  for (j = 0; j < 5; j++)
									  {
										  wshi[j] = 0;
										  cshi[j] = 0;
									  }

									  for (j = 1; j <= 4; j++)
									  {
										  nnumber = int(nijkl[i][j]);
										  if (nnumber == 0)continue;
										  if (h[nnumber] < height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
										  {
											  available1 = available1 + 1;//available1 = 1;
										  }   //h[nnumber]=0;改这里也可以，从降雨补给的分配节点上判断available有几个
										  if (h[nnumber] >= height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
										  {
											  available2 = available2 + 1;//
										  }

									  }

									  //--------------------------------------已经判断出来有一个、两个还是0个干邻点-----------------------------------
									  //--------------------------------------------------------------------------------------------------------------




									  //-----------------------------------------1 wet neighbour nodes------------------------------------------
									  //--------------------------------------------------------------------------------------------------------
									  if (available2 == 1)
									  {
										  newtonh0(i) = h[i];
										  for (j = 1; j <= 4; j++)
										  {
											  nnumber = int(nijkl[i][j]);
											  if (nnumber == 0)continue;
											  if (h[nnumber] >= height[nnumber])
											  {
												  newtonf(i) = h[i] - h[nnumber];
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber)
													  {
														  newtondif(i, kk) = -1;
													  }

													  else if (kk == i)
													  {
														  newtondif(i, kk) = 1;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }
												  }
											  }
										  }
									  }
									  //------------------------------------2 or more wet neighbour nodes---------------------------------------------
									  if (available2 >= 2)
									  {


										  //-----------------------------------------有1个高位的干邻点--------------------------------------------
										  if (available1 == 1)//中心节点是湿的，邻点是干的，有一个干节点就执行这个
										  { // available1
											  for (j = 1; j <= 4; j++)
											  { //186
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] <= h[nnumber])//这个是邻节点周围的湿节点,只有找到湿节点才判断层流紊流
												  {
													  b = nijkl[i][j + 2 * maj];
													  xl = nijkl[i][j + maj];
													  J = (h[i] - h[nnumber]) / xl;
													  kbl = 980.0 * b * b * b / 12.0 / 0.01;
													  Vl = kbl * J / b;
													  Re[i][nnumber] = 2 * b * Vl / 0.01;
													  if (fabs(Re[i][nnumber]) > 2300)
													  {
														  NRe = NRe + 1; z = z + 1; wshi[z] = j;
													  }
													  else
													  {
														  m = m + 1; cshi[m] = j;
													  }
													  //sumbl=sumbl+b*b*b/xl;//sumhbl=sumhbl+h[nnumber]*b*b*b/xl;		
												  }
											  } //! 186


											  for (j = 1; j <= 4; j++)
											  { //187
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] > h[nnumber])//如果是干节点
												  {
													  for (k = 1; k <= number_flow; k++)
														  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
														  {
															  //sumqbl=sumqbl+qab2[k][3]/c/sumbl;
															  //Vl=qab2[k][3]/b;
															  //Re[i][nnumber]=2*b*Vl/0.01;
															  gan = k;
														  }

												  }
											  } // 187


											  if (NRe < 1)
											  {
												  //h[i]=sumhbl/sumbl ;
												  fenn = 0;
												  qj1 = int(cshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i]; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);



												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];//这点比较特殊，只要把流量在i节点的水均衡计算中表示出来就可以了


												  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }

											  if (NRe == 1)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[1]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[2]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }
											  if (NRe == 2)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[1]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }

											  if (NRe == 3)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(wshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }





										  } // if available1>0


								  //-----------------------------------------有2个高位的干邻点--------------------------------------------
								  //---------------------------------有两个干节点的降雨情况-------------------------------------
								  //--------------------------两个干节点的高度如果都低于中心节点高度，不接受降雨补给------------
										  if (available1 == 2)
										  { // available1
											  for (j = 1; j <= 4; j++)
											  { //186
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] <= h[nnumber])
												  {
													  b = nijkl[i][j + 2 * maj];
													  xl = nijkl[i][j + maj];
													  J = (h[i] - h[nnumber]) / xl;
													  kbl = 980.0 * b * b * b / 12.0 / 0.01;
													  Vl = kbl * J / b;
													  Re[i][nnumber] = 2 * b * Vl / 0.01;
													  if (fabs(Re[i][nnumber]) > 2300)
													  {
														  NRe = NRe + 1; z = z + 1; wshi[z] = j;
													  }
													  else
													  {
														  m = m + 1; cshi[m] = j;
													  }

												  }
											  } // !186

											  for (j = 1; j <= 4; j++)
											  { //187
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] > h[nnumber])//干节点
												  {
													  for (k = 1; k <= number_flow; k++)
														  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
														  {

															  gan = k;
														  }
												  }
											  } //! 187

											  for (j = 1; j <= 4; j++)
											  { //188
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;//
												  if (height[nnumber] > h[nnumber])//是干节点
												  {
													  for (k = 1; k <= number_flow; k++)
														  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
														  {
															  if (k != gan)
																  gan2 = k;
														  }
												  }
											  } //! 188



											  if (NRe < 1)
											  {
												  //h[i]=sumhbl/sumbl ;
												  fenn = 0;
												  qj1 = int(cshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i]; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(gan);
												  nnumber3 = int(qab2[qj3][2]);
												  q3 = qab2[qj3][3];

												  qj4 = int(gan2);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }

											  if (NRe == 1)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[1]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(gan);
												  nnumber3 = int(qab2[qj3][2]);
												  q3 = qab2[qj3][3];

												  qj4 = int(gan2);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][2];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }
											  if (NRe == 2)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(gan);
												  nnumber3 = int(qab2[qj3][2]);
												  q3 = qab2[qj3][3];

												  qj4 = int(gan2);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][2];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) - q3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + 0 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }


										  } // if available1 = 2,有两个干节点，完成。





								  //------------------------------没有干节点，湿节点周围全都是湿节点-------------------

										  if (available1 == 0)
										  { // available1
											  for (j = 1; j <= maj; j++)
											  {//15
												  inode1 = i;
												  inode2 = int(nijkl[inode1][j]);
												  xl = nijkl[inode1][j + maj];
												  b = nijkl[inode1][j + 2 * maj];
												  if (b == 0.0) continue;
												  if (height[inode2] <= h[inode2])//available为0代表的是没有降雨，但仍有可能有低于中心点的干节点，
													  //如果不判断，这类干节点（所有干节点都低于中心点没有降雨补给的情况）被当做层流，会以水头0参与运算
												  {
													  J = (h[inode1] - h[inode2]) / xl;
													  kbl = 980 * b * b * b / 12 / 0.01;
													  Vl = kbl * J / b;
													  Re[inode1][inode2] = 2 * b * Vl / 0.01;

													  //m = m + 1; cshi[m] = j;
													  if (fabs(Re[inode1][inode2]) > 2300)
													  {
														  NRe = NRe + 1; z = z + 1; wshi[z] = j;
													  }
													  else { m = m + 1; cshi[m] = j; }

												  }//

											  }//15 continue
											  if (NRe < 1)
											  {
												  //h[i]=sumhbl/sumbl ;
												  fenn = 0;
												  qj1 = int(cshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);//np0 connected with np1，np2，np3，np4
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i]; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[4]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }

											  if (NRe == 1)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[1]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[2]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[3]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }
											  }

											  if (NRe == 2)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[1]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[2]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }

											  if (NRe == 3)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(wshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[1]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }

											  if (NRe == 4)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(wshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

												  qj4 = int(wshi[4]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i] + 0.00001; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c2 * pow(b4, 1.7143) * (h[i] - H4) / fabs(h[i] - H4) * pow(fabs((h[i] - H4) / xl4), 0.5714);
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }



										  } // ! if available1=0

									  }

								  } // i node, water budget equation and matrix
									//1951



							  }  // 255inner node i，head vector，water budget vector，and matrix。


					   //-------------------------------------------------inner node matrix completed----------------------------
					   //--------------------------------------------------------------------------------------------------------
					   
					   
					   
					   //------------------------------matrix computing---------------------------------------------

							  int iii;     iii = np1n2 ;
							  mat newtondifA(iii, iii), newtondifB(iii, iii), newtondifB0(iii, iii);

							  vec newtonh0A(iii);    vec newtonfA(iii);
							  vec newtonh0B(iii);    vec newtontol(iii);


							  for (i = 0; i <= (np1n2 - 1); i++)
							  {
								  int jjj;  jjj = i + 1;
								  newtonfA(i) = newtonf(jjj);
								  newtonh0A(i) = newtonh0(jjj);
								  for (j = 0; j <= (np1n2 - 1); j++)
								  {
									  int kkk;      kkk = j + 1;
									  newtondifA(i, j) = newtondif(jjj, kkk);
								  }
							  }
							  newtondifB0 = 1.0 * newtondifA;




							  {
								  
								  try
								  {
									  newtondifB = inv(newtondifB0);
								  }

								  catch (...)
								  {
									  cout << " error01catch ";
									  numit = 30;
									  newtondifB = newtondifB0;
									

								  }
								  newtontol = newtondifB * 1.0 * newtonfA;

							  }


							  //newtontol = solve(newtondifA,newtonfA);
							  newtonh0B = newtonh0A - newtontol;

							  for (i = 0; i <= (np1n2 - 1); i++)
							  {
								  if (fabs(newtontol(i)) > errmax) 	errmax = fabs(newtontol(i));
							  }

							  for (i = npb21-1; i <= (np1n2 - 1); i++)
							  {
								  h[i + 1] = newtonh0B(i);
							  }
							  




								  



							  //-----------------------------------------  water budget and error  ----------------------------


							  for (inode = 1; inode <= nodes; inode++)
							  {
								  qbalance[inode] = 0;
								  for (j = 1; j <= maj; j++)
								  {
									  inode1 = int(nijkl[inode][j]);
									  Qcell[inode][j] = 0.0; Re[inode][inode1] = 0;
									  if (fabs(h[inode]) < height[inode])continue;

									  if (inode1 < 1)continue;
									  if (inode < np1n1 && inode1 < np1n1)continue;

									  if (h[inode1] >= height[inode1])//h[inode1]!= 0
									  {
										  b = nijkl[inode][2 * maj + j];
										  xl = nijkl[inode][maj + j];
										  J = (h[inode] - h[inode1]) / xl;
										  kbl = 980.0 * b * b * b / 12.0 / 0.01;
										  Vl = kbl * J / b;
										  Re[inode][inode1] = 2 * b * Vl / 0.01;
										  if (fabs(Re[inode][inode1]) < 2300) { Qcell[inode][j] = kbl * J; }
										  else Qcell[inode][j] = 464.60533 * (h[inode] - h[inode1]) / fabs(h[inode] - h[inode1]) * pow(b, 1.7143) * pow(fabs(J), 0.5714);
									  }
									  else
									  {
										  for (k = 1; k <= number_flow; k++)
											  if (inode == int(qab2[k][1]) && inode1 == int(qab2[k][2]))
												  Qcell[inode][j] = -qab2[k][3];

									  }

									  qbalance[inode] = qbalance[inode] + Qcell[inode][j];
								  }
							  }
							  qqbudget = 0;

							  for (inode = np1n1; inode <= nodes; inode++)
	                         	{
	                         	  qqbudget = qqbudget + qbalance[inode];
	                         	}



						  } while ((errmax > tolmax && numit < 30) || (fabs(qqbudget) > 0.0001 && numit < 30));
						  
						  cout << " numit01   " << numit << " errmax01  " << errmax << "  qqbudget01  " << qqbudget ;

						  //-------------------------------------------------------------  inner iteration  ----------------------------------------------------
						  //------------------------------------------------------------   inner core   --------------------------------------------------------




						  double qqboundout1 = 0; double qqboundin1 = 0;
						  for (i = 1; i <= npbf2; i++)
						  {
							  qqboundout1 = qqboundout1 + qbalance[i];//Qcell[i][maj]node
						  }
						  for (i = 1; i <= number_flow; i++)
						  {
							  qqboundin1 = qqboundin1 + qab2[i][3];

						  }
						  cout << " qqboundout1 " << qqboundout1 << " qqboundin1 " << qqboundin1 << " qoutbudget1  " << (qqboundout1 + qqboundin1) << endl;


 

					  if (numit > 29)
					  {
						  if (layernodeA > 0)  h[layernodeA] = 0;

					  }

					  int layerwet2dry = 0; int layerdry2wet = 0;
					  for (i = 1; i <= layerwetnodei; i++)//原来是湿点如果变干了，不仅要升的干点不升了，还要把这些原来的湿点再恢复过来
					  {
						  
						  inode = int(layerwetnode[i][0]);
						  if (height[inode] > h[inode])
						  {
							  //h[inode] = layerwetnode[i][1];
							  layerwet2dry = layerwet2dry + 1;//
						  }
					  }
					  for (i = 1; i <= layerdrynodei; i++)//
					  {
						  inode = layerdrynode[i];                 //npbf1
						  if (height[inode] <= h[inode] &&  inode >= np1n1 && layernodeA != inode && layernodeA > 0)//变湿节点判断的时候要去掉渗出点，渗出点是有可能变湿的
						  {
							  
							  layerdry2wet = layerdry2wet + 1;
						  }
					  }
					  ////////////////如果造成别的湿节点变干,就把升上去的那个点再变干，恢复原来的形状/////////////////

					  if ((layerwet2dry > 0) && (layernodeA > 0))  h[layernodeA] = 0;
					  if ((layerdry2wet > 0) && (layernodeA > 0))  h[layernodeA] = 0;


					  ////////////////////////////////////下面是规定潜水面的发展形状，不能在潜水面以上还有湿节点//////////////////////////////////////////////
					  if (h[layernodeA] >= height[layernodeA]  && (layernodeA > 0))//&& h[layernodeA] <= (height[layernodeA]+10)
					  {

						  int jj4;
						  for (jj4 = 1; jj4 <= maj; jj4++)
						  {
							  inode4 = int(nijkl[layernodeA][jj4]);

							  double wetdryleft04=0; double wetdryright04=0;
							  if (h[inode4] >= height[inode4] && inode4 > 0 && xynp[inode4][1] >= xynp[layernodeA][1])
							  {
								  wetdryleft04 = xynp[layernodeA][1]; wetdryright04 = xynp[inode4][1];
							  }
							  if (h[inode4] >= height[inode4] && inode4 > 0 && xynp[inode4][1] < xynp[layernodeA][1])
							  {
								  wetdryleft04 = xynp[inode4][1];    wetdryright04 = xynp[layernodeA][1];
							  }

							

							  if (h[inode4] >= height[inode4] && inode4 > 0)
							  for (inode2 = np1n1; inode2 <= nodes; inode2++)
							  {           //npbf1
								  if (h[inode2] < height[inode2] && height[inode4] > height[inode2] && height[layernodeA] > height[inode2] && xynp[inode2][1] >= wetdryleft04 && xynp[inode2][1] <= wetdryright04)
									  h[layernodeA] = 0;
								  if (h[inode2] < height[inode2] && height[inode4] > height[inode2] && height[layernodeA] <= height[inode2] && xynp[inode2][1] >= wetdryleft04 && xynp[inode2][1] <= wetdryright04)
									  if(  fabs(xynp[inode4][2] - xynp[layernodeA][2]) / fabs(xynp[inode4][1] - xynp[layernodeA][1]) > fabs(xynp[inode2][2] - xynp[layernodeA][2]) / fabs(xynp[inode2][1] - xynp[layernodeA][1]))
									  h[layernodeA] = 0;

								  if (h[inode2] < height[inode2] && height[inode4] <= height[inode2] && height[layernodeA] > height[inode2] && xynp[inode2][1] >= wetdryleft04 && xynp[inode2][1] <= wetdryright04)
									  if (fabs(xynp[inode4][2] - xynp[layernodeA][2]) / fabs(xynp[inode4][1] - xynp[layernodeA][1]) > fabs(xynp[inode4][2] - xynp[inode2][2]) / fabs(xynp[inode4][1] - xynp[inode2][1]))
									  h[layernodeA] = 0;

							  }

							  if (h[inode4] >= height[inode4])///////解决湿裂隙下方没有干点，却有干裂隙的问题
								  for (inode2 = np1n1; inode2 <= nodes; inode2++)
								  {
									  if (h[inode2] >= height[inode2] && height[inode4] > height[inode2] && height[layernodeA] > height[inode2] && xynp[inode2][1] > wetdryleft04 && xynp[inode2][1] < wetdryright04)
										  for (j = 1; j <= maj; j++)
										  {
											  inode3 = int(nijkl[inode2][j]);//inode3就是干点，无论在左边还是右边
											  if (inode3 >= npbf1 && height[inode3] > h[inode3] && height[inode3] < height[inode2])
												  h[layernodeA] = 0;
										  }


									  if (h[inode2] >= height[inode2] && height[inode4] > height[inode2] && height[layernodeA] <= height[inode2] && xynp[inode2][1] > wetdryleft04 && xynp[inode2][1] < wetdryright04)
										  if (fabs(xynp[inode4][2] - xynp[layernodeA][2]) / fabs(xynp[inode4][1] - xynp[layernodeA][1]) >= fabs(xynp[inode2][2] - xynp[layernodeA][2]) / fabs(xynp[inode2][1] - xynp[layernodeA][1]))
											  for (j = 1; j <= maj; j++)
											  {
												  inode3 = int(nijkl[inode2][j]);//inode3就是干点，无论在左边还是右边
												  if (inode3 >= npbf1 && height[inode3] > h[inode3] && height[inode3] < height[layernodeA])
													  h[layernodeA] = 0;
											  }



									  if (h[inode2] >= height[inode2] && height[inode4] <= height[inode2] && height[layernodeA] > height[inode2] && xynp[inode2][1] > wetdryleft04 && xynp[inode2][1] < wetdryright04)
										  if (fabs(xynp[inode4][2] - xynp[layernodeA][2]) / fabs(xynp[inode4][1] - xynp[layernodeA][1]) >= fabs(xynp[inode4][2] - xynp[inode2][2]) / fabs(xynp[inode4][1] - xynp[inode2][1]))
											  for (j = 1; j <= maj; j++)
											  {
												  inode3 = int(nijkl[inode2][j]);//inode3就是干点，无论在左边还是右边，，因此要
												  if (inode3 >= npbf1 && height[inode3] > h[inode3] && height[inode3] < height[inode4])
													  h[layernodeA] = 0;
											  }

								  }

						  }
					  }
					  ////////////////////////////////////上面是规定潜水面的发展形状，在潜水面以上不能有湿节点，如果有就不升此次水头/////////////////////////






					  if ((h[layernodeA] == 0 && layernodeA > 0) || (numit == 30))
					  {

						  for (i = 1; i <= layerwetnodei; i++)
						  {
							  inode = int(layerwetnode[i][0]);
							  h[inode] = layerwetnode[i][1];
						  }

						  for (i = 1; i <= layerdrynodei; i++)
						  {
							  inode = layerdrynode[i];
							  h[inode] = 0;
						  }
						  


						  //------------------------------------  inner iteration of confined water table,  second time ----------------------------------------
						  //------------------------------------------------------------   inner core   --------------------------------------------------------


						  numit = 0;	      newtonh0.fill(0);      newtonf.fill(0);       newtondif.fill(0);	            newtondifA0.fill(0);
						  do
						  {

							  errmax1 = 0.0;
							  errmax2 = 0.0;
							  errmax = 0.0;
							  numit = numit + 1;

							  //-------------------------- Initiation for Water head calculation --------------------------

							  int n, n1; n = 0; n1 = 0;
							  h[0] = 0.0;

							  for (i = 1; i <= nboun1; i++)
							  {//81
								  n1 = icount[i];
								  for (j = 1; j <= n1; j++)
								  {//73
									  n = n + 1;
									  h[n] = hb1[i];
									  newtonf(n) = h[n] - hb1[i];
									  newtondif(n, n) = 1;
									  newtonh0(n) = h[n];
								  }//73 
							  }//81 



							  //---------------the following deals with the second boundary----------------------------------------
							  ib2 = 0;
							  inode = npb21 - 1;
							  inode1 = inode;

							  for (nb2 = nboun1 + 1; nb2 <= nboun1 + nboun2; nb2++)
							  {//620
								  ib2 = ib2 + 1;
								  allkb = 0.0;
								  for (j = 1; j <= icount[nb2]; j++)
								  {//605
									  inode1 = inode1 + 1;
									  inode2 = int(nijkl[inode1][maj]);
									  b3 = nijkl[inode1][3 * maj];

									  kbl = b3 * b3 * b3;

									  allkb = allkb + kbl;
								  }//605	continue

								  for (k = 1; k <= icount[nb2]; k++)
								  {//610
									  inode = inode + 1;
									  inode2 = int(nijkl[inode][maj]);
									  b3 = nijkl[inode][3 * maj];
									  xl = nijkl[inode][2 * maj];
									  qq = qb2[ib2] * b3 * b3 * b3 / allkb;

									  if (qq == 0)
									  {
										  h[inode] = h[inode2];
										  newtonf(inode) = h[inode] - h[inode2];
										  for (kk = 1; kk <= np1n2; kk++)
										  {
											  if (kk == inode)
											  {
												  newtondif(inode, kk) = 1;
											  }
											  else if (kk == inode2)
											  {
												  newtondif(inode, kk) = -1;
											  }
											  else
											  {
												  newtondif(inode, kk) = 0;
											  }

											  newtonh0(inode) = h[inode];
										  }
									  }
									  else
									  {
										  Vl = fabs(qq / b3);
										  reynold = b3 * Vl * 2.0 / 0.01;
										  if (reynold < 2300)
										  {
											  J = qq * 12 * 0.01 / b3 / b3 / b3 / 980.0;
											  h[inode] = h[inode2] + J * xl;
										  }
										  else {
											  qq1 = fabs(qq) / 4.70 / b3;
											  J = pow(qq1, 7.0) * 0.01 / b3 / b3 / b3 / b3 / b3 / 980.0 / 980.0 / 980.0 / 980.0;
											  J = pow(J, 0.25) * qq / fabs(qq);
											  h[inode] = h[inode2] + J * xl;
										  }


										  newtonf(inode) = h[inode] - h[inode2] - J * xl;
										  for (kk = 1; kk <= np1n2; kk++)
										  {
											  if (kk == inode)
											  {
												  newtondif(inode, kk) = 1;
											  }
											  else if (kk == inode2)
											  {
												  newtondif(inode, kk) = -1;
											  }
											  else
											  {
												  newtondif(inode, kk) = 0;
											  }

											  newtonh0(inode) = h[inode];
										  }
									  }
								  }//610 continue
							  }//620	continue

						  //--------------------the following deals with the third boundary condition-------------------------------------------------

							  for (i = npbf1; i <= npbf2; i++)
							  {
								  h[i] = 0.0;
								  newtonf(i) = h[i];
								  newtondif(i, i) = 1;
								  newtonh0(i) = h[i];
							  }


							  //--------------------------the following deals with the spring boundary--------------------------------------------------

							  Qspring = 0.0;
							  for (i = np1n1; i <= nodes; i++)
								  for (j = 1; j <= 4; j++)
								  {
									  nnumber = int(nijkl[i][j]);
									  b = nijkl[i][j + 2 * maj];
									  xl = nijkl[i][j + maj];
									  if (fabs(xynp[nnumber][1]) > 0.01 && xynp[nnumber][1] < (maxx - 0.01))continue;

									  if (nnumber > npb22 && nnumber < np1n1)
									  {
										  J = (h[i] - height[nnumber]) / xl;
										  kbl = 980.0 * b * b * b / 12.0 / 0.01;
										  Vl = kbl * J / b;
										  Re[i][nnumber] = 2 * b * Vl / 0.01;
										  if (Re[i][nnumber] > 2300) { Qspring = 464.60533 * pow(b, 1.7143) * pow(fabs(J), 0.5714); }
										  else Qspring = kbl * J;
										  // when Qsrping is greater than 0 , it means it can be regared as a spring
										  if (h[i] < (height[nnumber]))
										  {
											  h[nnumber] = 0;
											  newtonf(nnumber) = h[nnumber];
											  newtondif(nnumber, nnumber) = 1;
											  newtonh0(nnumber) = h[nnumber];
										  }

										  if (h[i] >= (height[nnumber]))
										  {
											  h[nnumber] = height[nnumber];
											  newtonf(nnumber) = h[nnumber] - height[nnumber];
											  newtondif(nnumber, nnumber) = 1;
											  newtonh0(nnumber) = h[nnumber];
										  }
									  }
								  }




							  //----------------------the following deals with the rainfall recharge on water table-------------------------------

							  for (i = 0; i < 500; i++)
								  for (j = 0; j < 4; j++)  qab2[i][j] = 0;


							  number_flow = 0;
							  for (i = np1n1; i <= nodes; i++)
							  {
								  nnumber1 = 0; nnumber2 = 0;
								  if (height[i] <= h[i])
								  {
									  for (j = 1; j <= 4; j++)
									  {
										  nnumber = int(nijkl[i][j]);
										  if (nnumber == 0)continue;
										  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

										  if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
										  {
											  nnumber1 = nnumber1 + 1;
										  }

										  if (h[nnumber] >= height[nnumber])
										  {
											  nnumber2 = nnumber2 + 1;
										  }
									  }

									  if (nnumber1 != 0 && nnumber2 >= 2)
										  for (j = 1; j <= 4; j++)
										  {
											  nnumber = int(nijkl[i][j]);
											  if (nnumber == 0)continue;
											  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

											  if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
											  {
												  number_flow = number_flow + 1;
												  qab2[number_flow][1] = i;
												  qab2[number_flow][2] = nnumber;
											  }

										  }


								  }
							  }



							  //-----------------the following calculates the flow rate of each fracture segment-------------
							  tran_coeff = 10.0 / 365.0 / 24.0 / 60.0 / 60.0;
							  //the unit of Qsurfcorro is mm/a, and the unit of Qsurfcorro1 is ml/s
							  //tran_coeff is the transferring coefficient.
							  Qsurfcorro1 = tran_coeff * (Qsurfcorro * (maxx - 0.0) * 1.0);//*(maxx-minx)

							  sumwidth = 0;
							  ninfil = 0;
							  for (i = 1; i <= number_flow; i++)
							  {
								  inode = int(qab2[i][1]);
								  for (k = 1; k <= 4; k++)
								  {
									  nnumber = int(nijkl[inode][k]);

									  if (nnumber == int(qab2[i][2]))
									  {
										  sumwidth = sumwidth + nijkl[inode][k + 2 * maj];
										  ninfil++;
									  }
								  }
							  }
							  for (i = 1; i <= number_flow; i++)
							  {
								  inode = int(qab2[i][1]);
								  for (k = 1; k <= 4; k++)
								  {
									  nnumber = int(nijkl[inode][k]);
									  if (nnumber == int(qab2[i][2]))//&& h[inode]<maxy
									  {

										  qab2[i][3] = Qsurfcorro1 / ninfil;//evenly distribute infiltration																	
									  }

								  }

							  }





							  //-------------------------------------------------inner node matrix------------------------------------
							  //------------------------------------------------------------------------------------------------------
							  //the following begins to calculate the water head of the inner nodes


							  c = 98000.0 / 12.0;   // c is pg/12u and the transferring of the unit
							  c2 = 464.60533;

							  for (i = np1n1; i <= nodes; i++)
							  { // 255内节点i                 

							//-----------------------------------------------inner dry nodes and their matrix-------------------------------------------

								  if (height[i] > h[i])
								  {
									  h[i] = 0.0;//eliminate the dry point

									  newtonf(i) = h[i];//

									  newtondif(i, i) = 1;
									  newtonh0(i) = h[i];
								  }

								  //-----------------------------------------inner wet nodes and their matrix-------------------------------

								  if (height[i] <= h[i])
								  { // available//1951



											  // situ1: the inner nodes with some unsatured a neigboring nodes
									  available1 = 0; available2 = 0;
									  sumbl = 0.0;
									  sumhbl = 0.0;
									  sumqbl = 0.0;
									  //errmax1=0.0 ;
									  err1 = 0.0;
									  NRe = 0;
									  m = 0;
									  z = 0;
									  gan = 0;
									  gan2 = 0;
									  nnumber1 = 0;
									  nnumber2 = 0;
									  nnumber3 = 0;
									  nnumber4 = 0;
									  b1 = 0;
									  b2 = 0;
									  b3 = 0;
									  b4 = 0;
									  xl1 = 0;
									  xl2 = 0;
									  xl3 = 0;
									  xl4 = 0;
									  H1 = 0;
									  H2 = 0;
									  H3 = 0;
									  H4 = 0;
									  qj1 = 0;
									  qj2 = 0;
									  qj3 = 0;
									  qj4 = 0;
									  fenx1 = 0;
									  fenx2 = 0;
									  fenx0 = 0;
									  aa1 = 0;
									  bb1 = 0;
									  cc1 = 0;
									  for (j = 0; j < 5; j++)
									  {
										  wshi[j] = 0;
										  cshi[j] = 0;
									  }

									  for (j = 1; j <= 4; j++)
									  {
										  nnumber = int(nijkl[i][j]);
										  if (nnumber == 0)continue;
										  if (h[nnumber] < height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
										  {
											  available1 = available1 + 1;//available1 = 1;
										  }   //h[nnumber]=0;改这里也可以，从降雨补给的分配节点上判断available有几个
										  if (h[nnumber] >= height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
										  {
											  available2 = available2 + 1;//
										  }

									  }

									  //--------------------------------------已经判断出来有一个、两个还是0个干邻点-----------------------------------
									  //--------------------------------------------------------------------------------------------------------------




									  //-----------------------------------------1 wet neighbour nodes------------------------------------------
									  //--------------------------------------------------------------------------------------------------------
									  if (available2 == 1)
									  {
										  newtonh0(i) = h[i];
										  for (j = 1; j <= 4; j++)
										  {
											  nnumber = int(nijkl[i][j]);
											  if (nnumber == 0)continue;
											  if (h[nnumber] >= height[nnumber])
											  {
												  newtonf(i) = h[i] - h[nnumber];
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber)
													  {
														  newtondif(i, kk) = -1;
													  }

													  else if (kk == i)
													  {
														  newtondif(i, kk) = 1;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }
												  }
											  }
										  }
									  }
									  //------------------------------------2 or more wet neighbour nodes---------------------------------------------
									  if (available2 >= 2)
									  {


										  //-----------------------------------------有1个高位的干邻点--------------------------------------------
										  if (available1 == 1)//中心节点是湿的，邻点是干的，有一个干节点就执行这个
										  { // available1
											  for (j = 1; j <= 4; j++)
											  { //186
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] <= h[nnumber])//这个是邻节点周围的湿节点,只有找到湿节点才判断层流紊流
												  {
													  b = nijkl[i][j + 2 * maj];
													  xl = nijkl[i][j + maj];
													  J = (h[i] - h[nnumber]) / xl;
													  kbl = 980.0 * b * b * b / 12.0 / 0.01;
													  Vl = kbl * J / b;
													  Re[i][nnumber] = 2 * b * Vl / 0.01;
													  if (fabs(Re[i][nnumber]) > 2300)
													  {
														  NRe = NRe + 1; z = z + 1; wshi[z] = j;
													  }
													  else
													  {
														  m = m + 1; cshi[m] = j;
													  }
													  //sumbl=sumbl+b*b*b/xl;//sumhbl=sumhbl+h[nnumber]*b*b*b/xl;		
												  }
											  } //! 186


											  for (j = 1; j <= 4; j++)
											  { //187
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] > h[nnumber])//如果是干节点
												  {
													  for (k = 1; k <= number_flow; k++)
														  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
														  {
															  //sumqbl=sumqbl+qab2[k][3]/c/sumbl;
															  //Vl=qab2[k][3]/b;
															  //Re[i][nnumber]=2*b*Vl/0.01;
															  gan = k;
														  }

												  }
											  } // 187


											  if (NRe < 1)
											  {
												  //h[i]=sumhbl/sumbl ;
												  fenn = 0;
												  qj1 = int(cshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i]; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);



												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];//这点比较特殊，只要把流量在i节点的水均衡计算中表示出来就可以了


												  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }

											  if (NRe == 1)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[1]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[2]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }
											  if (NRe == 2)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[1]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }

											  if (NRe == 3)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(wshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

												  qj4 = int(gan);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }





										  } // if available1>0


								  //-----------------------------------------有2个高位的干邻点--------------------------------------------
								  //---------------------------------有两个干节点的降雨情况-------------------------------------
								  //--------------------------两个干节点的高度如果都低于中心节点高度，不接受降雨补给------------
										  if (available1 == 2)
										  { // available1
											  for (j = 1; j <= 4; j++)
											  { //186
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] <= h[nnumber])
												  {
													  b = nijkl[i][j + 2 * maj];
													  xl = nijkl[i][j + maj];
													  J = (h[i] - h[nnumber]) / xl;
													  kbl = 980.0 * b * b * b / 12.0 / 0.01;
													  Vl = kbl * J / b;
													  Re[i][nnumber] = 2 * b * Vl / 0.01;
													  if (fabs(Re[i][nnumber]) > 2300)
													  {
														  NRe = NRe + 1; z = z + 1; wshi[z] = j;
													  }
													  else
													  {
														  m = m + 1; cshi[m] = j;
													  }

												  }
											  } // !186

											  for (j = 1; j <= 4; j++)
											  { //187
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (height[nnumber] > h[nnumber])//干节点
												  {
													  for (k = 1; k <= number_flow; k++)
														  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
														  {

															  gan = k;
														  }
												  }
											  } //! 187

											  for (j = 1; j <= 4; j++)
											  { //188
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;//
												  if (height[nnumber] > h[nnumber])//是干节点
												  {
													  for (k = 1; k <= number_flow; k++)
														  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
														  {
															  if (k != gan)
																  gan2 = k;
														  }
												  }
											  } //! 188



											  if (NRe < 1)
											  {
												  //h[i]=sumhbl/sumbl ;
												  fenn = 0;
												  qj1 = int(cshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i]; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(gan);
												  nnumber3 = int(qab2[qj3][2]);
												  q3 = qab2[qj3][3];

												  qj4 = int(gan2);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][3];

												  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }

											  if (NRe == 1)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[1]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(gan);
												  nnumber3 = int(qab2[qj3][2]);
												  q3 = qab2[qj3][3];

												  qj4 = int(gan2);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][2];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }
											  if (NRe == 2)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(gan);
												  nnumber3 = int(qab2[qj3][2]);
												  q3 = qab2[qj3][3];

												  qj4 = int(gan2);
												  nnumber4 = int(qab2[qj4][2]);
												  q4 = qab2[qj4][2];

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) - q3 - q4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = 0;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + 0 + 0;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }


										  } // if available1 = 2,有两个干节点，完成。





								  //------------------------------没有干节点，湿节点周围全都是湿节点-------------------

										  if (available1 == 0)
										  { // available1
											  for (j = 1; j <= maj; j++)
											  {//15
												  inode1 = i;
												  inode2 = int(nijkl[inode1][j]);
												  xl = nijkl[inode1][j + maj];
												  b = nijkl[inode1][j + 2 * maj];
												  if (b == 0.0) continue;
												  if (height[inode2] <= h[inode2])//available为0代表的是没有降雨，但仍有可能有低于中心点的干节点，
													  //如果不判断，这类干节点（所有干节点都低于中心点没有降雨补给的情况）被当做层流，会以水头0参与运算
												  {
													  J = (h[inode1] - h[inode2]) / xl;
													  kbl = 980 * b * b * b / 12 / 0.01;
													  Vl = kbl * J / b;
													  Re[inode1][inode2] = 2 * b * Vl / 0.01;

													  //m = m + 1; cshi[m] = j;
													  if (fabs(Re[inode1][inode2]) > 2300)
													  {
														  NRe = NRe + 1; z = z + 1; wshi[z] = j;
													  }
													  else { m = m + 1; cshi[m] = j; }

												  }//

											  }//15 continue
											  if (NRe < 1)
											  {
												  //h[i]=sumhbl/sumbl ;
												  fenn = 0;
												  qj1 = int(cshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);//np0 connected with np1，np2，np3，np4
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i]; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[4]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }

											  }

											  if (NRe == 1)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(cshi[1]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i]; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[2]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[3]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];
												  }
											  }

											  if (NRe == 2)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(cshi[1]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i]; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[2]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }

											  if (NRe == 3)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(wshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

												  qj4 = int(cshi[1]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i]; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c * b4 * b4 * b4 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }

											  if (NRe == 4)
											  {
												  fenn = 0;
												  qj1 = int(wshi[1]);
												  if (qj1 > 0)
												  {
													  nnumber1 = int(nijkl[i][qj1]);
													  b1 = nijkl[i][qj1 + 2 * maj];
													  xl1 = nijkl[i][qj1 + maj];
													  H1 = h[nnumber1];
												  }
												  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

												  qj2 = int(wshi[2]);
												  if (qj2 > 0)
												  {
													  nnumber2 = int(nijkl[i][qj2]);
													  b2 = nijkl[i][qj2 + 2 * maj];
													  xl2 = nijkl[i][qj2 + maj];
													  H2 = h[nnumber2];
												  }
												  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

												  qj3 = int(wshi[3]);
												  if (qj3 > 0)
												  {
													  nnumber3 = int(nijkl[i][qj3]);
													  b3 = nijkl[i][qj3 + 2 * maj];
													  xl3 = nijkl[i][qj3 + maj];
													  H3 = h[nnumber3];
												  }
												  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

												  qj4 = int(wshi[4]);
												  if (qj4 > 0)
												  {
													  nnumber4 = int(nijkl[i][qj4]);
													  b4 = nijkl[i][qj4 + 2 * maj];
													  xl4 = nijkl[i][qj4 + maj];
													  H4 = h[nnumber4];
												  }
												  else { H4 = h[i] + 0.00001; b4 = 0; xl4 = 1; }

												  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c2 * pow(b4, 1.7143) * (h[i] - H4) / fabs(h[i] - H4) * pow(fabs((h[i] - H4) / xl4), 0.5714);
												  for (kk = 1; kk <= np1n2; kk++)
												  {
													  if (kk == nnumber1)
													  {
														  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
													  }
													  else if (kk == nnumber2)
													  {
														  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
													  }
													  else if (kk == nnumber3)
													  {
														  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
													  }
													  else if (kk == nnumber4)
													  {
														  newtondif(i, kk) = -c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
													  }
													  else if (kk == i)
													  {
														  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
													  }
													  else
													  {
														  newtondif(i, kk) = 0;
													  }

													  newtonh0(i) = h[i];

												  }

											  }



										  } // ! if available1=0

									  }

								  } // i node, water budget equation and matrix
									//1951



							  }  // 255inner node i，head vector，water budget vector，and matrix。


					   //-------------------------------------------------inner node matrix completed----------------------------
					   //--------------------------------------------------------------------------------------------------------



					   //------------------------------matrix computing---------------------------------------------

							  int iii;     iii = np1n2;
							  mat newtondifA(iii, iii), newtondifB(iii, iii), newtondifB0(iii, iii);

							  vec newtonh0A(iii);    vec newtonfA(iii);
							  vec newtonh0B(iii);    vec newtontol(iii);


							  for (i = 0; i <= (np1n2 - 1); i++)
							  {
								  int jjj;  jjj = i + 1;
								  newtonfA(i) = newtonf(jjj);
								  newtonh0A(i) = newtonh0(jjj);
								  for (j = 0; j <= (np1n2 - 1); j++)
								  {
									  int kkk;      kkk = j + 1;
									  newtondifA(i, j) = newtondif(jjj, kkk);
								  }
							  }
							  newtondifB0 = 1.0 * newtondifA;




							  {

								  try
								  {
									  newtondifB = inv(newtondifB0);
								  }

								  catch (...)
								  {
									  cout << " error01catch ";
									  numit = 30;
									  newtondifB = newtondifB0;


								  }
								  newtontol = newtondifB * 1.0 * newtonfA;

							  }


							  //newtontol = solve(newtondifA,newtonfA);
							  newtonh0B = newtonh0A - newtontol;

							  for (i = 0; i <= (np1n2 - 1); i++)
							  {
								  if (fabs(newtontol(i)) > errmax) 	errmax = fabs(newtontol(i));
							  }

							  for (i = npb21 - 1; i <= (np1n2 - 1); i++)
							  {
								  h[i + 1] = newtonh0B(i);
							  }









							  //-----------------------------------------  water budget and error  ----------------------------


							  for (inode = 1; inode <= nodes; inode++)
							  {
								  qbalance[inode] = 0;
								  for (j = 1; j <= maj; j++)
								  {
									  inode1 = int(nijkl[inode][j]);
									  Qcell[inode][j] = 0.0; Re[inode][inode1] = 0;
									  if (fabs(h[inode]) < height[inode])continue;

									  if (inode1 < 1)continue;
									  if (inode < np1n1 && inode1 < np1n1)continue;

									  if (h[inode1] >= height[inode1])//h[inode1]!= 0
									  {
										  b = nijkl[inode][2 * maj + j];
										  xl = nijkl[inode][maj + j];
										  J = (h[inode] - h[inode1]) / xl;
										  kbl = 980.0 * b * b * b / 12.0 / 0.01;
										  Vl = kbl * J / b;
										  Re[inode][inode1] = 2 * b * Vl / 0.01;
										  if (fabs(Re[inode][inode1]) < 2300) { Qcell[inode][j] = kbl * J; }
										  else Qcell[inode][j] = 464.60533 * (h[inode] - h[inode1]) / fabs(h[inode] - h[inode1]) * pow(b, 1.7143) * pow(fabs(J), 0.5714);
									  }
									  else
									  {
										  for (k = 1; k <= number_flow; k++)
											  if (inode == int(qab2[k][1]) && inode1 == int(qab2[k][2]))
												  Qcell[inode][j] = -qab2[k][3];

									  }

									  qbalance[inode] = qbalance[inode] + Qcell[inode][j];
								  }
							  }
							  qqbudget = 0;

							  for (inode = np1n1; inode <= nodes; inode++)
							  {
								  qqbudget = qqbudget + qbalance[inode];
							  }



						  } while ((errmax > tolmax && numit < 30) || (fabs(qqbudget) > 0.0001 && numit < 30));

						  cout << " numit02   " << numit << " errmax02  " << errmax << "  qqbudget02  " << qqbudget ;

						  //-------------------------------------------------------------  inner iteration  ----------------------------------------------------
						  //------------------------------------------------------------   inner core   --------------------------------------------------------


						  double qqboundout2 = 0; double qqboundin2 = 0;
						  for (i = 1; i <= npbf2; i++)
						  {
							  qqboundout2 = qqboundout2 + qbalance[i];//Qcell[i][maj]node
						  }
						  for (i = 1; i <= number_flow; i++)
						  {
							  qqboundin2 = qqboundin2 + qab2[i][3];

						  }
						  cout << " qqboundout2 " << qqboundout2 << " qqboundin2 " << qqboundin2 << " qoutbudget2  " << (qqboundout2 + qqboundin2) << endl;

						  



						  if (numit == 30)
						  {
							  for (i = 1; i <= layerwetnodei; i++)
							  {
								  inode = int(layerwetnode[i][0]);
								  h[inode] = layerwetnode[i][1];
							  }

							  for (i = 1; i <= layerdrynodei; i++)
							  {
								  inode = layerdrynode[i];
								  h[inode] = 0;
							  }






							  //-----------------------------------  inner iteration of confined water table,   third time -----------------------------------------
							  //------------------------------------------------------------   inner core   --------------------------------------------------------

							  numit = 0;	      newtonh0.fill(0);      newtonf.fill(0);       newtondif.fill(0);	            newtondifA0.fill(0);
							  do
							  {

								  errmax1 = 0.0;
								  errmax2 = 0.0;
								  errmax = 0.0;
								  numit = numit + 1;

								  //-------------------------- Initiation for Water head calculation --------------------------

								  int n, n1; n = 0; n1 = 0;
								  h[0] = 0.0;

								  for (i = 1; i <= nboun1; i++)
								  {//81
									  n1 = icount[i];
									  for (j = 1; j <= n1; j++)
									  {//73
										  n = n + 1;
										  h[n] = hb1[i];
										  newtonf(n) = h[n] - hb1[i];
										  newtondif(n, n) = 1;
										  newtonh0(n) = h[n];
									  }//73 
								  }//81 



								  //---------------the following deals with the second boundary----------------------------------------
								  ib2 = 0;
								  inode = npb21 - 1;
								  inode1 = inode;

								  for (nb2 = nboun1 + 1; nb2 <= nboun1 + nboun2; nb2++)
								  {//620
									  ib2 = ib2 + 1;
									  allkb = 0.0;
									  for (j = 1; j <= icount[nb2]; j++)
									  {//605
										  inode1 = inode1 + 1;
										  inode2 = int(nijkl[inode1][maj]);
										  b3 = nijkl[inode1][3 * maj];

										  kbl = b3 * b3 * b3;

										  allkb = allkb + kbl;
									  }//605	continue

									  for (k = 1; k <= icount[nb2]; k++)
									  {//610
										  inode = inode + 1;
										  inode2 = int(nijkl[inode][maj]);
										  b3 = nijkl[inode][3 * maj];
										  xl = nijkl[inode][2 * maj];
										  qq = qb2[ib2] * b3 * b3 * b3 / allkb;

										  if (qq == 0)
										  {
											  h[inode] = h[inode2];
											  newtonf(inode) = h[inode] - h[inode2];
											  for (kk = 1; kk <= np1n2; kk++)
											  {
												  if (kk == inode)
												  {
													  newtondif(inode, kk) = 1;
												  }
												  else if (kk == inode2)
												  {
													  newtondif(inode, kk) = -1;
												  }
												  else
												  {
													  newtondif(inode, kk) = 0;
												  }

												  newtonh0(inode) = h[inode];
											  }
										  }
										  else
										  {
											  Vl = fabs(qq / b3);
											  reynold = b3 * Vl * 2.0 / 0.01;
											  if (reynold < 2300)
											  {
												  J = qq * 12 * 0.01 / b3 / b3 / b3 / 980.0;
												  h[inode] = h[inode2] + J * xl;
											  }
											  else {
												  qq1 = fabs(qq) / 4.70 / b3;
												  J = pow(qq1, 7.0) * 0.01 / b3 / b3 / b3 / b3 / b3 / 980.0 / 980.0 / 980.0 / 980.0;
												  J = pow(J, 0.25) * qq / fabs(qq);
												  h[inode] = h[inode2] + J * xl;
											  }


											  newtonf(inode) = h[inode] - h[inode2] - J * xl;
											  for (kk = 1; kk <= np1n2; kk++)
											  {
												  if (kk == inode)
												  {
													  newtondif(inode, kk) = 1;
												  }
												  else if (kk == inode2)
												  {
													  newtondif(inode, kk) = -1;
												  }
												  else
												  {
													  newtondif(inode, kk) = 0;
												  }

												  newtonh0(inode) = h[inode];
											  }
										  }
									  }//610 continue
								  }//620	continue

							  //--------------------the following deals with the third boundary condition-------------------------------------------------

								  for (i = npbf1; i <= npbf2; i++)
								  {
									  h[i] = 0.0;
									  newtonf(i) = h[i];
									  newtondif(i, i) = 1;
									  newtonh0(i) = h[i];
								  }


								  //--------------------------the following deals with the spring boundary--------------------------------------------------

								  Qspring = 0.0;
								  for (i = np1n1; i <= nodes; i++)
									  for (j = 1; j <= 4; j++)
									  {
										  nnumber = int(nijkl[i][j]);
										  b = nijkl[i][j + 2 * maj];
										  xl = nijkl[i][j + maj];
										  if (fabs(xynp[nnumber][1]) > 0.01 && xynp[nnumber][1] < (maxx - 0.01))continue;

										  if (nnumber > npb22 && nnumber < np1n1)
										  {
											  J = (h[i] - height[nnumber]) / xl;
											  kbl = 980.0 * b * b * b / 12.0 / 0.01;
											  Vl = kbl * J / b;
											  Re[i][nnumber] = 2 * b * Vl / 0.01;
											  if (Re[i][nnumber] > 2300) { Qspring = 464.60533 * pow(b, 1.7143) * pow(fabs(J), 0.5714); }
											  else Qspring = kbl * J;
											  // when Qsrping is greater than 0 , it means it can be regared as a spring
											  if (h[i] < (height[nnumber]))
											  {
												  h[nnumber] = 0;
												  newtonf(nnumber) = h[nnumber];
												  newtondif(nnumber, nnumber) = 1;
												  newtonh0(nnumber) = h[nnumber];
											  }

											  if (h[i] >= (height[nnumber]))
											  {
												  h[nnumber] = height[nnumber];
												  newtonf(nnumber) = h[nnumber] - height[nnumber];
												  newtondif(nnumber, nnumber) = 1;
												  newtonh0(nnumber) = h[nnumber];
											  }
										  }
									  }




								  //----------------------the following deals with the rainfall recharge on water table-------------------------------

								  for (i = 0; i < 500; i++)
									  for (j = 0; j < 4; j++)  qab2[i][j] = 0;


								  number_flow = 0;
								  for (i = np1n1; i <= nodes; i++)
								  {
									  nnumber1 = 0; nnumber2 = 0;
									  if (height[i] <= h[i])
									  {
										  for (j = 1; j <= 4; j++)
										  {
											  nnumber = int(nijkl[i][j]);
											  if (nnumber == 0)continue;
											  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

											  if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
											  {
												  nnumber1 = nnumber1 + 1;
											  }

											  if (h[nnumber] >= height[nnumber])
											  {
												  nnumber2 = nnumber2 + 1;
											  }
										  }

										  if (nnumber1 != 0 && nnumber2 >= 2)
											  for (j = 1; j <= 4; j++)
											  {
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (fabs(xynp[nnumber][1]) < (0.05) || xynp[nnumber][1] > (maxx - 0.05))continue;

												  if (height[i] < height[nnumber] && h[nnumber] < height[nnumber])
												  {
													  number_flow = number_flow + 1;
													  qab2[number_flow][1] = i;
													  qab2[number_flow][2] = nnumber;
												  }

											  }


									  }
								  }



								  //-----------------the following calculates the flow rate of each fracture segment-------------
								  tran_coeff = 10.0 / 365.0 / 24.0 / 60.0 / 60.0;
								  //the unit of Qsurfcorro is mm/a, and the unit of Qsurfcorro1 is ml/s
								  //tran_coeff is the transferring coefficient.
								  Qsurfcorro1 = tran_coeff * (Qsurfcorro * (maxx - 0.0) * 1.0);//*(maxx-minx)

								  sumwidth = 0;
								  ninfil = 0;
								  for (i = 1; i <= number_flow; i++)
								  {
									  inode = int(qab2[i][1]);
									  for (k = 1; k <= 4; k++)
									  {
										  nnumber = int(nijkl[inode][k]);

										  if (nnumber == int(qab2[i][2]))
										  {
											  sumwidth = sumwidth + nijkl[inode][k + 2 * maj];
											  ninfil++;
										  }
									  }
								  }
								  for (i = 1; i <= number_flow; i++)
								  {
									  inode = int(qab2[i][1]);
									  for (k = 1; k <= 4; k++)
									  {
										  nnumber = int(nijkl[inode][k]);
										  if (nnumber == int(qab2[i][2]))//&& h[inode]<maxy
										  {

											  qab2[i][3] = Qsurfcorro1 / ninfil;//evenly distribute infiltration																	
										  }

									  }

								  }





								  //-------------------------------------------------inner node matrix------------------------------------
								  //------------------------------------------------------------------------------------------------------
								  //the following begins to calculate the water head of the inner nodes


								  c = 98000.0 / 12.0;   // c is pg/12u and the transferring of the unit
								  c2 = 464.60533;

								  for (i = np1n1; i <= nodes; i++)
								  { // 255内节点i                 

								//-----------------------------------------------inner dry nodes and their matrix-------------------------------------------

									  if (height[i] > h[i])
									  {
										  h[i] = 0.0;//eliminate the dry point

										  newtonf(i) = h[i];//

										  newtondif(i, i) = 1;
										  newtonh0(i) = h[i];
									  }

									  //-----------------------------------------inner wet nodes and their matrix-------------------------------

									  if (height[i] <= h[i])
									  { // available//1951



												  // situ1: the inner nodes with some unsatured a neigboring nodes
										  available1 = 0; available2 = 0;
										  sumbl = 0.0;
										  sumhbl = 0.0;
										  sumqbl = 0.0;
										  //errmax1=0.0 ;
										  err1 = 0.0;
										  NRe = 0;
										  m = 0;
										  z = 0;
										  gan = 0;
										  gan2 = 0;
										  nnumber1 = 0;
										  nnumber2 = 0;
										  nnumber3 = 0;
										  nnumber4 = 0;
										  b1 = 0;
										  b2 = 0;
										  b3 = 0;
										  b4 = 0;
										  xl1 = 0;
										  xl2 = 0;
										  xl3 = 0;
										  xl4 = 0;
										  H1 = 0;
										  H2 = 0;
										  H3 = 0;
										  H4 = 0;
										  qj1 = 0;
										  qj2 = 0;
										  qj3 = 0;
										  qj4 = 0;
										  fenx1 = 0;
										  fenx2 = 0;
										  fenx0 = 0;
										  aa1 = 0;
										  bb1 = 0;
										  cc1 = 0;
										  for (j = 0; j < 5; j++)
										  {
											  wshi[j] = 0;
											  cshi[j] = 0;
										  }

										  for (j = 1; j <= 4; j++)
										  {
											  nnumber = int(nijkl[i][j]);
											  if (nnumber == 0)continue;
											  if (h[nnumber] < height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
											  {
												  available1 = available1 + 1;//available1 = 1;
											  }   //h[nnumber]=0;改这里也可以，从降雨补给的分配节点上判断available有几个
											  if (h[nnumber] >= height[nnumber])//height[i] < height[nnumber] && h[nnumber] < height[nnumber]存在高位干节点
											  {
												  available2 = available2 + 1;//
											  }

										  }

										  //--------------------------------------已经判断出来有一个、两个还是0个干邻点-----------------------------------
										  //--------------------------------------------------------------------------------------------------------------




										  //-----------------------------------------1 wet neighbour nodes------------------------------------------
										  //--------------------------------------------------------------------------------------------------------
										  if (available2 == 1)
										  {
											  newtonh0(i) = h[i];
											  for (j = 1; j <= 4; j++)
											  {
												  nnumber = int(nijkl[i][j]);
												  if (nnumber == 0)continue;
												  if (h[nnumber] >= height[nnumber])
												  {
													  newtonf(i) = h[i] - h[nnumber];
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber)
														  {
															  newtondif(i, kk) = -1;
														  }

														  else if (kk == i)
														  {
															  newtondif(i, kk) = 1;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }
													  }
												  }
											  }
										  }
										  //------------------------------------2 or more wet neighbour nodes---------------------------------------------
										  if (available2 >= 2)
										  {


											  //-----------------------------------------有1个高位的干邻点--------------------------------------------
											  if (available1 == 1)//中心节点是湿的，邻点是干的，有一个干节点就执行这个
											  { // available1
												  for (j = 1; j <= 4; j++)
												  { //186
													  nnumber = int(nijkl[i][j]);
													  if (nnumber == 0)continue;
													  if (height[nnumber] <= h[nnumber])//这个是邻节点周围的湿节点,只有找到湿节点才判断层流紊流
													  {
														  b = nijkl[i][j + 2 * maj];
														  xl = nijkl[i][j + maj];
														  J = (h[i] - h[nnumber]) / xl;
														  kbl = 980.0 * b * b * b / 12.0 / 0.01;
														  Vl = kbl * J / b;
														  Re[i][nnumber] = 2 * b * Vl / 0.01;
														  if (fabs(Re[i][nnumber]) > 2300)
														  {
															  NRe = NRe + 1; z = z + 1; wshi[z] = j;
														  }
														  else
														  {
															  m = m + 1; cshi[m] = j;
														  }
														  //sumbl=sumbl+b*b*b/xl;//sumhbl=sumhbl+h[nnumber]*b*b*b/xl;		
													  }
												  } //! 186


												  for (j = 1; j <= 4; j++)
												  { //187
													  nnumber = int(nijkl[i][j]);
													  if (nnumber == 0)continue;
													  if (height[nnumber] > h[nnumber])//如果是干节点
													  {
														  for (k = 1; k <= number_flow; k++)
															  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
															  {
																  //sumqbl=sumqbl+qab2[k][3]/c/sumbl;
																  //Vl=qab2[k][3]/b;
																  //Re[i][nnumber]=2*b*Vl/0.01;
																  gan = k;
															  }

													  }
												  } // 187


												  if (NRe < 1)
												  {
													  //h[i]=sumhbl/sumbl ;
													  fenn = 0;
													  qj1 = int(cshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i]; b1 = 0; xl1 = 1; }

													  qj2 = int(cshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i]; b2 = 0; xl2 = 1; }

													  qj3 = int(cshi[3]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i]; b3 = 0; xl3 = 1; }

													  qj4 = int(gan);



													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][3];//这点比较特殊，只要把流量在i节点的水均衡计算中表示出来就可以了


													  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];
													  }

												  }

												  if (NRe == 1)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(cshi[1]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i]; b2 = 0; xl2 = 1; }

													  qj3 = int(cshi[2]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i]; b3 = 0; xl3 = 1; }

													  qj4 = int(gan);
													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][3];

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];
													  }

												  }
												  if (NRe == 2)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(wshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

													  qj3 = int(cshi[1]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i]; b3 = 0; xl3 = 1; }

													  qj4 = int(gan);
													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][3];

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];

													  }

												  }

												  if (NRe == 3)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(wshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

													  qj3 = int(wshi[3]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

													  qj4 = int(gan);
													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][3];

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];

													  }

												  }





											  } // if available1>0


									  //-----------------------------------------有2个高位的干邻点--------------------------------------------
									  //---------------------------------有两个干节点的降雨情况-------------------------------------
									  //--------------------------两个干节点的高度如果都低于中心节点高度，不接受降雨补给------------
											  if (available1 == 2)
											  { // available1
												  for (j = 1; j <= 4; j++)
												  { //186
													  nnumber = int(nijkl[i][j]);
													  if (nnumber == 0)continue;
													  if (height[nnumber] <= h[nnumber])
													  {
														  b = nijkl[i][j + 2 * maj];
														  xl = nijkl[i][j + maj];
														  J = (h[i] - h[nnumber]) / xl;
														  kbl = 980.0 * b * b * b / 12.0 / 0.01;
														  Vl = kbl * J / b;
														  Re[i][nnumber] = 2 * b * Vl / 0.01;
														  if (fabs(Re[i][nnumber]) > 2300)
														  {
															  NRe = NRe + 1; z = z + 1; wshi[z] = j;
														  }
														  else
														  {
															  m = m + 1; cshi[m] = j;
														  }

													  }
												  } // !186

												  for (j = 1; j <= 4; j++)
												  { //187
													  nnumber = int(nijkl[i][j]);
													  if (nnumber == 0)continue;
													  if (height[nnumber] > h[nnumber])//干节点
													  {
														  for (k = 1; k <= number_flow; k++)
															  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
															  {

																  gan = k;
															  }
													  }
												  } //! 187

												  for (j = 1; j <= 4; j++)
												  { //188
													  nnumber = int(nijkl[i][j]);
													  if (nnumber == 0)continue;//
													  if (height[nnumber] > h[nnumber])//是干节点
													  {
														  for (k = 1; k <= number_flow; k++)
															  if (i == int(qab2[k][1]) && nnumber == int(qab2[k][2]))//将i对应到干湿节点上
															  {
																  if (k != gan)
																	  gan2 = k;
															  }
													  }
												  } //! 188



												  if (NRe < 1)
												  {
													  //h[i]=sumhbl/sumbl ;
													  fenn = 0;
													  qj1 = int(cshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i]; b1 = 0; xl1 = 1; }

													  qj2 = int(cshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i]; b2 = 0; xl2 = 1; }

													  qj3 = int(gan);
													  nnumber3 = int(qab2[qj3][2]);
													  q3 = qab2[qj3][3];

													  qj4 = int(gan2);
													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][3];

													  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];
													  }

												  }

												  if (NRe == 1)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(cshi[1]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i]; b2 = 0; xl2 = 1; }

													  qj3 = int(gan);
													  nnumber3 = int(qab2[qj3][2]);
													  q3 = qab2[qj3][3];

													  qj4 = int(gan2);
													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][2];

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 - q3 - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + 0 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];
													  }

												  }
												  if (NRe == 2)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(wshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

													  qj3 = int(gan);
													  nnumber3 = int(qab2[qj3][2]);
													  q3 = qab2[qj3][3];

													  qj4 = int(gan2);
													  nnumber4 = int(qab2[qj4][2]);
													  q4 = qab2[qj4][2];

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) - q3 - q4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = 0;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + 0 + 0;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];

													  }

												  }


											  } // if available1 = 2,有两个干节点，完成。





									  //------------------------------没有干节点，湿节点周围全都是湿节点-------------------

											  if (available1 == 0)
											  { // available1
												  for (j = 1; j <= maj; j++)
												  {//15
													  inode1 = i;
													  inode2 = int(nijkl[inode1][j]);
													  xl = nijkl[inode1][j + maj];
													  b = nijkl[inode1][j + 2 * maj];
													  if (b == 0.0) continue;
													  if (height[inode2] <= h[inode2])//available为0代表的是没有降雨，但仍有可能有低于中心点的干节点，
														  //如果不判断，这类干节点（所有干节点都低于中心点没有降雨补给的情况）被当做层流，会以水头0参与运算
													  {
														  J = (h[inode1] - h[inode2]) / xl;
														  kbl = 980 * b * b * b / 12 / 0.01;
														  Vl = kbl * J / b;
														  Re[inode1][inode2] = 2 * b * Vl / 0.01;

														  //m = m + 1; cshi[m] = j;
														  if (fabs(Re[inode1][inode2]) > 2300)
														  {
															  NRe = NRe + 1; z = z + 1; wshi[z] = j;
														  }
														  else { m = m + 1; cshi[m] = j; }

													  }//

												  }//15 continue
												  if (NRe < 1)
												  {
													  //h[i]=sumhbl/sumbl ;
													  fenn = 0;
													  qj1 = int(cshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);//np0 connected with np1，np2，np3，np4
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i]; b1 = 0; xl1 = 1; }

													  qj2 = int(cshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i]; b2 = 0; xl2 = 1; }

													  qj3 = int(cshi[3]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i]; b3 = 0; xl3 = 1; }

													  qj4 = int(cshi[4]);
													  if (qj4 > 0)
													  {
														  nnumber4 = int(nijkl[i][qj4]);
														  b4 = nijkl[i][qj4 + 2 * maj];
														  xl4 = nijkl[i][qj4 + maj];
														  H4 = h[nnumber4];
													  }
													  else { H4 = h[i]; b4 = 0; xl4 = 1; }

													  newtonf(i) = c * b1 * b1 * b1 * (h[i] - H1) / xl1 + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c * b1 * b1 * b1 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c * b1 * b1 * b1 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];
													  }

												  }

												  if (NRe == 1)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(cshi[1]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i]; b2 = 0; xl2 = 1; }

													  qj3 = int(cshi[2]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i]; b3 = 0; xl3 = 1; }

													  qj4 = int(cshi[3]);
													  if (qj4 > 0)
													  {
														  nnumber4 = int(nijkl[i][qj4]);
														  b4 = nijkl[i][qj4 + 2 * maj];
														  xl4 = nijkl[i][qj4 + maj];
														  H4 = h[nnumber4];
													  }
													  else { H4 = h[i]; b4 = 0; xl4 = 1; }

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c * b2 * b2 * b2 * (h[i] - H2) / xl2 + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c * b2 * b2 * b2 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c * b2 * b2 * b2 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];
													  }
												  }

												  if (NRe == 2)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(wshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

													  qj3 = int(cshi[1]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i]; b3 = 0; xl3 = 1; }

													  qj4 = int(cshi[2]);
													  if (qj4 > 0)
													  {
														  nnumber4 = int(nijkl[i][qj4]);
														  b4 = nijkl[i][qj4 + 2 * maj];
														  xl4 = nijkl[i][qj4 + maj];
														  H4 = h[nnumber4];
													  }
													  else { H4 = h[i]; b4 = 0; xl4 = 1; }

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c * b3 * b3 * b3 * (h[i] - H3) / xl3 + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c * b3 * b3 * b3 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c * b3 * b3 * b3 / xl3 + c * b4 * b4 * b4 / xl4;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];

													  }

												  }

												  if (NRe == 3)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(wshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

													  qj3 = int(wshi[3]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

													  qj4 = int(cshi[1]);
													  if (qj4 > 0)
													  {
														  nnumber4 = int(nijkl[i][qj4]);
														  b4 = nijkl[i][qj4 + 2 * maj];
														  xl4 = nijkl[i][qj4 + maj];
														  H4 = h[nnumber4];
													  }
													  else { H4 = h[i]; b4 = 0; xl4 = 1; }

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c * b4 * b4 * b4 * (h[i] - H4) / xl4;
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = -c * b4 * b4 * b4 / xl4;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c * b4 * b4 * b4 / xl4;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];

													  }

												  }

												  if (NRe == 4)
												  {
													  fenn = 0;
													  qj1 = int(wshi[1]);
													  if (qj1 > 0)
													  {
														  nnumber1 = int(nijkl[i][qj1]);
														  b1 = nijkl[i][qj1 + 2 * maj];
														  xl1 = nijkl[i][qj1 + maj];
														  H1 = h[nnumber1];
													  }
													  else { H1 = h[i] + 0.00001; b1 = 0; xl1 = 1; }

													  qj2 = int(wshi[2]);
													  if (qj2 > 0)
													  {
														  nnumber2 = int(nijkl[i][qj2]);
														  b2 = nijkl[i][qj2 + 2 * maj];
														  xl2 = nijkl[i][qj2 + maj];
														  H2 = h[nnumber2];
													  }
													  else { H2 = h[i] + 0.00001; b2 = 0; xl2 = 1; }

													  qj3 = int(wshi[3]);
													  if (qj3 > 0)
													  {
														  nnumber3 = int(nijkl[i][qj3]);
														  b3 = nijkl[i][qj3 + 2 * maj];
														  xl3 = nijkl[i][qj3 + maj];
														  H3 = h[nnumber3];
													  }
													  else { H3 = h[i] + 0.00001; b3 = 0; xl3 = 1; }

													  qj4 = int(wshi[4]);
													  if (qj4 > 0)
													  {
														  nnumber4 = int(nijkl[i][qj4]);
														  b4 = nijkl[i][qj4 + 2 * maj];
														  xl4 = nijkl[i][qj4 + maj];
														  H4 = h[nnumber4];
													  }
													  else { H4 = h[i] + 0.00001; b4 = 0; xl4 = 1; }

													  newtonf(i) = c2 * pow(b1, 1.7143) * (h[i] - H1) / fabs(h[i] - H1) * pow(fabs((h[i] - H1) / xl1), 0.5714) + c2 * pow(b2, 1.7143) * (h[i] - H2) / fabs(h[i] - H2) * pow(fabs((h[i] - H2) / xl2), 0.5714) + c2 * pow(b3, 1.7143) * (h[i] - H3) / fabs(h[i] - H3) * pow(fabs((h[i] - H3) / xl3), 0.5714) + c2 * pow(b4, 1.7143) * (h[i] - H4) / fabs(h[i] - H4) * pow(fabs((h[i] - H4) / xl4), 0.5714);
													  for (kk = 1; kk <= np1n2; kk++)
													  {
														  if (kk == nnumber1)
														  {
															  newtondif(i, kk) = -c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1;
														  }
														  else if (kk == nnumber2)
														  {
															  newtondif(i, kk) = -c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2;
														  }
														  else if (kk == nnumber3)
														  {
															  newtondif(i, kk) = -c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3;
														  }
														  else if (kk == nnumber4)
														  {
															  newtondif(i, kk) = -c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
														  }
														  else if (kk == i)
														  {
															  newtondif(i, kk) = c2 * pow(b1, 1.7143) * pow(fabs((h[i] - H1) / xl1), -0.4286) * 0.5714 / xl1 + c2 * pow(b2, 1.7143) * pow(fabs((h[i] - H2) / xl2), -0.4286) * 0.5714 / xl2 + c2 * pow(b3, 1.7143) * pow(fabs((h[i] - H3) / xl3), -0.4286) * 0.5714 / xl3 + c2 * pow(b4, 1.7143) * pow(fabs((h[i] - H4) / xl4), -0.4286) * 0.5714 / xl4;
														  }
														  else
														  {
															  newtondif(i, kk) = 0;
														  }

														  newtonh0(i) = h[i];

													  }

												  }



											  } // ! if available1=0

										  }

									  } // i node, water budget equation and matrix
										//1951



								  }  // 255inner node i，head vector，water budget vector，and matrix。


						   //-------------------------------------------------inner node matrix completed----------------------------
						   //--------------------------------------------------------------------------------------------------------



						   //------------------------------matrix computing---------------------------------------------

								  int iii;     iii = np1n2;
								  mat newtondifA(iii, iii), newtondifB(iii, iii), newtondifB0(iii, iii);

								  vec newtonh0A(iii);    vec newtonfA(iii);
								  vec newtonh0B(iii);    vec newtontol(iii);


								  for (i = 0; i <= (np1n2 - 1); i++)
								  {
									  int jjj;  jjj = i + 1;
									  newtonfA(i) = newtonf(jjj);
									  newtonh0A(i) = newtonh0(jjj);
									  for (j = 0; j <= (np1n2 - 1); j++)
									  {
										  int kkk;      kkk = j + 1;
										  newtondifA(i, j) = newtondif(jjj, kkk);
									  }
								  }
								  newtondifB0 = 1.0 * newtondifA;




								  {

									  try
									  {
										  newtondifB = inv(newtondifB0);
									  }

									  catch (...)
									  {
										  cout << " error01catch ";
										  numit = 30;
										  newtondifB = newtondifB0;


									  }
									  newtontol = newtondifB * 1.0 * newtonfA;

								  }


								  //newtontol = solve(newtondifA,newtonfA);
								  newtonh0B = newtonh0A - newtontol;

								  for (i = 0; i <= (np1n2 - 1); i++)
								  {
									  if (fabs(newtontol(i)) > errmax) 	errmax = fabs(newtontol(i));
								  }

								  for (i = npb21 - 1; i <= (np1n2 - 1); i++)
								  {
									  h[i + 1] = newtonh0B(i);
								  }









								  //-----------------------------------------  water budget and error  ----------------------------


								  for (inode = 1; inode <= nodes; inode++)
								  {
									  qbalance[inode] = 0;
									  for (j = 1; j <= maj; j++)
									  {
										  inode1 = int(nijkl[inode][j]);
										  Qcell[inode][j] = 0.0; Re[inode][inode1] = 0;
										  if (fabs(h[inode]) < height[inode])continue;

										  if (inode1 < 1)continue;
										  if (inode < np1n1 && inode1 < np1n1)continue;

										  if (h[inode1] >= height[inode1])//h[inode1]!= 0
										  {
											  b = nijkl[inode][2 * maj + j];
											  xl = nijkl[inode][maj + j];
											  J = (h[inode] - h[inode1]) / xl;
											  kbl = 980.0 * b * b * b / 12.0 / 0.01;
											  Vl = kbl * J / b;
											  Re[inode][inode1] = 2 * b * Vl / 0.01;
											  if (fabs(Re[inode][inode1]) < 2300) { Qcell[inode][j] = kbl * J; }
											  else Qcell[inode][j] = 464.60533 * (h[inode] - h[inode1]) / fabs(h[inode] - h[inode1]) * pow(b, 1.7143) * pow(fabs(J), 0.5714);
										  }
										  else
										  {
											  for (k = 1; k <= number_flow; k++)
												  if (inode == int(qab2[k][1]) && inode1 == int(qab2[k][2]))
													  Qcell[inode][j] = -qab2[k][3];

										  }

										  qbalance[inode] = qbalance[inode] + Qcell[inode][j];
									  }
								  }
								  qqbudget = 0;

								  for (inode = np1n1; inode <= nodes; inode++)
								  {
									  qqbudget = qqbudget + qbalance[inode];
								  }



							  } while ((errmax > tolmax && numit < 30) || (fabs(qqbudget) > 0.0001 && numit < 30));


							  cout << " numit03   " << numit << " errmax03  " << errmax << "  qqbudget03  " << qqbudget;

							  //-------------------------------------------------------------  inner iteration  ----------------------------------------------------
							  //------------------------------------------------------------   inner core   --------------------------------------------------------


							  double qqboundout3 = 0; double qqboundin3 = 0;
							  for (i = 1; i <= npbf2; i++)
							  {
								  qqboundout3 = qqboundout3 + qbalance[i];//Qcell[i][maj]node
							  }
							  for (i = 1; i <= number_flow; i++)
							  {
								  qqboundin3 = qqboundin3 + qab2[i][3];

							  }
							  cout << " qqboundout3 " << qqboundout3 << " qqboundin3 " << qqboundin3 << " qoutbudget3  "<< (qqboundout3 + qqboundin3) << endl;


							  if (numit == 30)
							  {
								  
								  for (i = 1; i <= layerwetnodei; i++)
								  {
									  inode = int(layerwetnode[i][0]);
									  h[inode] = layerwetnode[i][1];
								  }

								  for (i = 1; i <= layerdrynodei; i++)
								  {
									  inode = layerdrynode[i];
									  h[inode] = 0;
								  }
							  }
						  }


					  }

					  int layerwet22dry = 0;
					  for (i = 1; i <= layerwetnodei; i++)
					  {
						  inode = int(layerwetnode[i][0]);
						  if (height[inode] > h[inode])
						  {

							  layerwet22dry = layerwet22dry + 1;
						  }
					  }
					  int layerwet2max = 0;
					  for (i = 1; i <= nodes; i++)
					  {
						  
						  if (h[inode] > 2.0 * maxy)
						  {

							  layerwet2max = layerwet2max + 1;
						  }
					  }
					  if (layerwet22dry >= (layerwetnodei - 10) || layerwet2max > 0)
					  {

						  for (i = 1; i <= layerwetnodei; i++)
						  {
							  inode = int(layerwetnode[i][0]);
							  h[inode] = layerwetnode[i][1];
						  }

						  for (i = 1; i <= layerdrynodei; i++)
						  {
							  inode = layerdrynode[i];
							  h[inode] = 0;
						  }
					  }


					  if (h[layernodeA] > height[layernodeA] && (layernodeA > 0)) wettest01 = wettest01 + 1;




                    }


				  }



				  //-------------------------------------------------------- midlle layer iteration----------------------------------------
				  //------------------------------------------------------------   !   ----------------------------------------------------



			  }

		




		     {	wettest02 = 0;

				
				for (inode = np1n1; inode <= nodes; inode++)// npbf1
				{

					int wetnode = 0;
					for (j = 1; j <= maj; j++)/////////&& (fabs(h[inode1]) < 0.001)    (height[inode1]-2.0)
					{
						inode1 = int(nijkl[inode][j]);// 												
						if ((inode1 >= 1 ) && (height[inode] > h[inode]) && h[inode1] >= (height[inode] + 0))
							wetnode = wetnode + 1;
					}

					//(height[inode] < h[inode]) &&(inode2 >= np1n1) &&(height[inode] < height[inode2]) &&(fabs(h[inode2]) < 1.0) &&
					if (wetnode >= 1) wettest02 = wettest02 + 1;



					
                } 


             }
				numitouter = numitouter + 1;



        cout << "  numitouter  " << numitouter << "   wettest01   " <<wettest01<<"   wettest02  " << wettest02<< endl ;

} while ((wettest01 >= 1 && numitouter <100)  || (wettest02 >= 2 && numitouter <2) );//




	//---------------------------end------ 3 layers iterations of searching for the free surface-----------------------------
	//------------------------------------------------------------   !   ----------------------------------------------------
	//------------------------------------------------------------   !   ----------------------------------------------------




//-----------------------------------water budget and error--------------------------------
					qq1 =0.0;
					qq2=0.0;
					qoutmax=0.0;
					for(inode=1;inode<=npbf2;inode++)
					{//2000
					q[inode] = Qcell[inode][maj];

					if(qoutmax<fabs(q[inode]))qoutmax=fabs(q[inode]);
					
					qq1 = qq1 + q[inode];
					qq2 = qq2 + fabs(q[inode]);
					}//2000 continue
					qq1 = qq1 + Qsurfcorro1;
					qq2 = qq2 + fabs(Qsurfcorro1);
					qerror = qqbudget;//qerror = qq1/qq2*200; qq1
					qoutmax=qoutmax/qq2*200;

return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
