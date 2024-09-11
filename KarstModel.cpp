// KarstModel.cpp: implementation of the KarstModel class. 
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"

#include "KarstModel.h"
#include <math.h>
#include <armadillo>
using namespace arma;
//using namespace std;
#include <iostream>


#define	np_Max 19998 //nodes maximum

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

KarstModel::KarstModel()
{
	
}

KarstModel::~KarstModel()
{

}

//=====================================================================================
//=====================================================================================
//=====================================================================================
void KarstModel:: number(int jj1,int jj2,int nl,int &nw)
{
   int iiii2, ii1,ii2;
  nw=0;
  for(iiii2 =1; iiii2 <=nl; iiii2++)
  {
     ii1=nodel[iiii2][1];
	 ii2=nodel[iiii2][2];
	 if((ii1==jj1 && ii2==jj2) || (ii1==jj2 && ii2==jj1))
	 {
		 nw=2;
		 break;
	 } 
  }

return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void KarstModel:: fkarst(double krate,double xmidu,double qin0, double bo,double ca0,double x,double &ca,double &baverag,double dtt)
{
  double kc,CA1A,caeq,v,xmiddle,c0cmiddle,i;
  double dcca0,dca,gca,db,b;
  double dxdouble,dxmax,caup,dx,bb;
  double xxi,dt,rate,ca1,rat,c0c;

  caeq=2.00e-3;//3.35e-3mol/L, 134 mg/L
  //caeq=3.35e-4;//unit:mol/L，0.335mmol/L,13.4mg/L
  //caeq=2.0e-3;//unit:mol/L，2mmol/L,80mg/L
  kc=30.0;//unit:cm4/mol/s
  v=qin0/bo;//unit:cm/s
  

 //Method 0: ================================================================================
     // ++ rechaege fracture is unsoluable++++++
     //      if(nij.le.n1) kc=0.
     //----------- F=Kc(Ceq-c)**2-----
  
  if(krate==0)
  {  //if1                                 
    dt=x/v*100.0;//unit:t->s[T],x->m[L]

	dcca0=(caeq-ca0)/(1.0+2.0*kc*1.0e-3*dt*(caeq-ca0)/bo);
	ca=caeq-dcca0;
	dca=ca-ca0;

     //     vCa=vWater(cm+3)*dca(l)=Q*dtt=v*bo*H/1000.
	gca=31536.0*40*qin0*dca*dtt; // <<=>> gca=31536000.0*40.0*(qin0/1000)*dca*dtt   //86400*365=31356000
	db=gca/(x*100.0)/xmidu;  // <<=>> db=gca/(x*100.0*1.0)/xmidu
    //      db=31536000.*8170.*bo**3*dca*dtt/x**2/xmidu
	b=bo+db;
	baverag=b;
	return;
  } //end if0
  
 //Method 1:========================Laminar flow,1 and 2 order dissolution===========================================
  if(krate==1)
  {  //if 

	k1=4.0e-11;
	k2=4.0e-10;
	dca = 0;
	caup=ca0;//unit:mol/L caup
	c0c=caup/caeq;
	dt=x/v*100.0; //dt=x*100/(qin0/(bo*1))=x*bo/qin0*100//unit:s
	xmiddle = 0; int xnum = 0;

	if (c0c < 0.9)
	{
			dca = caeq - caup - (caeq - caup) * exp(-200000 * k1 * x / caeq / qin0);
			ca = caup + dca;
			c0c = ca / caeq;
		if(c0c >=0.9)
		{
			xmiddle = log(0.1 * caeq / (caeq - caup)) / (-200000 * k1 / caeq / qin0);
			caup = 0.9 * caeq;

			CA1A = 200000 * k2 * (x - xmiddle) * (caeq - caup)  + qin0 * caeq * caeq ;
			dca = caeq - caup - caeq *caeq *qin0* (caeq - caup) / CA1A;//unit:mol/cm2/s
			ca = caup + dca;
		}
	}	
	else 
		{
			CA1A = 200000 * k2 * (x - xmiddle) * (caeq - caup)  + qin0 * caeq * caeq ;
			dca = caeq - caup - caeq *caeq *qin0* (caeq - caup) / CA1A;//unit:mol/cm2/s
			ca = caup + dca;
		}
	
	gca=31536.0*100.0*qin0*dca*dtt; //  //86400*365*40/1000/100=12614.4   //86400*365=31356000
	db=gca/(x*100.0)/xmidu;//db=31536000.*8170.*bo**3*dca*dtt/x**2/xmidu
	
	b=bo+db;
	baverag=b;
	return;
  }//end if1

 //Method 2:========================Laminar flow,1 and 4 order dissolution===========================================
  if (krate == 2)
  {  //if 

	  k1 = 4.0e-12;
	  k4 = 4.0e-8;
	  dca = 0;
	  caup = ca0;//unit:mol/L caup
	  c0c = caup / caeq;
	  dt = x / v * 100.0; //dt=x*100/(qin0/(bo*1))=x*bo/qin0*100//unit:s
	  xmiddle = 0; int xnum = 0;

	  if (c0c < 0.9)
	  {
		  dca = caeq - caup - (caeq - caup) * exp(-200000 * k1 * x / caeq / qin0);
		  ca = caup + dca;
		  c0c = ca / caeq;
		  if (c0c >= 0.9)
		  {
			 /* i = 0;
			  do {
				  i = i + 1;
				  xmiddle = i * x / 1000;
				  dca = caeq - caup - (caeq - caup) * exp(-200000 * k1 * xmiddle / caeq / qin0);
				  ca = caup + dca;
				  c0c = ca / caeq;
				  xmiddle = log(0.1 * caeq / (caeq - caup)) / (-200000 * k1 / caeq / qin0);
			  } while (c0c < 0.9);
			  */
			  xmiddle = log(0.1 * caeq / (caeq - caup)) / (-200000 * k1 / caeq / qin0);
			  caup = 0.9 * caeq;

			  CA1A = pow(qin0 * caeq / (600000 * k4 * (x - xmiddle) * (caeq - caup) * (caeq - caup) * (caeq - caup) + qin0 * caeq * caeq * caeq * caeq), 1.0 / 3.0);
			  dca = caeq - caup - caeq * (caeq - caup) * CA1A;//unit:mol/cm2/s
			  ca = caup + dca;
		  }
	  }
	  else
	  {
		  CA1A = pow(qin0 * caeq / (600000 * k4 * x * (caeq - caup) * (caeq - caup) * (caeq - caup) + qin0 * caeq * caeq * caeq * caeq), 1.0 / 3.0);
		  dca = caeq - caup - caeq * (caeq - caup) * CA1A;//unit:mol/cm2/s
		  ca = caup + dca;
	  }

	  gca = 31536.0 * 100.0 * qin0 * dca * dtt; //
	  db = gca / (x * 100.0) / xmidu;//db=31536000.*8170.*bo**3*dca*dtt/x**2/xmidu

	  b = bo + db;
	  baverag = b;

	  if (isnan(b))
	  {
		  double deti222 = 0;
	  }
	  return;
  }
  //Method 3:========================Turbulent flow,1 and 4 order dissolution===========================================
 if(krate==3)
  {  

	k1=4.0e-11;
	k4=4.0e-7;
	dca = 0;
	caup=ca0;//unit:mol/L 
	c0c=caup/caeq;
	dt=x/v*100.0; //dt=x*100/(qin0/(bo*1))=x*bo/qin0*100//unit:s
	xmiddle = 0; int xnum = 0;

	if (c0c < 0.9)
	{
			dca = caeq - caup - (caeq - caup) * exp(-200000 * k1 * x / caeq / qin0);
			ca = caup + dca;
			c0c = ca / caeq;
		if(c0c >=0.9)
		{

			xmiddle = log(0.1 * caeq / (caeq - caup)) / (-200000 * k1 / caeq / qin0);
			caup = 0.9 * caeq;
			
			CA1A = pow(qin0 * caeq / (600000 * k4 * (x - xmiddle) * (caeq - caup) * (caeq - caup) * (caeq - caup) + qin0 * caeq * caeq * caeq * caeq), 1.0 / 3.0);
			dca = caeq - caup - caeq * (caeq - caup) * CA1A;//unit:mol/cm2/s
			ca = caup + dca;
		}
	}	
	else 
		{
			CA1A = pow(qin0 * caeq / (600000 * k4 * x * (caeq - caup) * (caeq - caup) * (caeq - caup) + qin0 * caeq * caeq * caeq * caeq), 1.0 / 3.0);
			dca = caeq - caup - caeq * (caeq - caup) * CA1A;//unit:mol/cm2/s
			ca = caup + dca;
		}
	
	gca=31536.0*100.0*qin0*dca*dtt; //
	db=gca/(x*100.0)/xmidu;//db=31536000.*8170.*bo**3*dca*dtt/x**2/xmidu
	
	b=bo+db;
	baverag=b;
	if (isnan(b))
	{
		double deti222 = 0;
	}
	return;
  }
//end if3 



  return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void KarstModel:: plane_KarstEvolve()
{	
    int i,j,k,i1,nw,jj1,jj2,nij;
	int maj,iterma,nup;
	double xmidu,b,qq,tol,err; 
	double caeq,cain,dcca,dh,dx,bo,bout,qin,qin1;
	double J,kbl,Vl,kbt,kb0,Re0;
	double hsort[5000][ 2];
	double dertac,newtoncm,newtondcm,sumcm,sumdcm;

	//----------------------------------
	FILE *fp1_input;  // Input karst evolvion timestep
    if((fp1_input=fopen("In_TimeStep.dat","r"))==NULL)
	{
		std::cout << "Fail to open In_TimeStep.dat\n";
     return;
	}
	fscanf(fp1_input,"%d%d%lf",&ndt,&iwr,&dtt); 
	fclose(fp1_input);
	
	caeq = 2.00e-3;//unit:mol/L
	maj=4 ;
	xmidu=2.5 ;

 
	t=0;

 for(idt=1;idt<=(ndt+1);idt++)//
{   //2500

	iwrite=idt % iwr ; 

	//-------------------------------water head computing------------------------------------
	this->plane_waterhead();



	iterma=1000;
	tol=1.0e-12;
	numit0=0; 


	for(i=1;i<=nodes;i++)
	{  //80
		for (j = 1; j <= 4; j++)
		{ 
			bijkl[i][j] = nijkl[i][j + 8];		

		}
	  ca[i] = 0.99 * caeq;
	}  //80 continue

	   //--------------------------------------concentration boundary and initial condition---------------------------------------
	for (i = 1; i <= nodes; i++)//i=n1+1
	{  //200
		nup = 0;
		for (j = 1; j <= 4; j++)
		{//100
			nij = int(nijkl[i][j]);
			if (nij > 0)			//
			{
				jj1 = j + 8;			
				dh = h[nij] - h[i];//只计算流入中心节点裂隙的溶蚀作用
				dx = nijkl[i][j + 4];
				bo = nijkl[i][jj1];
				qin1 = -Qcell[i][j];
			}//dh与Qcell[][]方向相反

			if (nij >= 1 && bo > 0.0 && qin1 > 1.0e-4)  
			{	
				nup = nup + 1;			

			}

		}

		if (nup >= 1)//内节点流量是负
		{//ca[i]=dcca/qq;
			int iuseless;
		}

		else//即非内节点，只有边界节点才满足这种情况，定流量的排泄边界点流量为负，是属于上面内节点的情况
		{
			//------------------定水头边界浓度赋值----------
			if (i <= npb12 && Qcell[i][4] >= 0)//高的定水头边界赋值浓度
			{
				ca[i] = 0.89 * caeq;
			}

			//--------------定流量边界节点----------

			if (i > npb12 && Qcell[i][4] > 0.001)
								
			{
				ca[i] = 0.89 * caeq;
			}


  
			 //--------------隔水边界----------
			if (i > npb12 && Qcell[i][4] <= 0.001)
			{
				ca[i] = caeq;
			}
		}

	}//结束节点i循环		
	   //--------------------------------------concentration boundary and initial condition---------------------------------------

	        for(i=1;i<=nodes;i++)	
			{        
				hsort[i][0]=i;      
				hsort[i][1]=h[i];	
			} 

			for(i=1;i<=nodes;i++)
			{ 
				for(j=i+1;j<=nodes;j++)
				{ 
					if ( hsort[j][1] > hsort[i][1] )
					{
						double temph,tempi;
						temph = hsort[i][1];
						tempi = hsort[i][0];
						hsort[i][1] = hsort[j][1];
						hsort[i][0] = hsort[j][0];
						hsort[j][1]=temph;
						hsort[j][0] = tempi;
					}
				}	
			}  

			//+++++++++++++++++++++++++++++++++ Major Calculation +++++++++++++++++++++++++++++++++

		int ireal;
	 for(i=1;i<=nodes;i++)//i=n1+1
	{  //200

		  nup=0;
		  qq=0.0;
		  dcca=0.0;
		  cain=0.0;
		  bout=0.0;
		  sumcm = 0;
		  sumdcm = 0;
		  
		  ireal =int( hsort[i][0]);
			for(j=1;j<=4;j++)
			{//100
				nij=int(nijkl[ireal][j]);
				if(nij>0)			
				{jj1=j+8;			
				dh=h[nij]-h[ireal];
				dx=nijkl[ireal][j+4];
				bo=nijkl[ireal][jj1];
				qin1=-Qcell[ireal][j];}
				
				if(nij>=1 && bo>0.0 && qin1> 1.0e-4) 
				{//开始if
				nup=nup+1;							

				if(fabs(Re[ireal][nij])<2300)
					{this->fkarst(2,xmidu,qin1,bo,ca[nij],dx,cain,bout,dtt);}//层流的溶蚀速率计算  cain是单裂隙出口浓度
				else
					{this->fkarst(3,xmidu,qin1,bo,ca[nij],dx,cain,bout,dtt);}//紊流的溶蚀速率

				  qin=qin1;

				  bijkl[ireal][j]=bout;					
				for(k=1;k<=4;k++)
				{
					i1=int (nijkl[nij][k]);
					if(i1==ireal)bijkl[nij][k]=bout;
				}

				dertac=cain-ca[nij];
				newtoncm=qin*ca[nij];
				

				newtondcm=qin*dertac;	
				sumcm=sumcm+newtoncm;
				sumdcm=sumdcm+newtondcm;
				qq=qq+qin;
				dcca=dcca+cain*qin;
			  }//结束if
		
			}//结束j循环

				if(nup>=1)
				{ca[ireal]=dcca/qq;}
				

	}//结束节点i循环
	



	//----------------------outputdata------------------------------
	this->outputdata();
	//---------------------      !    ------------------------------

	t = t + dtt;


	for(i=1;i<=nodes;i++)
	   for(j=1;j<=4;j++)
		   nijkl[i][j+8]=bijkl[i][j];


 } //2500 

  //+++++++++++++++++++++++++++++++++ Major Calculation +++++++++++++++++++++++++++++++++

 return;
}


//=====================================================================================
//=====================================================================================
//=====================================================================================
void KarstModel:: profile_KarstEvolve()
{

    int i,j,k,i1,nw,jj1,jj2,jj3,nij;
	int maj,iterma,nup;
	double xmidu,b,qq,tol,err;
	double hsort[5000][2];
	double caeq,cain,carain,dcca,dh,dx,bo,bout,qin,qin1;
	double J,kbl,Vl,kbt,kb0,Re0;
	//----------------------------------
	FILE *fp1_input;  // Input karst evolvion timestep
    if((fp1_input=fopen("In_TimeStep.dat","r"))==NULL)
	{
		std::cout << "Fail to open In_TimeStep.dat\n";
     return;
	}
	fscanf(fp1_input,"%d%d%lf",&ndt,&iwr,&dtt); 
	fclose(fp1_input);
	
	caeq = 2.00e-3;//unit:mol/L
	maj=4 ;
	xmidu=2.5 ;
	iwrite=0;
		
	double dertac,newtoncm,newtondcm,sumcm,sumdcm;
 
	t=0;
 

  for (idt = 1; idt <= (ndt + 1); idt++)
  {//2500

	  iwrite = idt % iwr;


	  //-------------------------------water head---------------------------------------------------------------
	  this->profile_waterhead();
	  	std::cout <<idt << " times evolve profile waterhead() is over" << endl ;
	  //--------------------------------------------------------------------------------------------------------

	  //----------------------outputdata------------------------------
		this->outputdata();
		//---------------------      !    ------------------------------




	  iterma = 1000;
	  tol = 1.0e-9;
	  numit0 = 0;
	  for (i = 1; i <= nodes; i++)
	  {  //80
		  for (j = 1; j <= 4; j++)bijkl[i][j] = nijkl[i][j + 8];
	  }  //80 continue



//------------------------------------------concentration boundary--------------------------------------------

	  for (i = 1; i <= nodes; i++)

	  { //220

		  ca[i] = 0.99 * caeq;
		  if (height[i] > h[i])

		  {
			  ca[i] = 0;
		  }

	  }  //220		   


	  for (i = 1; i <= nodes; i++)//i=n1+1
		  if (height[i] <= h[i])
		  {//678					

			

			  nup = 0;
			  qq = 0.0;
			  dcca = 0.0;
			  cain = 0.0;
			  bout = 0.0;
			  sumcm = 0;
			  sumdcm = 0;


			  for (j = 1; j <= 4; j++)
			  { //689
				  nij = int(nijkl[i][j]);
				  if (nij == 0)continue;


				  if (height[nij] <= h[nij])
				  {//704

					  nij = int(nijkl[i][j]);
					  if (nij > 0)			
					  {
						  jj1 = j + 8;			
						  dh = h[nij] - h[i];
						  dx = nijkl[i][j + 4];
						  bo = nijkl[i][jj1];
						  qin1 = -Qcell[i][j];
					  }//dh与Qcell[][]方向相反
					  if (nij >= 1 && bo > 0.0 && qin1 > 1.0e-4) 
					  {
						  nup = nup + 1;

					  }


				  }//704


				  if (height[nij] > h[nij])
					  for (k = 1; k <= number_flow; k++)
						  if (i == int(qab2[k][1]) && nij == int(qab2[k][2]))
						  {
							  ca[nij] = 0.98 * caeq;
							  
							  ////----------------------rain water near the water table---------------
						  }


			  }//结束j循环


			  if (nup >= 1)
			  {
				  int useless;//ca[i] = dcca / qq;
			  }

			  else
			  {
				  //--------------第一类定水头边界浓度赋值（流量大于0）----------
				  if (i <= npb12 && Qcell[i][4] >= 0)
				  {
					  ca[i] = 0.85 * caeq;
				  }
				  //--------------第二类给定流量边界节点----------
				  if (i > npb12 && i <= npb22 && Qcell[i][4] > 0.000001)
				  {
					  ca[i] = 0.99 * caeq;
				  }
  
				  //--------------第二类隔水边界----------
				  if (i > npb12 && i <= npb22 && fabs(Qcell[i][4]) <= 0.000001)
				  {
					  ca[i] = caeq;
				  }
			  }

		  }//678 



								  for (i = 1; i <= nodes; i++)
								  {
									  hsort[i][0] = i;
									  hsort[i][1] = h[i];
								  }
							
								  for (i = 1; i <= nodes; i++)
								  {
									  for (j = i + 1; j <= nodes; j++)
									  {
										  if (hsort[j][1] > hsort[i][1])
										  {
											  double temph, tempi;
											  temph = hsort[i][1];
											  tempi = hsort[i][0];
											  hsort[i][1] = hsort[j][1];
											  hsort[i][0] = hsort[j][0];
											  hsort[j][1] = temph;
											  hsort[j][0] = tempi;
										  }
									  }
								  }
							
								  //+++++++++++++++++++++++++++++++++ Major Calculation +++++++++++++++++++++++++++++++++


	  int ireal;
	  for (i = 1; i <= nodes; i++)//i=n1+1
	  {	 		  
		   ireal = int(hsort[i][0]);

		  if (height[ireal] <= h[ireal])
		  {//678					

		  

			  nup = 0;
			  qq = 0.0;
			  dcca = 0.0;
			  cain = 0.0;
			  bout = 0.0;
			  sumcm = 0;
			  sumdcm = 0;

			 
			  for (j = 1; j <= 4; j++)
			  { //689
				  nij = int(nijkl[ireal][j]);
				  if (nij == 0)continue;


				  if (height[nij] <= h[nij])
				  {//704

					  nij = int(nijkl[ireal][j]);
					  if (nij > 0)			
					  {
						  jj1 = j + 8;			
						  dh = h[nij] - h[ireal];
						  dx = nijkl[ireal][j + 4];
						  bo = nijkl[ireal][jj1];
						  qin1 = -Qcell[ireal][j];
					  }//dh与Qcell[][]方向相反

					  if (nij >= 1 && bo > 0.0 && qin1 > 1.0e-4)  
					  {//698开始if
						  nup = nup + 1;


						  if (fabs(Re[ireal][nij]) < 2300)
						  {
							  this->fkarst(2, xmidu, qin1, bo, ca[nij], dx, cain, bout, dtt);
						  }//层流的溶蚀速率计算
						  else
						  {
							  this->fkarst(3, xmidu, qin1, bo, ca[nij], dx, cain, bout, dtt);
						  }//紊流的溶蚀速率


						  bijkl[ireal][j] = bout;
						  qin = qin1;

						  for (k = 1; k <= 4; k++)
						  {//83
							  i1 = int(nijkl[nij][k]);
							  if (i1 == ireal)bijkl[nij][k] = bout;
						  }

						  dertac = cain - ca[nij];
						  newtoncm = qin * ca[nij];

						  newtondcm = qin * dertac;
						  sumcm = sumcm + newtoncm;
						  sumdcm = sumdcm + newtondcm;
						  qq = qq + qin;
						  dcca = dcca + cain * qin;
					  }//698结束if

				  }//704


				  if (height[nij] > h[nij])
					  for (k = 1; k <= number_flow; k++)
						  if (ireal == int(qab2[k][1]) && nij == int(qab2[k][2]))
						  {


							  newtoncm = qab2[k][3] * ca[nij];
							  newtondcm = 0;
							  sumcm = sumcm + newtoncm;
							  sumdcm = sumdcm + newtondcm;
							  qq = qq + qab2[k][3];
							  dcca = dcca + qab2[k][3] * ca[nij];

						  }



			  }//689 j


			  for (j = 1; j <= 4; j++)
			  {
				  nij = int(nijkl[ireal][j]);

				  if (height[nij] > h[nij])
					  for (k = 1; k <= number_flow; k++)
						  if (ireal == int(qab2[k][1]) && nij == int(qab2[k][2]))
						  {
							  ca[ireal] = dcca / qq;
						  }

			  }

			  if (nup >= 1)
			  {
				  ca[ireal] = dcca / qq;
			  }




		  }//678  

	  }//i



	for(i=1;i<=nodes;i++)
	   for(j=1;j<=4;j++)
		   nijkl[i][j+8]=bijkl[i][j];


	t = t + dtt;



  } //2500 continue
//+++++++++++++++++++++++++++++++++ Major Calculation +++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++  !  +++++++++++++++++++++++++++++++++++++++++
return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================
void KarstModel:: outputdata()
{
	int i,i1,i2,j,j1,idirect,maj(4),k;
	double x1,x2,y1,y2,ideta,detah,b,qs,sumqs(0.0),L,Rey;

	FILE *fp1_output;  // For: 'karstLine.dat' (output file)
	FILE *fp2_output;  // For: 'q.dat' (output file)
	FILE *fp3_output;  // For: 'Cca2.dat' (output file)
	FILE *fp4_output;  // temp (output file)
	FILE *fp5_output;  // Xsection flow line (output file)
	FILE *fp6_output;  // Ysection flow line (output file)
	FILE *fp8_output; 
	if(idt==1)
	{
		if((fp1_output=fopen("Out_KarstWaterFlux.dat","w"))==NULL)
		 {
			std::cout << "fail to open KarstLine.dat\n";
		 return;
		 }
   
		if((fp2_output=fopen("Out_KarstSolute.dat","w"))==NULL)
		 {
			std::cout << "fail to open KarstFlux.dat\n";
		return;
		 }
    
		if((fp3_output=fopen("Out_KarstLine.dat","w"))==NULL)
		 {
			std::cout << "fail to open KarstSolute.dat\n";
		return;
		 }
		
		if((fp4_output=fopen("Out_KarstTempOutput.dat","w"))==NULL)
		 {
			std::cout << "fail to open KarstSolute.dat\n";
		return;
		 }

		if((fp5_output=fopen("Out_KarstLineXSection.dat","w"))==NULL)
		 {
			std::cout << "fail to open KarstLineXSection.dat\n";
		return;
		 }

		if((fp6_output=fopen("Out_KarstLineYSection.dat","w"))==NULL)
		 {
			std::cout << "fail to open KarstLineYSection.dat\n";
		return;
		 }
         
         if(  (fp8_output=fopen("Out_Karstaperture.dat","w"))==NULL)
		 {
			std::cout << "fail to open Karstaperture.dat\n";
		return;
		 }
         
	}
    
	else
    
	{
		if((fp1_output=fopen("Out_KarstWaterFlux.dat","a"))==NULL)
		 {
			std::cout << "fail to open KarstLine.dat\n";
		 return;
		 }
   
		if((fp2_output=fopen("Out_KarstSolute.dat","a"))==NULL)
		 {
			std::cout << "fail to open KarstFlux.dat\n";
		return;
		 }
    
		if((fp3_output=fopen("Out_KarstLine.dat","a"))==NULL)
		 {
			std::cout << "fail to open KarstSolute.dat\n";
		return;
		 }
		
		if((fp4_output=fopen("Out_KarstTempOutput.dat","a"))==NULL)
		 {
			std::cout << "fail to open KarstSolute.dat\n";
		return;
		 }

		if((fp5_output=fopen("Out_KarstLineXSection.dat","a"))==NULL)
		 {
			std::cout << "fail to open KarstLineXSection.dat\n";
		return;
		 }

		if((fp6_output=fopen("Out_KarstLineYSection.dat","a"))==NULL)
		 {
			std::cout << "fail to open KarstLineYSection.dat\n";
		return;
		 }
         
         if(  (fp8_output=fopen("Out_Karstaperture.dat","a"))==NULL)
		 {
			std::cout << "fail to open Karstaperture.dat\n";
		return;
		 }
         
	}

			
	
	
	
	
	//-------------------------------water budget output------------------------------------	
    if(idt==1)
	{ 
	  fprintf(fp1_output,"FLOW INTO: Q[i](ml/s)\n");
      fprintf(fp1_output,"             ");
	  for(i=1;i<=npbf2;i++)
	  {
		  if(i<10)fprintf(fp1_output,"     Q[%d]  ",i);
			 else if(i<100)fprintf(fp1_output,"     Q[%d] ",i);
			    else if(i<1000)fprintf(fp1_output,"     Q[%d]",i);
			       else fprintf(fp1_output,"    Q[%d]   ",i);
	  }

	}

  {
		fprintf(fp1_output, "t=  %-8.1lf", t);	
		double qqboundout = 0; double qqboundin = 0;
		for (i = 1; i <= npbf2; i++)
		{
			fprintf(fp1_output, "  %8.5e", q[i]);
			qqboundout = qqboundout + q[i];
		}
		for (i = 1; i <= number_flow; i++)
		{
			qqboundin = qqboundin + qab2[i][3];

		}
    fprintf(fp1_output,"    qqboundout= %10.5e,qqboundin= %10.5e, 水头误差=%10.5e", qqboundout, qqboundin, errmax);

	fprintf(fp1_output, "      外迭代=%2d,内迭代=%2d,wettest01=%2d,wettest02=%2d\n",numitouter, numit,wettest01, wettest02);
	

  }




	//-------------------------------Ca2+----------------------------------------  
   //----------------------------------------------------------------------------
	 

	 if (iwrite == 0)
	 {
		 //if (t == 10000 || t==50000)
		 {
			 fprintf(fp2_output, "t=  %-8.1lf\n", t);
			 for (i = 1; i <= nodes; i++)
			 {
				 x1 = xynp[i][1];
				 y1 = xynp[i][2];
				 fprintf(fp2_output, " %10.5lf  %10.5lf  %12.8lf \n", x1, y1, ca[i]);
			 }//fprintf(fp2_output, "\n");
		 } //fprintf(fp2_output,"\n");
	 }
  


   //---------------------------------------------------------------------------------------------------------------   
   //-------------------fp3_output "Out_KarstLine.dat"--------------------------------------------------------------
  //---------------------------------------------------------------------------------------------------------------
   if(idt==1)fprintf(fp3_output,"%lf\n",ndt*dtt);
   if(iwrite==0)
   { // if

	// fprintf(fp5_output,"t=%-8.1lf\n",t);
	//fprintf(fp6_output,"t=%-8.1lf\n",t);
	 
	 

	 int  nw, jj1, jj2;	
	 double  b;

	 nl = 0; nw = 0;
	 for (i = 0; i <= 19999; i++)
	 {  
				 ijl1[i] = 0;
				 ijl2[i] = 0;
			 for (j = 0; j <= 2; j++)					 
				 {						 
					 nodel[i][j] = 0;				 
				 }
	 } 	 
	 
	 for (i = 1; i <= nodes; i++)
	 {  //220
		 for (j = 1; j <= 4; j++)
		 { //210

			 jj1 = i;
			 jj2 = (int)(nijkl[i][j]);
			 if (jj2 < 1) continue;
			 if (jj1 < np1n1 && jj2 < np1n1) continue;

			 b = nijkl[i][j + 8];
			 if (b <= 0) continue;
			 else number(jj1, jj2, nl, nw);

			 if (nw >= 1) continue;
			 else
			 {
				 nl = nl + 1;
				 nodel[nl][1] = jj1;
				 nodel[nl][2] = jj2;
				 ijl1[nl] = i;
				 ijl2[nl] = j;
			 }
		 }  //210	continue

	 }  //220	co	 	 
	 
	 
     fprintf(fp3_output,"%10.1lf %d\n",t,nl);
	 
	 
     ideta=0.00001;
	 for(i=1;i<=nl;i++)
	 {//260
	    i1=ijl1[i];
        j1=ijl2[i];
        i2=int (nijkl[i1][j1]);
		L=nijkl[i1][j1+4];
        x1=xynp[i1][1];
        y1=xynp[i1][2];
        x2=xynp[i2][1];
        y2=xynp[i2][2];
	    detah=h[i1]-h[i2];
		b= nijkl[i1][j1 + 8];
		if(h[i1]!=0 && h[i2]!=0)
		{
		Rey=detah/L*980*b*b*b/6/0.01/0.01;
		}
		else
		{
			Rey=0;
		}
		if(fabs(detah)<=ideta)idirect=0;
		  else 
		  {
		    if(detah>0)idirect=1;
			  else idirect=-1;
		  }

		
        fprintf(fp3_output,"%10.2lf %10.2lf %10.2lf %10.2lf %10.6lf %10.5lf %10.5lf %3d %15.4lf %3d %3d\n",x1,y1,x2,y2,b,h[i1],h[i2],idirect,Re[i1][i2],i1,i2);
		
	 }//260 continue
   } // end if
  //---------------------------------------------------------------------------------------------------------------   
   //-------------------------------------fp3_output "Out_KarstLine.dat"-------------------------------------------
  //---------------------------------------------------------------------------------------------------------------

  
  
  
   //----------------------------------------------------------------------------- 
	if(iwrite==0)
   {
     
	 fprintf(fp4_output,"t=%-8.1lf\n",t);
 
	 for(i=npbf1;i<=npbf2;i++)
	 {
		qs=Qcell[i][4];
		fprintf(fp4_output,"No%-10d%12.5lf\n ",i,qs);
		sumqs=sumqs+qs;//total flow of spring
	 }
	 fprintf(fp4_output," Sumqs=    %12.5lf ",sumqs);
	
	 //for(i=1;i<=nodes;i++)fprintf(fp2_output,"%12.5lf ",h[i]);
	 //fprintf(fp2_output,"\n");
	 fprintf(fp4_output,"\n");
   }


	if(iwrite==0)
   {
     
	// fprintf(fp5_output,"t=%-8.1lf\n",t);   96,236，404
		for (i =1; i <= 1; i++)
		{
			for (j = 1; j <= 4; j++)
			{
				if (h[i] == 0)continue;
				k = int(nijkl[i][j]);
				if (k < 1)continue;
				if (i < np1n1 && k < np1n1)continue;

				//L = nijkl[i][j + 4];
				//detah = h[i] - h[k];
				//b = bijkl[i][k];
				//if (h[i] != 0 && h[k] != 0){Re = detah / L * 980 * b * b * b / 6 / 0.01 / 0.01;}
				//else{ Re = 0;}

				fprintf(fp5_output, "t=  %-8.1lf %d %d  %10.4lf %10.4lf %10.4lf %12.8lf", t, i, k, fabs(Qcell[i][j]),fabs(Re[i][k]),nijkl[i][j+8], ca[i]);
				//fprintf(fp5_output,"t=%-8.1lf %d %d %10.4lf\n",t,i,k,Qcell[i][j]);
			}
		}
		for (i = 2; i <= 2; i++)
		{
			for (j = 1; j <= 4; j++)
			{
				if (h[i] == 0)continue;
				k = int(nijkl[i][j]);
				if (k < 1)continue;
				if (i < np1n1 && k < np1n1)continue;

				//L = nijkl[i][j + 4];
				//detah = h[i] - h[k];
				//b = bijkl[i][k];
				//if (h[i] != 0 && h[k] != 0) { Re = detah / L * 980 * b * b * b / 6 / 0.01 / 0.01; }
				//else { Re = 0; }

				fprintf(fp5_output, "   %d %d  %10.4lf %10.4lf %10.4lf %12.8lf", i, k, fabs(Qcell[i][j]), fabs(Re[i][k]), nijkl[i][j + 8], ca[i]);
				//fprintf(fp5_output,"t=%-8.1lf %d %d %10.4lf\n",t,i,k,Qcell[i][j]);
			}
		}
		for (i = 3; i <= 3; i++)
		{
			for (j = 1; j <= 4; j++)
			{
				if (h[i] == 0)continue;
				k = int(nijkl[i][j]);
				if (k < 1)continue;
				if (i < np1n1 && k < np1n1)continue;

				//L = nijkl[i][j + 4];
				//detah = h[i] - h[k];
				//b = bijkl[i][k];
				//if (h[i] != 0 && h[k] != 0) { Re = detah / L * 980 * b * b * b / 6 / 0.01 / 0.01; }
				//else { Re = 0; }

				fprintf(fp5_output, "   %d %d  %10.4lf %10.4lf %10.4lf %12.8lf", i, k, fabs(Qcell[i][j]), fabs(Re[i][k]), nijkl[i][j+8], ca[i]);
				//fprintf(fp5_output,"t=%-8.1lf %d %d %10.4lf\n",t,i,k,Qcell[i][j]);
			}
		}
		for (i = 4; i <= 21; i++)
		{
			for (j = 1; j <= 4; j++)
			{
				if (h[i] == 0)continue;
				k = int(nijkl[i][j]);
				if (k < 1)continue;
				if (i < np1n1 && k < np1n1)continue;

				//L = nijkl[i][j + 4];
				//detah = h[i] - h[k];
				//b = bijkl[i][k];
				//if (h[i] != 0 && h[k] != 0) { Re = detah / L * 980 * b * b * b / 6 / 0.01 / 0.01; }
				//else { Re = 0; }

				fprintf(fp5_output, "   %d %d  %10.4lf %10.4lf %10.4lf %12.8lf", i, k, fabs(Qcell[i][j]), fabs(Re[i][k]), nijkl[i][j+8], ca[i]);
				//fprintf(fp5_output,"t=%-8.1lf %d %d %10.4lf\n",t,i,k,Qcell[i][j]);// bijkl[i][j]
			}
		}
				fprintf(fp5_output, "\n");


	 }
	 
	 fprintf(fp8_output,"t=%-8.1lf\n",t);
	   for (i = 1; i <= nodes; i++)
	  {  //80
		  for (j = 1; j <= 4; j++)  fprintf(fp8_output, "  %10.5lf ", nijkl[i][j + 8]);
		  
		  fprintf(fp8_output, "\n");
	  }
	 
	//---------------------------------------------------------------------------
	//-------------------------------    !   ------------------------------------
	
		fclose(fp1_output);
		fclose(fp2_output);
		fclose(fp3_output);
		fclose(fp4_output);
		fclose(fp5_output);
		fclose(fp6_output);
				fclose(fp8_output);
return;
}
//=====================================================================================
//=====================================================================================
//=====================================================================================