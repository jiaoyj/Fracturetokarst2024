// RandomGenerator.cpp: implementation of the RandomGenerator class.
//
//////////////////////////////////////////////////////////////////////


#include <cstdlib>

#include "RandomGenerator.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RandomGenerator::RandomGenerator()
{

}

RandomGenerator::~RandomGenerator()
{

}
//---------------------------------------------
void RandomGenerator::init1(double seed_factor){
dmodul = (double) RAND_MAX*seed_factor;
return;
}
//---------------------------------------------
double RandomGenerator::rn(){
    double random_data;
  	double rrn;
    random_data = rand();
    rrn = fmod(random_data,dmodul);
    rrn = rrn/dmodul;
	return rrn;
	}
//-------------------- --------------------------- 
double RandomGenerator::unifrm(double a,double b){
	double random,rrn;
	rrn = rn();
	random = a + (b - a)* rrn;
	return random;
	}
//------------------------------------------
double RandomGenerator::gauss(double mean,double sigma,double min,double max){
    double v,r1,r2,random;
g100:;
	r1 = rn();
	r2 = rn();
	v=sqrt((-2.0*log(r1)))*cos(6.2832 * r2);
	random =  v* sigma + mean;
	if(random<min || random>max) goto g100;
	return random;
	}
//----------------------------------------------
double RandomGenerator::lognor(double mean,double sigma,double min,double max){
	double sigmx,sigmx2,r,v,random,meanx;
	sigmx2=log(sigma*sigma/(mean*mean) + 1.0);
	meanx=log(mean)-0.5*sigmx2;
	sigmx=sqrt(sigmx2);
l100:;
    r=rn();
	if(r==0.) goto l100;
	v=sqrt(-2.0*log(r))*cos(6.2832*r);
	random=exp(v*sigmx+meanx);
	if(random<min || random>max) goto l100;
	return random;
	}
//-------------------------------------------------
double RandomGenerator::erlang(double mean,double min,double max){
	double r,random,alpha;
        alpha=1./mean;
e100:;
    r=1.0;
	r=r*rn();
	if(r==0.) goto e100;
	random = -1.0/alpha*log(r);
	if(random<min || random>max) goto e100;
	return random;
	}
//----------------------------------------------