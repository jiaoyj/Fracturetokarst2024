// KarstModel.h: interface for the KarstModel class.
//
//////////////////////////////////////////////////////////////////////



#include "FractureModel.h"

class KarstModel : public FractureModel  
{
public:
	KarstModel();
	virtual ~KarstModel();

private:
	int nl,ndt,idt,iwrite,iwr;  
	double t,dtt,k1,k2,k4,xmiddle,camiddle;   
	int nodel[20000][3],ijl1[20000],ijl2[20000];
	double ca[5000],bijkl[5000][5];
    void number(int jj1,int jj2,int nl,int &nw);
    void fkarst(double krate,double xmidu,double qin0, double bo,double ca0,double x,double &ca,double &baverag,double dtt);
	
public:
	void plane_KarstEvolve();
	void profile_KarstEvolve();
	void outputdata();
	int numit0; double errcmax;


};
