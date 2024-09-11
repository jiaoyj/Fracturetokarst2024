// FractureModel.h: interface for the FractureModel class.
//
//////////////////////////////////////////////////////////////////////



#define pii 57.29578  //原180.0/3.14159              (修改过)
#define rnum 12
#define macoun 499 //单个结点相连裂隙的最大值
#define	npmax 4999 //裂隙结点的最大值

class FractureModel  
{
public:
	FractureModel();
	virtual ~FractureModel();

	//平面流计算用到的变量及函数
	//double drn,dfact,dmodul,dconst;
private:
	double meanl[5],minl[5],maxl[5],meanb[5],minb[5],maxb[5];
	int    ndx,ndy;
	int    ndisl[5],ndisb[5],ndisd[5];
	double meand[5],mind[5],maxd[5],meanc[5],len;
	double sigmal[5],sigmab[5],sigmad[5];
	double minx,miny,maxx,maxy;
	double xy12[5000][5],xjoin[5000][500],b[5000],yjoin[5000][500],xbp[100],ybp[100];
	int    nupoin[5000][3],ijpoin[5000][500];
	int    nj[10],iiw[5000];
	double Qsurfcorro,Qrain,con_recharg; 
public:
	int    ngroup,napt,nconec,nodes;
	int    npb11,npb12,npb21,npb22,np1n1,np1n2;
	int    nboun1,nboun2,nbound; 
	int    numit,icount[5000],numitouter; int wettest01; int wettest02;
	double Remax,qoutmax,qerror,errmax;
	double xynp[5000][3],nijkl[5000][13], bnijkl[5000][13], q[1000],Qcell[5000][5],Re[5000][5000], Re0[5000][5000];
	double alpha[10],Tb[10],h[5000],hb1[100],qb2[100];
	double height[5000],qab2[1001][6];
	double newtondifAA[5000][5000]; int newtonii=0;
	
	int number_flow;

private:
	void poin();
	double amin1(double x1,double x2);
	double amax1(double x1,double x2);
	int outin(double x,double y);
	void queue();
	int given_no(int ni,int nj);
	void cell_parameter_calculate();
public:
	void inputdata();
	void generate_fracture();
	void connectedFracture();
	void plane_waterhead();
	void profile_waterhead();
    virtual void outputdata();

//剖面流计算用到的变量及函数 public
	int nbounf,npbf1,npbf2;

};


