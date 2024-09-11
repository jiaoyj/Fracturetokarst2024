// RandomGenerator.h: interface for the RandomGenerator class.
//
//////////////////////////////////////////////////////////////////////


class RandomGenerator  
{
private:
	double dmodul;

public:
	RandomGenerator();
	virtual ~RandomGenerator();

	void init1(double seed_factor); // change the "seed_factor" to generate different data series 
                                    // " 0.0 < seed_factor < 1.0"
	double rn();// 0-1 uniform random data                     
	double unifrm(double a,double b); //a-b uniform random data 
	double gauss(double mean,double sigma,double min,double max); // normal distributed random data
	double lognor(double mean,double sigma,double min,double max);// log_normal distributed data
	double erlang(double mean,double min,double max); // exponential distributed data
};
