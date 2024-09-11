
//#include "stdafx.h"

#include <ctime>
#include "KarstModel.h"

#include <iostream>
#include <armadillo>
using namespace arma;
//using namespace std;
int main()
{
    cout << "Karst Fracture Free Flow and Evolution" << endl << endl;
//std::std::std::

	time_t now_time, now_time01, now_time02;
	now_time = time(NULL);
//	cout << "now time(s): " << now_time << endl << endl;

	KarstModel *pKarstModel2 = new KarstModel;

	pKarstModel2->inputdata();
	cout << "KarstModel.inputdata() is over." << endl << endl;
	
	pKarstModel2->generate_fracture();
	cout << "KarstModel.generate_fracture() is over" << endl << endl;

	pKarstModel2->connectedFracture();
	std::cout << "KarstModel.ConnectedFracture() is over" << endl << endl;

	
	now_time01 = time(NULL);
	cout << "用时(s): " << (now_time01- now_time)<< endl << endl;


	pKarstModel2->profile_KarstEvolve();
	std::cout << "KarstModel.profile_KarstEvolve() is over" << endl << endl;

	now_time02 = time(NULL);
	cout << "用时(hour): " << (now_time02- now_time01)/3600.0<< endl << endl;

	delete pKarstModel2;
}

