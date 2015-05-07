#include "molecule.h"
static void outputWelcome()
{
	cout << "##################################################" << endl;
	cout << "#      Welcome to use Molecule Recognition       #" << endl;
	cout << "# Many stable configurations will be produced    #" << endl;
	cout << "##################################################" << endl;
}
int main()
{
    outputWelcome();
    ifstream thisfile("task.txt",ios::in);
        string ifile1,ifile2;
        thisfile>>ifile1>>ifile2;
        double StepLength = 3.5;
	double StepPrecision = 0.15;
	double MaxRotTime = 10000;
	double RotPrecision = 20;
        double dT=10;
        thisfile>>StepLength>>StepPrecision>>MaxRotTime>>RotPrecision>>dT;
        thisfile.close();
	Molecule a, b;
	a.ReadFromXYZfile("InitiConfig/"+ifile1);
	a.output();
	a.PerformRandomRotEuler(20);
	b.ReadFromXYZfile("InitiConfig/"+ifile2);
	b.output();
	b.PerformRandomRotEuler(20);
	string basis = "6-31g";
        /* This part is suitable for single test*/
	string MinConfigName = "SaveConfigs/Min_";
	string RelativeMinConfigName = "SaveConfigs/RelaMin_";
	MonteCarlo(a, b, StepLength, StepPrecision,MaxRotTime, RotPrecision,basis,MinConfigName, RelativeMinConfigName,dT);
	return 0;
}
