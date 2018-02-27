#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

void AggregateTxt(std::string inFileName, std::string outFileName, int AggregationLimit)
{
	std::ifstream             infile(inFileName.c_str());
	std::ofstream             outfile;
	std::map<long int,double> mTimeVd;
	std::vector<long int>     avTime;
	std::vector<double>       avVD;
	long int t,avT=0;
	double   vd,avV=0.;

	while (infile >> t >> vd)
		{
			mTimeVd[t] = vd;
		}

	outfile.open(outFileName.c_str());
  
	for (std::map<long int,double>::iterator it=mTimeVd.begin(); it!=mTimeVd.end(); it++)
		{
			// outfile << it->first << " " << it->second << std::endl;
			if (avTime.size()<AggregationLimit){
				avTime.push_back(it->first);
				avVD.push_back(it->second);
			}
			if (avTime.size()==AggregationLimit){
				for(std::vector<long int>::iterator it1 = avTime.begin(); it1 != avTime.end(); ++it1){
					avT += *it1;
				}
				for(std::vector<double>::iterator it2 = avVD.begin(); it2 != avVD.end(); ++it2){
					avV += *it2;
				}
				avT /= avTime.size();
				avV /= avVD.size();

				outfile << avT << " " << avV << std::endl;

				avT = 0;
				avV = 0.;
				avTime.clear();
				avVD.clear();
			}
		}
	outfile.close();
  
}

int main(int argc, char** argv)
{
	std::string inFileName;
	std::string outFileName;
	int         AggLim = 60;
  
	//-------------------------------USAGE-----------------------------------------//
	if (argc<3){
		return 1;
	}

	for(int i=1; i<argc; i++){
		if(std::string(argv[i]) != "-i" &&
		   std::string(argv[i]) != "-o"){
			std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " <<  argv[i] << std::endl;
			std::cout << "Aggregator finished with exit code " << "2" << std::endl;
			return 2;
		} else {
			if(std::string(argv[i]) == "-i" && i!=argc-1) {
				inFileName = argv[++i];
				continue;
			}
			if(std::string(argv[i]) == "-i" && i==argc-1) {
				std::cerr << "\n[ERROR]: Input file name was not specified " << std::endl;
				std::cout << "Aggregator finished with exit code 3" << std::endl;
				return 3;
			}
			if(std::string(argv[i]) == "-o" && i!=argc-1) {
				outFileName = argv[++i];
				continue;
			}
			if(std::string(argv[i]) == "-o" && i==argc-1) {
				std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
				std::cout << "Aggregator finished with exit code 4" << std::endl;
				return 4;
			}
		}
	}
	//-------------------------------USAGE-----------(end)-------------------------//

	AggregateTxt(inFileName,outFileName,AggLim);

}
