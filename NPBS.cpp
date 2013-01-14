#include <cstdlib>
#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>


// to run: ./NPBS.exe (# of reads in OUTPUT_TEXT_FILE_FROM_Preprocess_Map_Files)
// example: ./NPBS.exe 1000000

using namespace std;

int main(int argc, char * argv[])
{
  
	char buf[BUFSIZ];
	if (argc < 2) {
		std::cerr << "num reads expected as an argument.\n";
		return 1;
	}
	long size = 0;
	size = std::atol(argv[1]);
	if (size < 1) {
		std::cerr << "Expecting num reads to be > 0.\n";
		return 1;
	}
	char* test;
	test = strtok(buf, " ");
	
	std::vector<int> array_a(size, 0);
	
	string line;
	ifstream infile // ("OUTPUT_TEXT_FILE_FROM_Preprocess_Map_Files");
	if (infile.is_open())
	{
		int i = 0;
		int j;
		
		cout << "Reading raw reads to memory." << endl;
		
		while ( infile.good() )
		{
			getline (infile, line);
			j = atoi ( line.c_str() );
			if (j != 0)
			{
				array_a[i] = j;
				i++;
			}
		}
		infile.close();
	}
	
	else cout << "Unable to open file";
	
	int sim;
	const int NUM_SIMS = 10;
	const int NUM_REPS = 100;
	const int NUM_GENES = 96943;
	for (sim=0; sim<NUM_SIMS; sim++)
	{
		std::stringstream ss;
		ss << //"LIBRARY_NAME_Bootstrap" << sim << ".csv";
		std::string name = ss.str();
		
		
		std::vector<int> counts(NUM_REPS, 0);
		std::vector<std::vector<int> >	array_b(NUM_GENES, counts);
		
		
		int h;
		srand ( time(NULL) );
		for (h=0; h<NUM_REPS; h++)
		{
			int g = 0;
			while (g < size)
			{
				cout << g << endl;
				array_b[array_a[(rand() % size)]-1][h] += 1.0;
				g++;
			}
		}
		
		int l;
		ofstream outfile;
		outfile.open( name.c_str() );
		for (l = 0; l < NUM_GENES; l++)
		{
			int b;
			for (b=0; b<NUM_REPS; b++)
			{
				outfile << array_b[l][b] << ",";
			}
			outfile << "\n";
		}
		outfile.close();
		printf("%i of %d Non-parametric bootstrap simulations completed. \n", (sim+1), NUM_SIMS);
	}
	return 0;
}


