#include "hmm.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int cti(char input){
	return (int(input)-65);
}

int main(int argc, char *argv[]){
	string modellist = argv[1];
	string testing_file_name = argv[2];
	string result_file_name = argv[3];
	const int model_num = 5;
	const int T = 50;
	const int state = 6;


	HMM hmms[model_num];
	load_models( "modellist.txt", hmms, model_num);
	// dump_models( hmms, model_num);

	double count = 0;
	string line;
	string line2;
	ifstream myfile(testing_file_name);
	ifstream answerFile("testing_answer.txt");

	ofstream resultFile;
	resultFile.open(result_file_name);

	if(myfile.is_open()){
		while(getline(myfile, line)){
			// cout<<line<<endl;
			double m_prob = -1;
			int m_best = -1;
			for(int m=0 ; m<model_num ; m++){
				double decode[T][state];
				// initial of decode (t = 0)
				for(int i=0 ; i<state ; i++){
					decode[0][i] = hmms[m].initial[i] * hmms[m].observation[cti(line[0])][i];
				}
				// induction of decode
				for(int t=1 ; t<T ; t++){
					// j 是後面
					for(int j=0 ; j<state ; j++){
						double max = -1;
						for(int i=0 ; i<state ; i++){
							double a = decode[t-1][i] * hmms[m].transition[i][j];
							if(a>max)
								max = a;
						}
						decode[t][j] = max * hmms[m].observation[cti(line[t])][j];
					}
				}

				double mmmp = 0;
				for(int i=0 ; i<state ; i++){
					mmmp += decode[T-1][i];
				}
				if(mmmp>m_prob){
					m_prob = mmmp;
					m_best = m+1;
				}
			}
			getline(answerFile, line2);
			string result = "model_0" + to_string(m_best) + ".txt";
			if(line2.compare(result)==0)
				count+=1;
			// cout << result << endl;
			// cout << line2 <<endl<<endl;
			resultFile << result << " " << m_prob << endl;
			// cout << "model_0" << m_best << ".txt" << endl;
		}
		// cout<<count<<endl;
		// cout<<count/2500<<endl;
		myfile.close();
	}	

	return 0;
}