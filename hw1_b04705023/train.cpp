#include "hmm.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int cti(char input){
	return (int(input)-65);
}

int main(int argc, char *argv[]){
	int iteration = atoi(argv[1]);
	string model = argv[2];
	string seq_model = argv[3];
	string save_model = argv[4];

	HMM hmm;
	loadHMM( &hmm, model.c_str());
	dumpHMM( stderr, &hmm );

	const int T = 50;
	const int state = hmm.state_num;
	const int observ = hmm.observ_num;



int count = 0;
while(count<iteration){
	count++;
	cout<<"===================="<<count<<"================"<<endl;

	double updatePi[state];
	double updateAup[state][state];
	double updateAdown[state];
	double updateBup[state][state];
	double updateBdown[state];
	for(int i=0 ; i<state ; i++){
		updatePi[i] = 0;
		updateAdown[i] = 0;
		updateBdown[i] = 0;
		for(int j=0 ; j<state ; j++){
			updateAup[i][j] = 0;
			updateBup[i][j] = 0;
		}
	}

	string line;
	ifstream myfile(seq_model);
	if(myfile.is_open()){
		while(getline(myfile, line)){
			// cout<<line<<endl;
			// declare alpha, beta, prob
			double alpha[T][state];
			double beta[T][state];
			double prob[T][state];
			double gamma[T][state];
			double eps[T-1][state][state];

			// ==========ALPHA==========
			// alpha initialization
			for(int j=0 ; j<state ; j++){
				alpha[0][j] = hmm.initial[j] * hmm.observation[cti(line[0])][j];
			}
			// alpha induction
			for(int i=1 ; i<line.length() ; i++){
				for(int j=0 ; j<state ; j++){
					alpha[i][j] = 0;
					for(int k=0 ; k<state ; k++){
						alpha[i][j] += (alpha[i-1][k] * hmm.transition[k][j]);
					}
					alpha[i][j] *= hmm.observation[cti(line[i])][j];
				}
			}
			// cout<<endl;

			// ==========BETA==========
			// beta initialization
			for(int j=0 ; j<state ; j++){
				beta[T-1][j] = 1;
			}
			// beta induction
			for(int i=T-2 ; i>=0 ; i--){
				for(int j=0 ; j<state ; j++){
					beta[i][j] = 0;
					for(int k=0 ; k<state ; k++){
						beta[i][j] += (beta[i+1][k] * hmm.transition[j][k] * hmm.observation[cti(line[i+1])][k] );
					}
				}
			}

			// ==========PROB==========
			for(int i=0 ; i<T ; i++){
				for(int j=0 ; j<state ; j++){
					prob[i][j] = alpha[i][j] * beta[i][j];
				}
			}

			// ==========GAMMA==========
			for(int i=0 ; i<T ; i++){
				double sum_state = 0;
				for(int j=0 ; j<state ; j++){
					sum_state += prob[i][j];
				}
				for(int j=0 ; j<state ; j++){
					gamma[i][j] = (prob[i][j]/sum_state);
				}
			}

			// ==========pre_EPSILON==========
			double pre_eps[T-1][state][state];
			for(int t=0 ; t<T-1 ; t++){
				for(int i=0 ; i<state ; i++){
					for(int j=0 ; j<state ; j++){
						pre_eps[t][i][j] = alpha[t][i]
										* hmm.transition[i][j] //Aij
										* hmm.observation[cti(line[t+1])][j] // Bi(Ot+1)
										* beta[t+1][j];
					}
				}
			}
			// ==========EPSILON==========
			for(int t=0 ; t<T-1 ; t++){
				double sum_eps = 0;
				for(int i=0 ; i<state ; i++){
					for(int j=0 ; j<state ; j++){
						sum_eps += pre_eps[t][i][j];
					}
				}
				for(int i=0 ; i<state ; i++){
					for(int j=0 ; j<state ; j++){
						eps[t][i][j] = (pre_eps[t][i][j]/sum_eps);
					}
				}
			}

			// // ==========check check==========
			// if(line=="ABCCCACDCFFACFFBCCDCCCBDEFBCCCCDCCBFACCCCCDBFFFFFF"){
			// 	cout<<"==========alpha"<<endl;
			// 	for(int i=0 ; i<state ; i++){
			// 		cout<<alpha[0][i]<<" "<<alpha[1][i]<<endl;
			// 	}
			// 	cout<<"==========beta"<<endl;
			// 	for(int i=0 ; i<state ; i++){
			// 		cout<<beta[0][i]<<" "<<beta[1][i]<<endl;
			// 	}
			// 	cout<<"==========prob"<<endl;
			// 	for(int i=0 ; i<state ; i++){
			// 		cout<<prob[0][i]<<" "<<prob[1][i]<<endl;
			// 	}
			// 	cout<<"==========gamma"<<endl;
			// 	for(int i=0 ; i<state ; i++){
			// 		cout<<gamma[0][i]<<" "<<gamma[1][i]<<endl;
			// 	}
			// 	cout<<"==========epsilon"<<endl;
			// 	for(int i=0 ; i<state ; i++){
			// 		cout<<eps[0][0][i]<<" "<<eps[0][1][i]<<endl;
			// 	}
			// }

			// ==========update Pi===========
			for(int i=0 ; i<state ; i++){
				updatePi[i] += gamma[0][i];
			}

			// ==========update A==========
			for(int i=0 ; i<state ; i++){
				// up
				for(int j=0 ; j<state ; j++){
					for(int t=0 ; t<T-1 ; t++){
						updateAup[i][j] += eps[t][i][j];
					}
				}
				// down
				for(int t=0 ; t<T-1 ; t++){
					updateAdown[i] += gamma[t][i];
				}
			}

			// ==========update B==========
			for(int i=0 ; i<state ; i++){
				double p[observ];
				for(int g=0 ; g<observ ; g++)
					p[g] = 0;

				double p2 = 0;
				for(int t=0 ; t<T ; t++){
					p[cti(line[t])] += gamma[t][i];
					p2 += gamma[t][i];
					// updateBup[cti(line[t])][i] += gamma[t][i];
				}

				updateBdown[i] += p2;
				for(int k=0 ; k<observ ; k++){
					updateBup[k][i] += p[k];
				}
			}

		}
		myfile.close();
	}

	// update after iterating all the data
	for(int i=0 ; i<state ; i++){
		hmm.initial[i] = updatePi[i] / 10000;
	}
	for(int i=0 ; i<state ; i++){
		for(int j=0 ; j<state ; j++){
			hmm.transition[i][j] = updateAup[i][j] / updateAdown[i];
		}
	}
	for(int i=0 ; i<observ ; i++){
		for(int j=0 ; j<state ; j++){
			hmm.observation[i][j] = updateBup[i][j] / updateBdown[j];
		}
	}


	dumpHMM( stderr, &hmm );
}



	ofstream modelfile;
	modelfile.open(save_model);
	modelfile << "initial: " << state <<endl;
	for(int i=0 ; i<state ; i++){
		if(i!=state-1)
			modelfile<<hmm.initial[i]<<" ";
		else
			modelfile<<hmm.initial[i]<<endl;
	}
	modelfile << endl;

	modelfile << "transition: " << state << endl;
	for(int i=0 ; i<state ; i++){
		for(int j=0 ; j<state ; j++){
			if(j!= state-1)
				modelfile << hmm.transition[i][j] <<" ";
			else
				modelfile << hmm.transition[i][j] <<endl;
		}
	}
	modelfile << endl;
	modelfile << "observation: " << observ << endl;
	for(int i=0 ; i<observ ; i++){
		for(int j=0 ; j<state ; j++){
			if(j!= state-1)
				modelfile << hmm.observation[i][j]<<" ";
			else
				modelfile << hmm.observation[i][j]<<endl;
		}
	}




	return 0;
}