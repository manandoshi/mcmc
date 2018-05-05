#include <iostream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

int getTimeSeries(int p, float d, int q, float variance, float *phi, float *nu){

	// initialize the series

	default_random_engine generator;
	//normal dist takes mean and std as an input
    normal_distribution<double> distribution(0,pow(variance, 0.5));

	float series[1000], errors[1000], diff_coeff[1000];
	int iter=0;

	// assign values to the first p terms. Typically, p is a single digit value

	for(int j=0; j<p; j++){

		series[j] = distribution(generator);
	}

	// populate the errors and diff coefficients

	for(int j=0; j<1000; j++){
		errors[j] = distribution(generator);
		diff_coeff[j] = (tgamma(d+1)/(tgamma(j+2)*tgamma(d-j)))*(pow(-1, j+1));
		cout<<tgamma(d+j)<<endl;
		//0th value corresponds to k=1 in the formula  
	}


	for(int i=p; i<1000; i++){

		// calculate MA terms for ith value

		float MA = errors[i];
		iter = fmin(i,q);
		for(int j=1; j<=iter; j++){
			MA += nu[j-1]*errors[i-j];
		}

		// Calculate AR and Differential Terms for ith value 

		 float firstterm=0, secondterm=0, thirdterm=0;

		 for(int j=1; j<=p; j++){
		 	firstterm += phi[j-1]*series[i-j];
		 	
		 }
		 for(int j=1; j<=i; j++){
		 	secondterm += diff_coeff[j-1]*series[i-j];
		 }
		 for(int j=1; j<=p; j++){
		 	for(int k=1; k<=i; k++){
		 		if(i-j-k >= 0){
		 			thirdterm += phi[j-1]*diff_coeff[j-1]*series[i-j-k];
		 		}
		 	}
		 }

		//assign the ith value 
		series[i] = firstterm + secondterm + thirdterm + MA;
	}
	//for(int i=0; i<1000; i++)
	//	cout << series[i];
	return 0;

}

int main(){
	float phi[3]={3,4,5}, nu[3]={1,2,3};
	getTimeSeries(3,3,3,4,phi, nu);
	cout<<phi[2]<<endl;
	return 0;
}
