#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

#define H 3
#define D_STD 1
#define MU_STD 1
#define RHO_STD 1
#define SIGMA_STD 1


using namespace std;

struct params{
    int id;
    double d;
    double mu;
    double rho;
    double sigma;
    bool is_generated;
};

params generateStateProposal(params prev_params, int new_state_id){

    params state_proposal;

    int seedval = chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator;
    generator.seed(seedval);
    
    normal_distribution<double> d_distr(prev_params.d, D_STD);
    normal_distribution<double> mu_distr(prev_params.mu, MU_STD);
    normal_distribution<double> rho_distr(prev_params.rho, RHO_STD);
    normal_distribution<double> sigma_distr(prev_params.sigma, SIGMA_STD);

    state_proposal.id = new_state_id;
    state_proposal.d = d_distr(generator);
    state_proposal.mu = mu_distr(generator);
    state_proposal.rho = rho_distr(generator);
    state_proposal.sigma = sigma_distr(generator);
    state_proposal.is_generated = true;

    cout << "Proposal for state " << state_proposal.id << " starting from state " << prev_params.id << endl;
    cout << state_proposal.d << " " << prev_params.d << endl;
    //cout << state_proposal.mu << " " << prev_params.mu << endl;
    //cout << state_proposal.rho << " " << prev_params.rho << endl;
    //cout << state_proposal.sigma << " " << prev_params.sigma << endl;
    cout << endl;

    return state_proposal;
}


int main()
{

    params X_0;

    int total_states = pow(2, H);
    params *generated_states = new params [total_states];
    for(int j = 0; j < total_states; j++){
        generated_states[j].is_generated = false;
        generated_states[j].id = j;
    }

    generated_states[0].d = 1;
    generated_states[0].mu = 2;
    generated_states[0].rho = 4;
    generated_states[0].sigma = 1;
    generated_states[0].is_generated = true;

    for(int j = 0; j < total_states; j++)
    {

        if(j % 2 == 1)
            continue;

        for(int k = 1; j+k < total_states; k*=2)
        {
            if(generated_states[j+k].is_generated == true)
                break;
            else
            {
                generated_states[j+k] = generateStateProposal(generated_states[j], j+k);
            }
        }

    }

}
