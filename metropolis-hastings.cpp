#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <stdlib>
#include <time>


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
    double potential;
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

    srand(time(NULL));

    int total_states = pow(2, H);
    params *proposals = new params [total_states];
    for(int j = 0; j < total_states; j++){
        proposals[j].is_generated = false;
        proposals[j].id = j;
    }

    //Define starting state
    proposals[0].d = 0.5;
    proposals[0].mu = 0;
    proposals[0].rho = 0;
    proposals[0].sigma = 1;
    proposals[0].is_generated = true;

    //Generate state proposals for the pre-fetching tree
    for(int j = 0; j < total_states; j++)
    {

        if(j % 2 == 1)
            continue;

        for(int k = 1; j+k < total_states; k*=2)
        {
            if(proposals[j+k].is_generated == true)
                break;
            else
            {
                proposals[j+k] = generateStateProposal(generated_states[j], j+k);
            }
        }
    }

    //Calculate state potentials for finding acceptance probabilities
    for(int j = 0; j < total_states; j++)
    {
        
    }

    params* selected_states = new params [H];

    //Select states according to acceptance probability
    int same_selections = 0;
    double accept_prob;
    selected_states[0] = proposals[0];

    for(int j = 1; j < H; j++){
        int next_prop = selected_states[j-1].id + pow(2, same_selections);
        accept_prob = proposals[next_prop].potential/selected_states[j-1].potential;
        if(accept_prob >= 1)
            selected_states[j] = proposals[next_prop];
        else{
            if(rand() <= accept_prob)
                selected_states[j] = proposals[next_prop];
            else{
                selected_states[j] = selected_states[j-1];
                same_selections++;
            }
        }
    }

}
