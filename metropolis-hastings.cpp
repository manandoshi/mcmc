#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include "mle.h"

#define H 3
#define D_STD 0.02
#define PHI_STD 0.05
#define VAR_STD 0.8


using namespace std;

struct params{
    int id;
    float d;
    float phi;
    float var;
    bool is_generated;
    float negln_potential;
    bool is_valid;
};

struct mcmc_state{

    float d;
    float phi;
    float var;
};

mcmc_state copyToState(params p){

    mcmc_state new_state;
    new_state.d = p.d;
    new_state.phi = p.phi;
    new_state.var = p.var;

    return new_state;
}

params generateStateProposal(params prev_params, int new_state_id){

    params state_proposal;

    int seedval = chrono::system_clock::now().time_since_epoch().count();

    default_random_engine generator;
    generator.seed(seedval);
    
    normal_distribution<float> d_distr(prev_params.d, D_STD);
    normal_distribution<float> phi_distr(prev_params.phi, PHI_STD);
    normal_distribution<float> var_distr(prev_params.var, VAR_STD);

    //Get state parameters from normal distribution
    state_proposal.d = d_distr(generator);
    state_proposal.phi = phi_distr(generator);
    state_proposal.var = var_distr(generator);

    //Set other values for state proposal
    state_proposal.id = new_state_id;
    state_proposal.is_generated = true;

    //cout << "Proposal for state " << state_proposal.id << " starting from state " << prev_params.id << endl;
    //cout << state_proposal.d << " " << prev_params.d << endl;
    //cout << endl;

    return state_proposal;
}

void printState(mcmc_state state)
{

    fstream fs_d, fs_phi, fs_var;
    fs_d.open("d.txt", ios::app);
    fs_phi.open("phi.txt", ios::app);
    fs_var.open("var.txt", ios::app);

    cout << "d: " << state.d << "\tphi: " << state.phi << "\tvar: " << state.var << endl;

    fs_d << state.d << endl;
    fs_phi << state.phi << endl;
    fs_var << state.var << endl;

    fs_d.close();
    fs_phi.close();
    fs_var.close();
}

params* generateStateTree(mcmc_state start_state, double* y, int n)
{

    srand(time(NULL));

    int total_states = pow(2, H);
    params *proposals = new params [total_states];
    for(int j = 0; j < total_states; j++){
        proposals[j].is_generated = false;
        proposals[j].id = j;
    }

    //Define starting state
    proposals[0].d = start_state.d;
    proposals[0].phi = start_state.phi;
    proposals[0].var = start_state.var;
    proposals[0].id = 0;
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
                proposals[j+k] = generateStateProposal(proposals[j], j+k);
            }
        }
    }
    omp_set_num_threads(total_states);

    //Calculate state potentials for finding acceptance probabilities
    #pragma omp parallel for
    for(int j = 0; j < total_states; j++)
    {
        //cout << "up top\n";
        if(proposals[j].d <= -0.5 || proposals[j].d >= 0.5 || proposals[j].var <= 0)
        {
            proposals[j].negln_potential = 0;
            proposals[j].is_valid = false;
        }
        else
        {
            proposals[j].negln_potential = calc_MLE(y, n, proposals[j].var, proposals[j].d, proposals[j].phi);
            proposals[j].is_valid = true;
        }
            
    }
    params* selected_states = new params [H+1];

    //Select states according to acceptance probability
    int same_selections = H-1;
    float accept_prob;
    selected_states[0] = proposals[0];

    for(int j = 0; j < H; j++)
    {
        //Get next proposal index
        int next_prop = selected_states[j].id + pow(2, H-1-j);
        //Get acceptance probability
        if(proposals[next_prop].is_valid == true){
            accept_prob = exp(-proposals[next_prop].negln_potential + selected_states[j].negln_potential);
        }else{
            accept_prob = 0;
        }
        //Select the j^th state as either new proposal or same as previous state
        if(accept_prob >= 1)
            selected_states[j+1] = proposals[next_prop];
        else{
            if((float) rand() / RAND_MAX <= accept_prob)
                selected_states[j+1] = proposals[next_prop];
            else{
                selected_states[j+1] = selected_states[j];
                same_selections++;
            }
        }
    }
    delete proposals;

    return selected_states;

}

int main(int argc, char **argv)
{
    string fname = argv[1];
    int size = stoi(argv[2]);
    int n = stoi(argv[3]);
    mcmc_state state[n];

    //Getting the data
    ifstream ifile(fname);
    string value;
    double y[n];
    int i = 0;
    while(getline(ifile, value))
    {
        y[i] = stod(value);
        i++;
    }


    state[0].d = 0.1;
    state[0].phi = 0.1;
    state[0].var = 0.1;

    params* next_state;
    for(int new_start_state = 0; new_start_state <= n; new_start_state += H)
    {
        if(new_start_state % 60 == 0)
            cout << "Finished state " << new_start_state << endl;
        next_state = generateStateTree(state[new_start_state], y, size);
        for(int j = 1; j <= H; j++){
            state[new_start_state+j] = copyToState(next_state[j]);
            cout<<new_start_state + j<<" "<<next_state[j].negln_potential<<endl;
        }
        delete next_state;
    }

    for(int j = 0; j < n; j++)
    {
        printState(state[j]);
    }
}
