//
//  helper_nnets.hpp
//  Declaration file of the main functions / classes used for the nnets part of
//  MixMHC2pred.
//
//  Created by Julien Racle on 09.07.21.
//  Copyright © 2021 CCB. All rights reserved.
//

#ifndef helper_nnets_hpp
#define helper_nnets_hpp

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <numeric>
#include <vector>
#include <string>
#include <climits>
#include "helper_general.hpp"
#include "helper_MixMHC2pred.hpp"

using namespace std;


/* We define a tempalte class valWithRange that will be used in the
 nn-encoding to tell the value of a given parameter as well as its allowed
 ranges. */
template< class T > class valWithRange {
public:
    valWithRange();
    valWithRange(string name, T min, T max, int n_elements=0, bool multi1=false);
    T min() const;
    T max() const;
    string get_name() const;
    bool empty() const;
    size_t size() const;
    void push_back(const T &val);
    T& operator[](int i);
    const T& operator[](int i) const;
    bool use_multi1() const;

private:
    vector< T > values;
    T val_min;
    T val_max;
    string name;
    bool multi1;
    // multi1 is used for the one-hot encoding cases only to give values of 1
    // in all columns "before the current one" and 0 in all columns after it
    // instead of having a single 1 per row.
};


// ################################
// ## Declaration of different functions used for the nnets.
// ################################
vector< vector< double > > encode_nn_in(size_t nData, const vector<string> &nn_in_var,
    const vector<string> &aaSeq,
    const valWithRange<double> &scores = valWithRange<double>(),
    const vector< valWithRange<int> > &one_hot=vector< valWithRange<int> >());
/* Function taking different parameters and putting these in a square matrix
    of dimension nData x nn_in_var.size() of the encoded values (aaSeq will be
    blosum-encoded, scores will be min-max normaliued and "one_hot" will be
    one-hot encoded). */

vector< vector< double > > nn_relu(const vector< vector<double> > &x,
    const vector< vector<double> > &w, const vector<double> &b);
/* Returns the rectified linear unit (ReLU) of the weights inputs and biases
    This function returns
        y_ki = max(0, sum_j(w_ji * x_kj) + b_i).
        (dim: n_obs * n_nodes)
    Where:
    - x A matrix of the input values (dim: n_obs * n_inputs).
    - w A matrix of the weights (dim: n_inputs * n_nodes).
    - b A vector of the biases (with n_nodes elements).
*/

vector< vector< double > > nn_logistic(const vector< vector<double> > &x,
    const vector< vector<double> > &w, const vector<double> &b);
/* Returns the logistic (sigmoid) function of the weights inputs and biases
    Ths function returns
        y_ki = 1 / (1 + exp(-(sum_j(w_ji * x_kj) + b_i))).
        (dim: n_obs * n_nodes).
    Where:
    - x A matrix of the input values (dim: n_obs * n_inputs).
    - w A matrix of the weights (dim: n_inputs * n_nodes).
    - b A vector of the biases (with n_nodes elements).
*/

vector< vector< double > > nn_softmax(const vector< vector<double> > &x,
    const vector< vector<double> > &w, const vector<double> &b);
/* Returns the softmax function of the weights inputs and biases
    This function first computes
    y_ki = exp(sum_j(w_ji * x_kj) + b_i)
    and it then returns
    z_ki = y_ki / (sum_j(y_kj)).
    (dim: n_obs * n_nodes).
    Where:
    - x A matrix of the input values (dim: n_obs * n_inputs).
    - w A matrix of the weights (dim: n_inputs * n_nodes).
    - b A vector of the biases (with n_nodes elements).
*/

void pred_neural_networks(string nnetsDefFile, vector<string> alleles,
    vector<Peptide> &Peptides, bool noContext);
/* Function doing the predictions of the 2nd neural networks.
*/

double transf_toPrank (double cScore, int pepL,
    const vector<double> &rank_thr,
    const map<int, vector<double>> &sc_thr,
    double pRank_T, const map<int, double> &f_L);
/* Returns the %Rank based on pre-computed scores-to-rank thresholds, making
    the follow the pepL distribution encoded by f_L.
    - scores The score of current peptide to transform to %Rank.
    - pepL Current peptide size.
    - rank_thr and sc_thr: The thresholds used for the transformation (rank_thr
        tells to which ranks the thresholds correspond and sc_thr tells for each
        pepL what are the threshold scores).
    - pRank_T the limit until which we want the %Ranks to follow the pepL
        distribution
    - f_L The pepL distribution to target in the transformation.
*/

// ################################
// ## Inline functions definitions
// ################################
template< class T > valWithRange<T>::valWithRange(){}
template< class T > valWithRange<T>::valWithRange(string c_name, T min, T max,
    int n_elements, bool c_multi1) :
    val_min(min), val_max(max), name(c_name), multi1(c_multi1) {
        if (n_elements > 0){
            values.resize(n_elements);
        }
    }
template< class T > inline T valWithRange<T>::min() const {
    return val_min;
}
template< class T > inline T valWithRange<T>::max() const {
    return val_max;
}
template< class T > inline string valWithRange<T>::get_name() const {
    return name;
}
template< class T > inline bool valWithRange<T>::empty() const {
    return values.empty();
}
template< class T > inline size_t valWithRange<T>::size() const {
    return values.size();
}
template< class T > inline void valWithRange<T>::push_back(const T &val){
    values.push_back(val);
}
template< class T > inline const T& valWithRange<T>::operator[](int i) const {
    return values[i];
}
template< class T > inline T& valWithRange<T>::operator[](int i) {
    return values[i];
}
template< class T > inline bool valWithRange<T>::use_multi1() const {
    return multi1;
}


#endif /* helper_nnets_hpp */
