//
// Definition of the alleles available with MixMHC2pred
//

#ifndef alleles_def_hpp
#define alleles_def_hpp

#include <vector>
#include <set>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <climits>
using namespace std;

string allele_def_file(string cAllele, const vector<string> &allelesDefFolders);
/* It'll return the full path to the file where the allele is defined searching
    for the allele in the different allelesDefFolders (using first match when
    the allele is defined in multiple folders). It returns "" when the allele
    isn't found in these folders. */

void get_allele_def(string cAllele, const vector<string> &allelesDefFolders,
    const map<char, int> &aa_to_ind, int &c_n_specif,
    vector<vector<vector<double>>> &c_PWM, vector<double> &c_w_k,
    vector<int> &c_spec_id, map<int, double> &c_alpha,
    vector<map<int, map<int, double>>> &c_w_ks,
    int &pepL_min_def, int &pepL_max_def, int &CS_upper_def, int &L_core);

void get_nnets_def(string nnetsDefFile, vector<double> &PWM_range,
    vector<int> &pepL_range, vector<int> &P1_range, double &log_min,
    vector<string> &nn_in_var, vector<vector<vector<double>>> &nnetsWeights,
    vector<vector<double>> &nnetsBiases, vector<double> &rank_thr,
    map<int, vector<double>> &sc_thr, 
    double &pRank_T, map<int, double> &f_L, bool noContext);

#endif /* alleles_def_hpp */