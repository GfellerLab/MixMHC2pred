//
// Definition of the alleles available with MixMHC2pred
//

#ifndef alleles_def_hpp
#define alleles_def_hpp

#include <vector>
#include <set>
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <climits>
#include "helper_general.hpp"
using namespace std;


map<string, string> alleles_def_map(const vector<string> &allelesDefFolders,
    const vector<string> &alleles);
/* It'll return a map between allele names and the full path to the file to use
    for the definition of this allele (could correspond to the same allele name
    or the name of another allele with identical binding site sequence);
    Note that it'll be an error if a same allele name is defined multiple times
    (which could happen e.g. in different allelesDefFolders or if same name is
    present multiple times in the alleles_map file).
*/

string allele_def_file(string cAllele, const map<string, string> &allelesMap, const vector<string> &allelesDefFolders);
/* It'll return the full path to the file where the given allele is defined.
    It returns "" when the allele isn't found in these folders. */

void get_allele_def(string cAllele, const map<string, string> &allelesMap,
    const vector<string> &allelesDefFolders,
    const map<char, int> &aa_to_ind, int &c_n_specif,
    vector<vector<vector<double>>> &c_PWM, vector<double> &c_w_k,
    vector<int> &c_spec_id, map<int, double> &c_alpha,
    vector<map<int, map<int, double>>> &c_w_ks,
    int &pepL_min_def, int &pepL_max_def, int &CS_upper_def, int &L_core);

vector<string> nnets_def_files(string &nnetsFolder, string nnetsDefRoot,
    string execFolder);
/* To return the full path where the nnets definition files are defined. It'll
    also check if the version is defined in these files (defining pars::version
    with it) and verifying that all files have the same version. */

void get_nnets_def(string nnetsDefFile, vector<double> &PWM_range,
    vector<int> &pepL_range, vector<int> &P1_range, double &log_min,
    vector<string> &nn_in_var, vector<vector<vector<double>>> &nnetsWeights,
    vector<vector<double>> &nnetsBiases, vector<double> &rank_thr,
    map<int, vector<double>> &sc_thr, 
    double &pRank_T, map<int, double> &f_L, bool noContext);

#endif /* alleles_def_hpp */