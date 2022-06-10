//
// Definition of the alleles available with MixMHC2pred
//

#include "alleles_def.hpp"
#include "helper_MixMHC2pred.hpp"


bool allele_is_defined(string cAllele, string allelesDefFolder){
    ifstream alleleFile(allelesDefFolder + cAllele + ".txt");
    return (alleleFile.good());
}

void get_allele_def(string cAllele, string allelesDefFolder,
    const map<char, int> &aa_to_ind, int &c_n_specif,
    vector<vector<vector<double>>> &c_PPM, vector<double> &c_w_k,
    vector<int> &c_spec_id, map<int, double> &c_alpha,
    vector<map<int, map<int, double>>> &c_w_ks,
    int &pepL_min_def, int &pepL_max_def, int &CS_upper_def, int &L_core) {
    
    vector< map<int, map<int, double>>> w_s;
    bool CS_perL(false);
    int x;
    string cAlleleFile(allelesDefFolder + cAllele + ".txt"), trash;
    ifstream alleleFile(cAlleleFile);
    char cAA;
    vector< int > x_vect;

    if (alleleFile.fail()){
        throw(string("Couldn't open allele definition: \n") + cAlleleFile);
    }
    alleleFile >> trash >> trash; // Name of allele
    alleleFile >> trash >> trash;      // Type of predictor definition
    alleleFile >> trash >> c_n_specif;
    // Get w_k values:
    alleleFile >> trash; // npep_u text
    c_w_k = vector<double>(c_n_specif);
    for (int j(0); j < c_n_specif; j++) {
        alleleFile >> c_w_k[j];
    }
    double npep_subSpec_tot(accumulate(c_w_k.begin(), c_w_k.end(), 0.0));
    for (int j(0); j < c_n_specif; j++) {
        c_w_k[j] = c_w_k[j] / npep_subSpec_tot;
    }
    // There is first e.g. the line for npep_t that isn't used.
        while ((trash != "L_core") && (trash != "spec_ids")) {
            alleleFile >> trash;  
        }
    c_spec_id = vector<int>(c_n_specif);
    if (trash == "spec_ids"){
        for (int j(0); j < c_n_specif; j++) {
            alleleFile >> c_spec_id[j];
        }
    } else {
        cerr << "#######\n#######\nWarning, spec_ids wasn't found in the "
            << "PPM definition for " << cAllele
            << " will assume all sub-spec are forward binders." << endl;
        for (int j(0); j < c_n_specif; j++) {
            c_spec_id[j] = j+1;
        }
    }

    /* Next information are some information that we assume to be the
        same between all alleles and don't re-read everything for the moment
        as these definition files will likely change in a near future as
        needed for the pan-specific predictor.  */
    while (trash != "L_core") {
        alleleFile >> trash;  
    }
    alleleFile >> L_core;
    do {
        alleleFile >> trash;  
    } while((trash != "CS_perL") && (trash != "Core_shift_upper"));
    if (trash == "CS_perL"){
        /* The older definition files didn't include CS_perL header and
            didn't define the CS per peptide length. */
        alleleFile >> trash;
        CS_perL = (trash == "1");
    }
    do {
        alleleFile >> trash;
    } while (trash != "PWM_norm.1");

    // Get the PPM sub specificities
    c_PPM = vector<vector<vector<double>>>(
        c_n_specif,
        vector<vector<double>>(L_core, vector<double>(pars::n_aa, 0.0)));
    for (int j(0); j < c_n_specif; j++) {
        alleleFile >> trash;
        getline(alleleFile, trash); // PPM 'header'
        for (int k(0); k < pars::n_aa; k++) {
            alleleFile >> cAA; // The name of the aa
            for (int l(0); l < L_core; l++) {
                alleleFile >> c_PPM[j][l].at(aa_to_ind.at(cAA));
                // In the def file, the columns correspond to positions and
                // rows to different aa.
            }
        }
    }

    // Get the w_s and alpha values
    w_s.assign(c_n_specif, map<int, map<int, double>>());
    if (CS_perL) { // CS depend on the peptide sizes
        alleleFile >> trash;
        getline(alleleFile, trash); // just 'header' defining these.
        pepL_min_def = INT_MAX;
        pepL_max_def = 0;
        while (!alleleFile.eof()) {
            alleleFile >> x >> c_alpha[x];
            if (pepL_min_def > x)
                pepL_min_def = x;
            if (pepL_max_def < x)
                pepL_max_def = x;
            for (int j(0); j < c_n_specif; j++){
                for (int CS_max(ceil((x - L_core) / 2.0)), CS(-CS_max);
                    CS <= CS_max; CS++) {
                    alleleFile >> w_s[j][x][CS];
                }
            }
        }
        CS_upper_def = ceil((pepL_max_def - L_core) / 2.0);

    } else { // CS are the same for all peptide sizes.
        alleleFile >> trash >> trash; // Core_shift test and 1st offset used
        while (trash != "w_s.1"){
            x_vect.push_back(stoi(trash)); // Defining the different offsets used
            alleleFile >> trash;
        }
        CS_upper_def = *(max_element(x_vect.begin(), x_vect.end()));
        for (int j(0); j < c_n_specif; j++){
            for (auto CS : x_vect){
                alleleFile >> w_s[j][0][CS];
            }
            alleleFile >> trash; // Will contain w_s.2 or pepL
        }

        // And getting the alpha values
        if (trash != "pepL"){
            throw(string("Expecting pepL for allele ") + cAllele +
                " but it was " + trash + ".");
        }
        x_vect.clear();
        alleleFile >> trash; // First pepL used
        while (trash != "alpha") {
            x_vect.push_back(stoi(trash)); // Defining the different pepL used
            alleleFile >> trash;
        }
        for (auto pepL : x_vect){
            alleleFile >> c_alpha[pepL];
        }
        auto pepL_range = minmax_element(x_vect.begin(), x_vect.end());
        pepL_min_def = *pepL_range.first;
        pepL_max_def = *pepL_range.second;
    }

    /* Multiplying the w_s by w_k to get w_ks because we won't need the w_s
        anymore and we'd be doing this multiplication often. */
    c_w_ks = vector<map<int, map<int, double>>>(c_n_specif);
    for (int j(0); j < c_n_specif; j++) {
        for (auto c_a : c_alpha){
            if (CS_perL){
                x = c_a.first;
            } else {
                x = 0;
                /* When we don't use CS per L, the w_s above had all been
                    assigned to a value 0 in the map, while with CS per L,
                    the map needs the pepL as first element. */
            }
            for (auto c_ws : w_s[j][x]){
                c_w_ks[j][c_a.first][c_ws.first] = c_w_k[j] * c_ws.second;
            }
        }
    }
}

void get_nnets_def(string nnetsDefFile, vector<double> &PPM_range,
    vector<int> &pepL_range, vector<int> &P1_range, double &log_min,
    vector<string> &nn_in_var, vector<vector<vector<double>>> &nnetsWeights,
    vector<vector<double>> &nnetsBiases, vector<double> &rank_thr,
    map<int, vector< double > > &sc_thr,
    double &pRank_T, map<int, double> &f_L, bool noContext) {
    
    ifstream nnetsDef(nnetsDefFile);
    string cStr;
    int nVal_1, nVal_2, nLayers;

    if (nnetsDef.fail()){
        throw(string("Couldn't open the nnets definition: \n") + nnetsDefFile);
    }
    
    /* We'll read the nnets params from the given definition file. cStr are the
        "header" telling what is the next data in the definition file.
    */
    PPM_range.resize(2);
    pepL_range.resize(2);
    P1_range.resize(2);
    nnetsDef >> cStr >> PPM_range[0] >> PPM_range[1];
    if (cStr != "PWM_range") {
        throw(string("Expecting 'PWM_range' in the nnetsDefFile ") +
              nnetsDefFile + " but read '" + cStr + "'.");
    }
    nnetsDef >> cStr >> pepL_range[0] >> pepL_range[1]
        >> cStr >> P1_range[0] >> P1_range[1]
        >> cStr >> log_min >> cStr >> nVal_1;
    nn_in_var.resize(nVal_1);
    for (int i(0); i < nVal_1; i++){
        nnetsDef >> nn_in_var[i];
    }
    nnetsDef >> cStr >> nLayers;
    if (cStr != "N_layers") {
        throw(string("Expecting 'N_layers' in the nnetsDefFile ") +
              nnetsDefFile + " but read '" + cStr + "'.");
    }
    for (int i(0); i < nLayers; i++){
        nnetsDef >> cStr >> nVal_1 >> nVal_2;
        if (cStr != "nnetsWeights") {
            throw(string("Expecting 'nnetsWeights' in the nnetsDefFile ") +
                  nnetsDefFile + " but read '" + cStr + "'.");
        }
        nnetsWeights.push_back(vector<vector<double>>(nVal_1,
            vector<double>(nVal_2)));
        nnetsBiases.push_back(vector<double>(nVal_2));
        for (int j(0); j < nVal_1; j++){
            for (int k(0); k < nVal_2; k++){
                nnetsDef >> nnetsWeights[i][j][k];
            }
        }
        nnetsDef >> cStr;
        for (int k(0); k < nVal_2; k++){
            nnetsDef >> nnetsBiases[i][k];
        }
    }

    nnetsDef >> cStr >> nVal_1;
    if (cStr != "rank_thr") {
        throw(string("Expecting 'rank_thr' in the nnetsDefFile ") + nnetsDefFile +
              " but read '" + cStr + "'.");
    }
    rank_thr.resize(nVal_1);
    for (int j(0); j < nVal_1; j++) {
        nnetsDef >> rank_thr[j];
    }
    nnetsDef >> cStr;
    if (cStr != "sc_thr") {
        throw(string("Expecting 'sc_thr' in the nnetsDefFile ") + nnetsDefFile +
              " but read '" + cStr + "'.");
    }
    for (int cL(pepL_range[0]); cL <= pepL_range[1]; cL++){
        nnetsDef >> nVal_2;
        if (nVal_2 != cL) {
            throw(string("Expecting a pepL '") + to_string(cL) +
                "' in the nnetsDefFile " + nnetsDefFile +
                " but read '" + to_string(nVal_2) + "'.");
        }
        sc_thr[cL].resize(nVal_1);
        int k(0);
        for (int j(0); j < 2*nVal_1; j++) {
            if (noContext != (j < nVal_1)){
                /* The first nVal_1 values correspond to sc_thresholds including
                    the context, while the next nVal_1 values correspond to
                    sc_thresholds when context was replaced by X. Here we only
                    need to keep the sc_thresholds corresponding to currently
                    asked case. */
                nnetsDef >> sc_thr[cL][k];
                k++;
            } else {
                nnetsDef >> cStr;
            }
        }
    }
    nnetsDef >> cStr >> pRank_T >> cStr;
    if (cStr != "f_L") {
        throw(string("Expecting 'f_L' in the nnetsDefFile ") + nnetsDefFile +
              " but read '" + cStr + "'.");
    }
    for (int cL(pepL_range[0]); cL <= pepL_range[1]; cL++){
        nnetsDef >> f_L[cL];
    }
}
