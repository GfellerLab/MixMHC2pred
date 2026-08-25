//
// Definition of the alleles available with MixMHC2pred
//

#include "alleles_def.hpp"
#include "helper_MixMHC2pred.hpp"


map<string, string> alleles_def_map(const vector<string> &allelesDefFolders,
    const vector<string> &alleles) {
    map<string, string> cmap;
    string cFolder;

    /* Internal function to check if the allele is part of the requested
        input allele, check if it is already present in cmap,
        and add it if needed. */
    auto check_and_insert = [&cmap, &alleles](string alleleName, string alleleMap,
        const string& cFolder) {
        alleleMap = cFolder + alleleMap + ".txt";
        for (int i(0); i <2; i++){
            if (i == 1){
                /* We determine if first letters of alleleName are 'HLA-' to
                    also allow having the names without HLA- part (i.e., to
                    allow giving the allele as HLA-DRB1*01:01, DRB1*01:01 or
                    DRB1_01_01). */
                if (alleleName.substr(0, 4) == "HLA-"){
                    alleleName = (alleleName.substr(4));
                } else {
                    break;
                    /* alleleName didn't start with 'HLA-', no need to reinsert
                        the same value as 2nd time. */
                }
            }
            
            if (find(alleles.begin(), alleles.end(), alleleName) != alleles.end()){
                /* Check if the given alleleName was part of the user-input alleles. */
                if ((cmap.find(alleleName) != cmap.end()) &&
                (cmap[alleleName] != alleleMap)) {
                    /* 2nd check is in case current alleleName was defined multiple
                    times pointing to the same file, we don't want to throw an error.
                    (This could be the case if the 2 allele columns were identical for
                    example for user-defined PWM files). */
                    throw(string("The allele '") + alleleName +
                    "' is defined multiple times in the alleles mapping files.");
                }
                cmap[alleleName] = alleleMap;
            }
        }        
    };

    /* We loop over the folders, reading the Alleles_map files found in these
        to determine the mapping between allele names. */
    for (int i(0); i < allelesDefFolders.size(); i++){
        cFolder = allelesDefFolders[i];
        vector<string> mapFiles = list_matching_files(cFolder, R"(Alleles_map.*\.txt)");
        // Search for files Alleles_map.....txt in the folder.
        if (mapFiles.size() == 0){
            throw(string("Couldn't find any Alleles_map file in the folder '") + cFolder + "'. Please make sure  to have at least one such file per folder containing the allele definitions files.");
        }

        for (const auto& cMapFile : mapFiles){
            ifstream allelesMapFile(cFolder + cMapFile);
            if (!allelesMapFile.is_open()){
                throw(string("Couldn't open the alleles mapping file '") +
                    cFolder + cMapFile + "'.");
            }
            string trash, alleleName1, alleleName2, alleleMap;
            while (allelesMapFile.peek() == '#'){
                safeGetline(allelesMapFile, trash);
                // Removing first comment lines
            }
            allelesMapFile >> alleleName1 >> alleleName2 >> alleleMap;
            /* We have two columns of allele names in the file, the first
                        one containing "official names" (i.e., IMGT-HLA format),
                        and the 2nd column with the name in my format with
                        underscores, etc. */
            if (alleleName1 != "OfficialName" || alleleName2 != "AlleleName" || alleleMap != "mapName"){
                throw(string("The alleles mapping file '") + cFolder + cMapFile +
                    "' doesn't have the expected header line (OfficialName, AlleleName, Allele_map, followed by additional columns).");
            }
            safeGetline(allelesMapFile, trash);

            while (allelesMapFile.good()){
                allelesMapFile >> alleleName1 >> alleleName2 >> alleleMap;
                safeGetline(allelesMapFile, trash);
                // Remove additional columns not used.

                if (!allelesMapFile.eof()){
                    /* We do this check here because if the last line was empty
                        the 'while' above still went in but it then didn't
                        get new values in alleleName1, etc.
                    */
                   check_and_insert(alleleName1, alleleMap, cFolder);
                   check_and_insert(alleleName2, alleleMap, cFolder);
                }
            }
        }
    }
    return cmap;
}


string allele_def_file(string cAllele, const map<string, string> &allelesMap, const vector<string> &allelesDefFolders) {

    map<string, string>::const_iterator it(allelesMap.find(cAllele));
    string cFile("");
    // Empty string used in case allele not found in the data.
    if (it != allelesMap.end()) {
        // Allele found in the map, return the corresponding file path.
        cFile = it->second;
        ifstream alleleFile(cFile);
        if (!alleleFile.is_open()) {
            cerr << "Warning: couldn't find the definition file '" << cFile << endl << " needed for allele '" << cAllele << "'." << endl;
            cFile = "";
            // Reset to empty string if the file cannot be opened (because file
            // didn't exist for example).
        }
    }
    return (cFile);
}

void get_allele_def(string cAllele, const map<string, string> &allelesMap,
    const vector<string> &allelesDefFolders,
    const map<char, int> &aa_to_ind, int &c_n_specif,
    vector<vector<vector<double>>> &c_PWM, vector<double> &c_w_k,
    vector<int> &c_spec_id, map<int, double> &c_alpha,
    vector<map<int, map<int, double>>> &c_w_ks,
    int &pepL_min_def, int &pepL_max_def, int &CS_upper_def, int &L_core) {
    
    vector< map<int, map<int, double>>> w_s;
    bool CS_perL;
    int x;
    string cAlleleFile(allele_def_file(cAllele, allelesMap, allelesDefFolders)),
        trash;
    ifstream alleleFile(cAlleleFile);
    char cAA;
    vector< int > x_vect;

    if (alleleFile.fail()){
        throw(string("Couldn't open allele definition: \n") + cAlleleFile);
    }
    alleleFile >> trash >> trash;
    /* Name of allele (the name is already in file name, we don't check it in
        case the user had renamed its alleles in user-defined allele files). */
    alleleFile >> trash >> trash; // Version number
    if ((!pars::ignore_version) && (trash != pars::version)){
        throw(string("Running MixMHC2pred v" + pars::version +
            " but allele definition file of " + cAllele + " (" +
            cAlleleFile + ") is version " + trash +
            ".\nPlease download new allele definition files."));
    }
    /* There are then some lines that aren't used ('Type', also the allele
        sequence when PWM build on user sequence) */
    while (trash != "n_specificities") {
        alleleFile >> trash;
    }
    alleleFile  >> c_n_specif;
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
    // There are then some lines that aren't used (e.g. npep_t)
    while (trash != "spec_ids") {
        alleleFile >> trash;  
    }
    c_spec_id = vector<int>(c_n_specif);
    for (int j(0); j < c_n_specif; j++) {
        alleleFile >> c_spec_id[j];
    }
    alleleFile >> trash >> L_core;
    while (trash != "CS_perL") {
        alleleFile >> trash;  
    }
    alleleFile >> trash;
    CS_perL = (trash == "1");
    /* Next lines are not used (e.g. pepL_lower/upper obtained later, n_diff_aa
        assumed to be the same for all alleles).
    */
    while (trash != "PWM_norm.1") {
        alleleFile >> trash;  
    }

    // Get the PWM sub specificities
    c_PWM = vector<vector<vector<double>>>(
        c_n_specif,
        vector<vector<double>>(L_core, vector<double>(pars::n_aa, 0.0)));
    for (int j(0); j < c_n_specif; j++) {
        alleleFile >> trash;
        safeGetline(alleleFile, trash); // PWM 'header'
        for (int k(0); k < pars::n_aa; k++) {
            alleleFile >> cAA; // The name of the aa
            for (int l(0); l < L_core; l++) {
                alleleFile >> c_PWM[j][l].at(aa_to_ind.at(cAA));
                // In the def file, the columns correspond to positions and
                // rows to different aa.
            }
        }
    }

    // Get the w_s and alpha values
    w_s.assign(c_n_specif, map<int, map<int, double>>());
    if (CS_perL) { // CS depend on the peptide sizes
        alleleFile >> trash;
        safeGetline(alleleFile, trash); // just 'header' defining these.
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

vector<string> nnets_def_files(string &nnetsFolder, string nnetsDefRoot,
                               string execFolder) {
    string cStr;
    vector<string> nnetsDefFiles;
    ifstream nnetsDef;
    if (nnetsFolder == "default") {
        nnetsFolder = execFolder + "/nnetsDef/";
    } else if (nnetsFolder.substr(0, 5) == "exec:") {
        nnetsFolder = execFolder + "/" + nnetsFolder.substr(5) + "/";
    } else {
        nnetsFolder += "/";
        // Make sure path ends with a '/'.
    }
    for (size_t c_nnet(1); c_nnet <= pars::n_nnets; c_nnet++) {
        nnetsDefFiles.push_back(nnetsFolder + nnetsDefRoot + "_" +
                                to_string(c_nnet) + ".txt");
        nnetsDef.open(nnetsDefFiles.back());
        if (nnetsDef.fail()){
            throw(string("Couldn't open the nnets definition: \n") + nnetsDefFiles.back());
        }
        if (!pars::ignore_version){
            nnetsDef >> cStr;
            if (cStr == "Version"){
                nnetsDef >> cStr; // cStr corresponds to version name.
            } else {
                throw(string("Expecting 'Version' in the nnetsDefFile ") +
                      nnetsDefFiles.back() + " but read '" + cStr + "'.");
            }
            if (pars::version == "Undefined"){
                pars::version = cStr;
                // Version hadn't yet been defined.
            } else if (pars::version != cStr){
                // Version was already defined but doesn't match current file.
                throw(string("Running MixMHC2pred v" + pars::version +
                    " but nnets definition file '" +
                    nnetsDefFiles.back() + "' is version " + cStr +
                    ".\nPlease make sure all nnets definition files are " +
                    "updated simultaneously."));
            }
        }
        nnetsDef.close();
    }
    return (nnetsDefFiles);
}


void get_nnets_def(string nnetsDefFile, vector<double> &PWM_range,
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
    PWM_range.resize(2);
    pepL_range.resize(2);
    P1_range.resize(2);
    nnetsDef >> cStr;
    while (cStr != "PWM_range") {
        nnetsDef >> cStr;
        // Possibly remove some first lines defining the version of the nnetsDef.
    }
    nnetsDef >> PWM_range[0] >> PWM_range[1];
    // if (cStr != "PWM_range") {
    //     throw(string("Expecting 'PWM_range' in the nnetsDefFile ") +
    //           nnetsDefFile + " but read '" + cStr + "'.");
    // }
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
