//
//  helper_MixMHC2pred.hpp
//  Definition file of the main functions / classes used for MixMHC2pred.
//
//  Created by Julien Racle on 30.01.19.
//  Copyright © 2019 CCB. All rights reserved.
//

#include "helper_MixMHC2pred.hpp"
#include "alleles_def.hpp"

// #########################
// ## Definition of the class methods and variables
// #########################
// ***** Class of pars containing parameter values to use as global variables.
/* Definintion of the default parameter values. */
vector<string> pars::alleles;
string pars::version;
bool pars::ignore_version;
size_t pars::context_L = 3, pars::n_aa = 20, pars::L_core = 9, pars::n_nnets = 5;
map<char, int> pars::aa_to_ind = { {'A', 0}, {'C', 1}, {'D', 2}, {'E', 3}, {'F', 4}, {'G', 5}, {'H', 6}, {'I', 7}, {'K', 8}, {'L', 9}, {'M', 10}, {'N', 11}, {'P', 12}, {'Q', 13}, {'R', 14}, {'S', 15}, {'T', 16}, {'V', 17}, {'W', 18}, {'Y', 19}, {'X', 20}, {'-', 20} };

vector<char> pars::blosumColumns = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-'};
map< char, vector<double> > pars::blo_plus_1 = { 
	{'A', vector<double>{ 1.2902, 0.0216, 0.0297, 0.0405, 0.0216, 0.0783, 0.0148, 0.0432, 0.0445, 0.0594, 0.0175, 0.0256, 0.0297, 0.0256, 0.031, 0.085, 0.0499, 0.0688, 0.0054, 0.0175, 0} },
	{'C', vector<double>{ 0.065, 1.4837, 0.0163, 0.0163, 0.0203, 0.0325, 0.0081, 0.0447, 0.0203, 0.065, 0.0163, 0.0163, 0.0163, 0.0122, 0.0163, 0.0407, 0.0366, 0.0569, 0.0041, 0.0122, 0} },
	{'D', vector<double>{ 0.041, 0.0075, 1.3974, 0.0914, 0.0149, 0.0466, 0.0187, 0.0224, 0.0448, 0.028, 0.0093, 0.069, 0.0224, 0.0299, 0.0299, 0.0522, 0.0354, 0.0243, 0.0037, 0.0112, 0} },
	{'E', vector<double>{ 0.0552, 0.0074, 0.0902, 1.2965, 0.0166, 0.035, 0.0258, 0.0221, 0.0755, 0.0368, 0.0129, 0.0405, 0.0258, 0.0645, 0.0497, 0.0552, 0.0368, 0.0313, 0.0055, 0.0166, 0} },
	{'F', vector<double>{ 0.0338, 0.0106, 0.0169, 0.019, 1.3869, 0.0254, 0.0169, 0.0634, 0.019, 0.1142, 0.0254, 0.0169, 0.0106, 0.0106, 0.019, 0.0254, 0.0254, 0.055, 0.0169, 0.0888, 0} },
	{'G', vector<double>{ 0.0783, 0.0108, 0.0337, 0.0256, 0.0162, 1.5102, 0.0135, 0.0189, 0.0337, 0.0283, 0.0094, 0.0391, 0.0189, 0.0189, 0.0229, 0.0513, 0.0297, 0.0243, 0.0054, 0.0108, 0} },
	{'H', vector<double>{ 0.042, 0.0076, 0.0382, 0.0534, 0.0305, 0.0382, 1.355, 0.0229, 0.0458, 0.0382, 0.0153, 0.0534, 0.0191, 0.0382, 0.0458, 0.042, 0.0267, 0.0229, 0.0076, 0.0573, 0} },
	{'I', vector<double>{ 0.0471, 0.0162, 0.0177, 0.0177, 0.0442, 0.0206, 0.0088, 1.271, 0.0236, 0.1679, 0.0368, 0.0147, 0.0147, 0.0133, 0.0177, 0.025, 0.0398, 0.1767, 0.0059, 0.0206, 0} },
	{'K', vector<double>{ 0.057, 0.0086, 0.0415, 0.0708, 0.0155, 0.0432, 0.0207, 0.0276, 1.2781, 0.0432, 0.0155, 0.0415, 0.0276, 0.0535, 0.1071, 0.0535, 0.0397, 0.0328, 0.0052, 0.0173, 0} },
	{'L', vector<double>{ 0.0445, 0.0162, 0.0152, 0.0202, 0.0547, 0.0213, 0.0101, 0.1154, 0.0253, 1.3754, 0.0496, 0.0142, 0.0142, 0.0162, 0.0243, 0.0243, 0.0334, 0.0962, 0.0071, 0.0223, 0} },
	{'M', vector<double>{ 0.0522, 0.0161, 0.0201, 0.0281, 0.0482, 0.0281, 0.0161, 0.1004, 0.0361, 0.1968, 1.1606, 0.0201, 0.0161, 0.0281, 0.0321, 0.0361, 0.0402, 0.0924, 0.008, 0.0241, 0} },
	{'N', vector<double>{ 0.0427, 0.009, 0.0831, 0.0494, 0.018, 0.0652, 0.0315, 0.0225, 0.0539, 0.0315, 0.0112, 1.3169, 0.0202, 0.0337, 0.0449, 0.0697, 0.0494, 0.027, 0.0045, 0.0157, 0} },
	{'P', vector<double>{ 0.0568, 0.0103, 0.031, 0.0362, 0.0129, 0.0362, 0.0129, 0.0258, 0.0413, 0.0362, 0.0103, 0.0233, 1.4936, 0.0207, 0.0258, 0.0439, 0.0362, 0.031, 0.0026, 0.0129, 0} },
	{'Q', vector<double>{ 0.0559, 0.0088, 0.0471, 0.1029, 0.0147, 0.0412, 0.0294, 0.0265, 0.0912, 0.0471, 0.0206, 0.0441, 0.0235, 1.2147, 0.0735, 0.0559, 0.0412, 0.0353, 0.0059, 0.0206, 0} },
	{'R', vector<double>{ 0.0446, 0.0078, 0.031, 0.0523, 0.0174, 0.0329, 0.0233, 0.0233, 0.1202, 0.0465, 0.0155, 0.0388, 0.0194, 0.0484, 1.345, 0.0446, 0.0349, 0.031, 0.0058, 0.0174, 0} },
	{'S', vector<double>{ 0.1099, 0.0175, 0.0489, 0.0524, 0.0209, 0.0663, 0.0192, 0.0297, 0.0541, 0.0419, 0.0157, 0.0541, 0.0297, 0.0332, 0.0401, 1.2199, 0.082, 0.0419, 0.0052, 0.0175, 0} },
	{'T', vector<double>{ 0.073, 0.0178, 0.0375, 0.0394, 0.0237, 0.0434, 0.0138, 0.0533, 0.0454, 0.0651, 0.0197, 0.0434, 0.0276, 0.0276, 0.0355, 0.0927, 1.2465, 0.071, 0.0059, 0.0178, 0} },
	{'V', vector<double>{ 0.07, 0.0192, 0.0178, 0.0233, 0.0357, 0.0247, 0.0082, 0.1646, 0.0261, 0.1303, 0.0316, 0.0165, 0.0165, 0.0165, 0.0219, 0.0329, 0.0494, 1.2688, 0.0055, 0.0206, 0} },
	{'W', vector<double>{ 0.0303, 0.0076, 0.0152, 0.0227, 0.0606, 0.0303, 0.0152, 0.0303, 0.0227, 0.053, 0.0152, 0.0152, 0.0076, 0.0152, 0.0227, 0.0227, 0.0227, 0.0303, 1.4924, 0.0682, 0} },
	{'Y', vector<double>{ 0.0405, 0.0093, 0.0187, 0.028, 0.1308, 0.0249, 0.0467, 0.0436, 0.0312, 0.0685, 0.0187, 0.0218, 0.0156, 0.0218, 0.028, 0.0312, 0.028, 0.0467, 0.028, 1.3179, 0} },
	{'-', vector<double>{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2} },
	{'X', vector<double>{ 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476, 0.0476} } };

// ***** Class of Peptide
Peptide::Peptide(int nAlleles, string seq, string context, bool innerContext)
    : sequence(seq), score_isNA(false) {
        PWMrank_L = PWMrawScore = vector<double>(nAlleles, nan(""));
        nnets_rank = nnets_rawScore = vector<double>(nAlleles, 0.0);
        /* Need to set these to 0 for nnets instead of Nan as we'll add together
           the values from the various repetitions (without checking if current
           repet is the 1st one or not). Instead, we have a boolean indicating
           if this nnets_rank/rawScore is NA. */
        P1 = vector<int>(nAlleles, -1);
        best_subSpec_ids = vector<int>(nAlleles, 0);
        if (context.empty()){
            if (innerContext){
                context = seq.substr(0, pars::context_L)
                    + seq.substr(seq.size() - pars::context_L, pars::context_L);
                // Only inner context, based on first and last AAs of the peptide.
            } else {
                context.assign(4*pars::context_L, 'X');
                // Used the implementation without context corresponding to
                // having 'X' for all context AAs.
            }
        } else if (context.size() != 4*pars::context_L){
                throw(string("The context (" + context + ") of peptide " + seq
                    + " doesn't have the expected size of "
                    + to_string(4*pars::context_L) + " (please make sure that "
                    + "your input file contains two columns, with the first "
                    + "column corresponding to the peptide sequence and the "
                    + "second column giving the peptide's context sequence)."));
        }
        pepCon = context;
    }

int Peptide::get_P1_centered(int Aind) const {
    int cP1(P1[Aind]);
    cP1 = cP1 - ceil(double(this->get_pepL() - pars::L_core) / 2.0);
    if ((cP1 >= 0) && (((this->get_pepL() - pars::L_core) % 2) != 0))
        cP1++;
        // Check if the value 0 is absent for current peptide based on its size
    return cP1;
}

// ################################
// ## Helper functions definitions
// ################################
void pepData_import(string pepFile, vector<Peptide> &peptides,
    set<int> &pepL_used, string contextType, int nAlleles) {
    ifstream inFile(pepFile);
    string cLine, cPep(""), cContext("");
    bool fasta(false);

    if (inFile.fail())
        throw(string("Couldn't open the input file ") + pepFile);


    while (inFile.good()) {
        if (contextType == "full"){
            inFile >> cPep >> cContext;
            if (!inFile.fail()){
                if (cPep != "") {
                        peptides.push_back(Peptide(nAlleles, cPep, cContext));
                        pepL_used.insert(peptides.back().get_pepL());
                    }
            }
        } else {
            safeGetline(inFile, cLine);
            // We read the file with the peptides, line by line.
            if ((cLine[0] == '>') && (!fasta)) {
                if (peptides.empty()) {
                    fasta = true;
                } else {
                    throw(string("Fasta files should start with '>'. If ") +
                          "this is not a fasta file, it should not contain any "
                          "'>' " +
                          "character.");
                }
            }
            if (inFile.eof() && fasta && (cLine[0] != '>')) {
                // Reached last line of file in fasta format; we need to tell
                // this is also the end of the current peptide.
                cPep = cPep + cLine;
                cLine = "> last";
                inFile.clear(ios::eofbit);
                /* Put the eofbit only instead of possibly also the failbit
                    that was present if last line was empty to make this peptide
                    is indeed included in the list. */
            }
            if (!inFile.fail()) {
                if (cLine[0] != '>')
                    cPep += cLine;
                /* For fasta files, when the line starts with ">" it is the
                    name/id of the next peptide and so won't add this line to
                   the peptide sequence, we'll use it to tell below that we
                    reached the end of a peptide (in fasta files, all lines
                    between two lines with '>' correspond to the aa sequence of
                    a same peptide/protein, so need to first read all those
                    lines before saving the peptide sequence). Note that we
                   don't use the name of peptide or other characteristics that
                   are given after the '>' symbol from fasta files.
                */
                if ((!fasta) || (cLine[0] == '>')) {
                    /* We have a full peptide sequence (after each line when not
                        in fasta format or when current line was describing a
                       new peptide. */
                    if (cPep != "") {
                        peptides.push_back(Peptide(nAlleles, cPep, "", 
                            (contextType == "inner")));
                        pepL_used.insert(peptides.back().get_pepL());
                    }
                    cPep = ""; // Put back to empty sequence for next peptide
                }
            }
        }
    }
    inFile.close();
}

void comp_pepScores(vector<Peptide> &peptides,
                    int ws_scoring, bool needAllAlleles,
                    string rPepFolder, set<int> &pepL_used,
                    const vector<string> &allelesDefFolders) {
    string tString;
    vector<string> alleles(pars::alleles), allelesMissing, predDefFiles;
    int pepL, pepL_min_def, npep(peptides.size()), nAlleles,
        CS_up_cpep, cP1, nRandPep(0), nRandPep_perL, L_core, CS_upper,
        pepL_max_def, c_L_core, c_CS_upper_def, c_pepL_min_def, c_pepL_max_def;
    /* Peptides shorter than pepL_min_def will have a score of NA.
        nRandPep is the total number of random
       peptides used (will be obtained below), while nRandPep_perL is the
       approximate number of random peptides we'll use per pepL to compute
       score of the current peptides relative to random ones (for longer
       peptides we use less random peptides as anyway there are only few
       binders, so we'll be a bit less precise for those, in order to not
       explode the size of our package with random sequences) - this value is
       used to reserve space needed for the random peptide sequences.
    */
    bool with_w_s_0;

    // Definition of some parameters describing the system
    nRandPep_perL = 10000;

    // Check first if there are some missing alleles
    for (unsigned int i(0); i < alleles.size();) {
        if (allele_def_file(alleles[i], allelesDefFolders) != "") {
            i++;
        } else {
            allelesMissing.push_back(alleles[i]);
            alleles.erase(alleles.begin() + i);
        }
    }

    if (allelesMissing.size() > 0) {
        tString = allelesMissing[0];
        for (unsigned int i(1); i < allelesMissing.size(); i++)
            tString += ", " + allelesMissing[i];
        string alleles_def_folers_str("");
        for (string cFolder : allelesDefFolders) {
            alleles_def_folers_str += "'" + cFolder + "', ";
        }
        if (needAllAlleles || (alleles.size() == 0)){
            throw(string("The allele(s) '") + tString + "' aren't defined in " +
                alleles_def_folers_str + " cannot perform the computations.");
        } else {
            cerr << "Warning: The allele(s) '" << tString
                << "' aren't defined in " <<  alleles_def_folers_str
                << " doing the computation on the remaining alleles." << endl;
        }
    }
    nAlleles = alleles.size();
    if (nAlleles == 0) {
        throw(string("Need at least one allele to do the predictions."));
    }

    // Load the parameters of the alleles asked for the predictions.
    vector<int> n_specif(nAlleles);
    vector<vector<vector<vector<double>>>> PWM(nAlleles);
    vector<map<int, double>> alpha(nAlleles);
    vector<vector<double>> w_k(nAlleles);
    vector<vector<int>> spec_ids(nAlleles);
    vector<vector<map<int, map<int, double>>>> w_ks(nAlleles);
    /* w_ks has 4 dimensions: (nAlleles, n_spec, n_pepL, n_CS) with the 2 last
        dimensions being maps to reach them based on the pepL and CS values
        directly.
    */
    for (int i(0); i < nAlleles; i++) {
        get_allele_def(alleles[i], allelesDefFolders, pars::aa_to_ind,
            n_specif[i], PWM[i], w_k[i], spec_ids[i], alpha[i], w_ks[i],
            c_pepL_min_def, c_pepL_max_def, c_CS_upper_def, c_L_core);
        
        if (i == 0){
            L_core = c_L_core;
            pepL_min_def = c_pepL_min_def;
            pepL_max_def = c_pepL_max_def;
            CS_upper = c_CS_upper_def;
        } else {
            if ((L_core != c_L_core) ||
                (pepL_min_def != c_pepL_min_def) ||
                (pepL_max_def != c_pepL_max_def) || (CS_upper != c_CS_upper_def)){
                throw(string("It seems that not all alleles have the same ") +
                    "L_core, pepL_min/max or CS_upper.");
            }
        }
        // Adding the unspec AA to the PWMs (i.e. we consider an average PWM
        // for these cases)
        for (int j(0); j < n_specif[i]; j++) {
            for (int l(0); l < L_core; l++) {
                PWM[i][j][l].push_back(
                    accumulate(PWM[i][j][l].begin(), PWM[i][j][l].end(), 0.0));
                PWM[i][j][l].back() /= double(pars::n_aa);
            }
        }
    }

    // Computing the scores from each peptide.
    vector<int> cPep, bestP1_perAllele(nAlleles), best_subspec_perAllele(nAlleles);
    vector<double> score_perAllele(nAlleles), prctileScore_perL_perA(nAlleles);
    string cPepSeq, coreSeq;
    double subScore, sScore, prctileScore_perL, bestAlignScore_main,
        bestAlignScore_sub, ca, cb;
    int bestP1_main, bestP1_sub, Lmin, Lmax, nDiffL, prctileInd;
    

    // First, define some additional variables to load random peptides
    // and get their scores to be able to give percentile scores to the true
    // peptides based on the random peptides.
    ifstream rPepFile;
    auto pepL_range = minmax_element(pepL_used.begin(), pepL_used.end());
    Lmax = min(pepL_max_def, *pepL_range.second);
    // Won't add random peptides of size bigger than pepL_max_def nor bigger
    // than the longest peptide from the input data.
    Lmin = max(max(pepL_min_def, L_core), *pepL_range.first);
    /* Peptides shorter than that will not be predicted (we compare our
        threshold to the core-length in case the core_length was bigger than
        our threshold). We also don't need to take shorter random peptides
        than the sizes present in the input peptides.
    */
    nDiffL = Lmax - Lmin + 1;
    if (nDiffL < 0){
        /* This happens when all input peptides are shorter than pepL_min_def
            (or L_core) or all peptides are longer than pepL_max_def. In this
            case we won't need to do the computations for
            random peptides because the score for such short peptides will
            simply be a "NA" value. Therefore we set nDiffL = 0 (making no
            random peptides are used) - otherwise, if nDiffL < 0 it makes that
            the code crashed when reading a negative index from a vector. */
        nDiffL = 0;
    }

    map<int, vector<vector<double>>> rPepScores;
    /* rPepScores is a matrix [nDiffL][nAlleles][nRandPep_perL] that will
        contain the scores from each random peptide towards each allele. */
    vector<double> *rPepScores_pt;
    /* rPepScores_pt is a pointer that will be used when comparing current pep
        score against random peptides so that we don't always need to call the
        first two dimensions of the rPepScores 'matrix', which can be a bit more
        efficient. */

    // Now we first load all the peptides and random peptides into a common
    //  vector (we don't do this anymore directly in the loop of the scores
    //  because we might need to have the full list first in case N/Cterm_score
    //  are computed).
    vector<Peptide> allPeps;
    allPeps.reserve(nDiffL * nRandPep_perL + npep);
    if (nRandPep_perL > 0) {
        int nRand_cL;
        for (int pepL(Lmin); pepL <= Lmax; pepL++) {
            rPepFile.open(rPepFolder + "rand_L_" + to_string(pepL) + ".txt");
            if (rPepFile.fail()) {
                throw(string("There was an error reading the random peptide file ")
                    + rPepFolder + "rand_L_" + to_string(pepL) + ".txt.");
            }
            nRand_cL = 0;
            while (rPepFile.good()){
                rPepFile >> cPepSeq;
                if (!rPepFile.fail()){
                    /* It fails when we reach the end of the file, trying to
                        get the sequence of an inexisting additional peptide. */
                    allPeps.push_back(Peptide(nAlleles, cPepSeq));
                    nRand_cL++;
                }
            }
            rPepFile.close();
            nRandPep += nRand_cL;
            rPepScores[pepL] = vector<vector<double>>(nAlleles,
                vector<double>(nRand_cL));
            /* Make that this rPepScores has good size to accomodate all
                peptides from this length. */
        }
    }
    allPeps.insert(allPeps.end(), peptides.begin(), peptides.end());

    // And finally loop over all peptides to compute their scores.
    size_t pRand(0), pAll(0);
    // Declares them outside of the for loop as not of same type than p.
    for (int p(-nRandPep); p < npep; p++, pAll++) {
        // Negative p will be for the random peptides for which we'll first
        // compute the scores and p >= 0 will be for the peptides from the
        // input.
        cPepSeq = allPeps[pAll].getSequence();
        pepL = cPepSeq.size();
        if ((p >= 0) && ((pepL < Lmin) || (pepL > Lmax))) {
            peptides[p].set_na_score();
            goto undef_pep_score;
            /* Peptide is too short or too long, it cannot fit in the alleles,
                it's score will remain NA (was set when building all the
                peptides), to indicate the score of this peptide cannot really
                be computed. */
        }
        cPep.assign(pepL, 0);
        for (int l(0); l < pepL; l++) {
            if (pars::aa_to_ind.count(cPepSeq[l]) == 0){
                if (p < 0){
                    throw(string("Random peptides cannot contain undefined ") +
                        "amino acids!");
                }
                peptides[p].set_na_score();
                goto undef_pep_score;
                /* Some AA not defined in the PWMs is present: this peptide
                    will have a NA score. */
            } else 
                cPep[l] = pars::aa_to_ind[cPepSeq[l]];
        }
        with_w_s_0 = ((pepL - L_core) % 2) == 0;
        /* Depending on peptide and core sizes, the position s=0 will be
            available or not, in order to make the allowed s values fully
            symmetric. */
        CS_up_cpep = ceil(double(pepL - L_core) / 2.0);

        for (int i(0); i < nAlleles; i++) {
            score_perAllele[i] = 0;
            prctileScore_perL_perA[i] = 0;
            bestAlignScore_main = -1;
            bestP1_main = -1;
            // These values are per allele to determine within each allele where
            // is the best peptide alignment.
            for (int j(0); j < n_specif[i]; j++) {
                subScore = 0;
                bestAlignScore_sub = 0;
                bestP1_sub = -1;
                // Similarly but within the subspecificities
                cP1 = 0;
                for (int s(-CS_up_cpep); s <= CS_up_cpep; s++) {
                    if ((s == 0) && (!with_w_s_0)) {
                        // This alignment value isn't possible based on peptide
                        // length (was made to be fully symmetric around 0).
                        continue;
                    }
                    if ((ws_scoring == 2) || (ws_scoring == 3)) {
                        sScore = w_k[i][j];
                        /* We don't give weight to the offset with these
                            ws_scoring schemes. Still need to multiply by the
                            "weight" of the given subspecificity. */
                    } else if (s <= -CS_upper) {
                        sScore = w_ks[i][j][pepL][-CS_upper] /
                                 (CS_up_cpep - CS_upper + 1);
                        // To share these extreme w_ks values among all extreme
                        // core shifts from N-term side.
                    } else if (s >= CS_upper) {
                        sScore = w_ks[i][j][pepL][CS_upper] /
                                 (CS_up_cpep - CS_upper + 1);
                        // To share these extreme w_ks values among all extreme
                        // core shifts from C-term side.
                    } else {
                        sScore = w_ks[i][j][pepL][s];
                    }
                    for (int l(0); l < L_core; l++) {
                        sScore *= PWM[i][j][l][cPep[l + cP1]];
                    }
                    subScore += sScore;
                    if (sScore > bestAlignScore_sub) {
                        bestP1_sub = cP1;
                        bestAlignScore_sub = sScore;
                    }
                    cP1++;
                }
                if ((ws_scoring == 1) || (ws_scoring == 2)){
                    subScore = bestAlignScore_sub;
                    // Count the best alignment instead of summing over all.
                }
                if (ws_scoring == 4){
                    /* We just keep the score from best alignment, multiplied
                        by the weight of the given subspecifity, without
                        weighting the alignment. Need thus to recompute the score
                    */
                   subScore = w_k[i][j];
                   for (int l(0); l < L_core; l++) {
                        subScore *= PWM[i][j][l][cPep[l + bestP1_sub]];
                    }
                }

                if (subScore > bestAlignScore_main) {
                    bestAlignScore_main = subScore;
                    bestP1_main = bestP1_sub;
                    best_subspec_perAllele[i] = spec_ids[i][j];
                    // To indicate which is the best sub-specificity and
                    // alignment, we compare the full scores of each
                    // sub-specificities (i.e. from all its alignments), and we
                    // consider the best alignment within this best sub-spec.
                }
                score_perAllele[i] += subScore;
            }
            bestP1_perAllele[i] = bestP1_main;

        }

        if (p >= 0) {
            // It is an input peptide, need to save it's score for the output
            // We first need to transform the raw scores to percentile scores to
            // have a comparable scale between the alleles.
            for (int i(0); i < nAlleles; i++) {
                sScore = score_perAllele[i];

                rPepScores_pt = &rPepScores[pepL][i];

                if (sScore > rPepScores_pt->back()) {
                    prctileScore_perL = 100.0 * double(rPepScores_pt->size()) / (rPepScores_pt->size()+1.0);
                    /* Didn't have any random peptide with a better score
                        than this peptide. Instead of giving a score of
                        exactly 100, I multiply by nRand / (nRand+1) to
                        account for the "precision" of our %Ranks (i.e. this
                        peptide is likely not the best of best binders but
                        we cannot distinguish between nearby values and if
                        adding one extra random peptide, maybe this peptide
                        would have been less good than best random one).
                        With this a very long peptide that has better score
                        than all similar length peptide will still show a
                        less good score than if it had been slightly shorter.
                        */
                } else {
                    prctileInd = 0;
                    while (sScore > rPepScores_pt->at(prctileInd)) {
                        prctileInd++;
                    }
                    if (prctileInd == 0) {
                        ca = 0;
                    } else {
                        ca = rPepScores_pt->at(prctileInd - 1);
                    }
                    cb = rPepScores_pt->at(prctileInd);
                    prctileScore_perL =
                        (prctileInd + (sScore - ca) / (cb - ca)) /
                        (rPepScores_pt->size()+1.0) * 100.0;
                    /* We use the equation with ca and cb here to linearly
                            interpolate the percentiles between the values
                        available from the random peptides, instead of just
                        using direct threshold values. We divide by nRand+1
                        for the same reason as above. */
                }

                prctileScore_perL = 100 - prctileScore_perL;
                /* This is to make the best scores have a value of 0 and
                    worst binders have a score of 100 - but in my
                    computations for transforming between perctile_perL to
                    perctile I was considering best cases had a score of 100.
                    It was easier to transform like that instead of deriving
                    again the equations of transformations (which would need
                    recreating new PWM def files with other alpha values as
                    well), but in the future it might still be good to also
                    modify these equations to make 1-2 less computations in
                    the code.
                */
                prctileScore_perL_perA[i] = prctileScore_perL;
            }
            peptides[p].set_PWM_res(prctileScore_perL_perA, score_perAllele,
                bestP1_perAllele, best_subspec_perAllele);
        } else { // The peptide was a random peptide
            for (int i(0); i < nAlleles; i++) {
                rPepScores[pepL][i][pRand] = score_perAllele[i];
            }
            pRand++;
            if (pRand == rPepScores[pepL][0].size()) {
                /* We have computed the scores from all random peptides of
                    current pepL; now we need to order these scores and will
                    start next iteration for the next peptide length. */
                for (int i(0); i < nAlleles; i++) {
                    sort(rPepScores[pepL][i].begin(),
                         rPepScores[pepL][i].end());
                }
                pRand = 0;
            }
        }
        undef_pep_score:;
        /* This is used by the goto from above, to skip all the computations
            for the peptides that have undefined scores (either too short or
            with an unknown aa). */
    }
}

void printPepResults(ofstream &outStream, const vector<Peptide> &Peptides,
    bool noContext, bool print_extra_results) {
    int bestInd, n_NA;
    const vector<double> *score_perAllele_pt;
    // A pointer to the scores, for easier reuse below.
    outStream << "Peptide\tContext\tBestAllele\t%Rank_best\tCore_best\tCoreP1_best\tSubSpec_best";

    for (auto cAllele : pars::alleles) {
        outStream << "\t%Rank_" << cAllele << "\tCoreP1_" << cAllele
                 << "\tSubSpec_" << cAllele;
    }
    n_NA = 5 + 3 * pars::alleles.size();
    // Tells number of NA values to use for peptides that couldn't be scored.

    if (print_extra_results) {
        for (auto cAllele : pars::alleles) {
            outStream << "\tScore_" << cAllele << "\t%RankPWM_" << cAllele
                << "\tScorePWM_" << cAllele;
        }
        /* When this info is asked we add it after the usual columns so that
            a script that would be based on column indices would still work. */
        n_NA += 3 * pars::alleles.size();
    }
    outStream << endl;

    for (auto cPep : Peptides) {
        outStream << cPep.getSequence() << "\t";
        if (!noContext){
            outStream << cPep.getContext();
            /* Only output the context string when we were including the
                context (otherwise it's always XXXXXXXXXXXX). Note that we'll
                still keep this column as empty when no context so that the
                output file always has the same columns. */
        }
        if (!cPep.na_score()) {
            // First determine which was the best allele.
            score_perAllele_pt = &cPep.get_nnets_rank();
            bestInd = min_element(score_perAllele_pt->begin(),
                                  score_perAllele_pt->end()) -
                      score_perAllele_pt->begin();
            /* min_element is a pointer to the element with min value from
                the given vector, and when removing the start pointer position
                this returns the index of this mîn element. */
            
            outStream << "\t" << pars::alleles[bestInd] << "\t" << setprecision(3)
                      << score_perAllele_pt->at(bestInd) << "\t"
                      << cPep.getSequence(cPep.get_P1()[bestInd], pars::L_core)
                      << "\t" << cPep.get_P1()[bestInd] + 1
                      << "\t" << cPep.get_subSpec_ID()[bestInd];
            /* We add 1 to bestP1 here and below to make the first AA has a
             * value of 1 and not 0 contrary to C++. */
            for (size_t i(0); i < pars::alleles.size(); i++) {
                outStream << "\t" << setprecision(3) << score_perAllele_pt->at(i)
                          << "\t" << cPep.get_P1()[i] + 1
                          << "\t" << cPep.get_subSpec_ID()[i];
            }
            if (print_extra_results){
                for (size_t i(0); i < pars::alleles.size(); i++) {
                    outStream
                        << "\t" << setprecision(3) << cPep.get_nnetsScore(i)
                        << "\t" << setprecision(3) << cPep.get_PWMrank_L(i)
                        << "\t" << setprecision(3) << cPep.get_PWMscore(i);
                }
            }
        } else {
            for (size_t i(0); i < n_NA; i++) {
                outStream << "\tNA";
                /* Repeat NA values the number of times needed to have same
                 * size than for the peptides that get scored. */
            }
        }
        outStream << endl;
    }
}


