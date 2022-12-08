//
//  helper_MixMHC2pred.hpp
//  Declaration file of the main functions / classes used for MixMHC2pred.
//
//  Created by Julien Racle on 30.01.19.
//  Copyright © 2019 CCB. All rights reserved.
//

#ifndef helper_mixMHC2pred_hpp
#define helper_mixMHC2pred_hpp

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

using namespace std;

/* First defining a "class" that will contain various parameter values. This
	class is used to define kind of global variables, through static class
	variables. */
class pars {
public:
	/* Class containing various parameters describing the system. */
	static size_t n_aa, context_L, L_core, n_nnets;
    static map<char, int> aa_to_ind;
    static map< char, vector<double> > blo_plus_1;
    static vector<char> blosumColumns;
    static vector<string> alleles;
    static string version;
    static bool ignore_version;
	/* 
         context_L is the number of AA per "part of the context" (i.e. # before
            N-term, # just after it, # just before C-term and # just after it),
            so that the total context length will be 4 times this number.
         n_aa is the number of different amino acids in this alphabet (including
            only the 'true' AA, i.e. without the allowed unknown AA (X) or gap
            symbol (-).
        n_nnets is the number of nnets repetitions that have been used.
        aa_to_ind is used in the PWM-based computations to know which row of
            the PWM to use in function of which is the current AA (for both
            X and '-' we'll add an extra row to the PWM corresponding to an
            average value).
        blo_plus_1 is the blosum62 + one-hot transformation matrix, given as a map, so that
            we can call blo_plus_1["F"] for example to have all elements to encode
            this F AA in blosum.
        blosumColumns tells which AAs are present in the 'columns' of this
            blo_plus_1 matrix (it includes the true AAs and the gap).
        alleles gives the list of alleles asked in the input.
        version stores the current version number of MixMHC2pred.
        ignore_version is used to avoid checking correct version number.
	 */
private:
	pars(){}
	/* Make the class constructor private to forbid the creation of objects
	from this class as it only contains static elements. */
};

class Peptide {
    /* Definition of a class that will store the results from each peptide
    and that can then possibly output these results. */
  public:
    Peptide(int nAlleles, string sequence, string context="", bool innerContext=false);
    string getSequence(int from = 0, int length = INT_MAX) const;
    /* To return a string with the aa sequence of the peptide, starting
    at position from and of the given length (or max length). */
    string getContext() const;
    int get_pepL() const;
    
    const vector< int >& get_P1() const;
    int get_P1_centered(int Aind) const;
    const vector< int >& get_subSpec_ID() const;
    const vector< double >& get_nnets_rank() const;
    double get_PWMrank_L(int Aind) const;
    double get_PWMscore(int Aind) const;
    double get_nnetsScore(int Aind) const;

    void set_PWM_res(const vector< double > &PWMrank_L,
        const vector< double > &PWMrawScore, const vector< int > &P1,
        const vector< int > &best_subSpec_ids);
    void set_nnets_rank(double c_nnets_rank, double c_nnets_raw, int allele_ind);
    void set_na_score();
    const bool na_score();


  protected:
    string sequence, pepCon;
    vector<int> P1, best_subSpec_ids;
    // best_subSpec_ids tells which is the best sub-specificity ID per allele.
    vector<double> nnets_rank, nnets_rawScore, PWMrank_L, PWMrawScore;
    /* nnets_rank is the %Rank after the nnets, PWMrank_L is the peptide
        %Rank_perL after doing the part of the PWMs (which is used as input in
        the 2nd neural network) and the ..._rawScores are before %Rank
        transformations). */
    bool score_isNA;
    // Tells that the peptide couldn't get its score (e.g. bad size or bad AA)
};

void pepData_import(string pepFile, vector<Peptide> &peptides,
    set<int> &pepL_used, string contextType, int nAlleles);
/* To import the amino acid sequences of each peptide found in the given input
 file and return a set containing the peptide lengths present in this file
 to only do random computations against these ones.*/

void comp_pepScores(vector<Peptide> &peptides,
                    int ws_scoring, bool needAllAlleles,
                    string rpepFolder, set<int> &pepL_used,
                    const vector<string> &allelesDefFolders);
/* This is the main function called to compute the scores of the peptides.
    This scores are directly updated in the 'peptides'. */

void printPepResults(ofstream &outStream, const vector<Peptide> &Peptides,
    bool noContext, bool print_extra_results);
/* Print into the outStream the score, core binding region and other
    information from the current peptide. */

// ################################
// ## Definition of the inline class methods
// ################################
inline string Peptide::getSequence(int from, int length) const{
    return sequence.substr(from, length);
}
inline string Peptide::getContext() const{
    return pepCon;
}
inline int Peptide::get_pepL() const {
    return sequence.size();
}
inline const vector< int >& Peptide::get_P1() const {
    return P1;
}
inline const vector< int >& Peptide::get_subSpec_ID() const {
    return best_subSpec_ids;
}
inline const vector< double >& Peptide::get_nnets_rank() const{
    return nnets_rank;
}
inline double Peptide::get_PWMrank_L(int Aind) const{
    return PWMrank_L[Aind];
}
inline double Peptide::get_PWMscore(int Aind) const{
    return PWMrawScore[Aind];
}
inline double Peptide::get_nnetsScore(int Aind) const{
    return nnets_rawScore[Aind];
}

inline void Peptide::set_PWM_res(const vector< double > &c_PWMrank_L,
    const vector< double > &c_PWMrawScore, const vector< int > &cP1,
    const vector< int > &c_best_subSpec_ids){
    PWMrank_L = c_PWMrank_L;
    PWMrawScore = c_PWMrawScore;
    P1 = cP1; // These indices start at 0 (aren't centered values).
    best_subSpec_ids = c_best_subSpec_ids;
    nnets_rank.resize(PWMrank_L.size());
}

inline void Peptide::set_nnets_rank(double c_nnets_rank, double c_nnets_raw,
    int allele_ind){
    nnets_rank[allele_ind] += c_nnets_rank / double(pars::n_nnets);
    nnets_rawScore[allele_ind] += c_nnets_raw / double(pars::n_nnets);
    /* Each nnets repetition will contribute to a fraction of the score (i.e.
        the final result corresponds to the average ranks of the different
        repetitions). */
}

inline void Peptide::set_na_score() {
    score_isNA = true;
}
inline const bool Peptide::na_score() {
    return score_isNA;
}


#endif /* helper_mixMHC2pred_hpp */