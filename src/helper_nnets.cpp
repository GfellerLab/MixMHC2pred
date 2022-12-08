//
//  helper_nnets.hpp
//  Definition file of the main functions / classes used for the nnets part of
//  MixMHC2pred.
//
//  Created by Julien Racle on 09.07.21.
//  Copyright © 2021 CCB. All rights reserved.
//

#include "helper_nnets.hpp"
#include "helper_MixMHC2pred.hpp"
#include "alleles_def.hpp"

// #########################
// ## Definition of the class methods.
// #########################

// ################################
// ## Helper functions definitions
// ################################
vector< vector< double > > encode_nn_in(size_t nData, const vector<string> &nn_in_var,
    const vector<string> &aaSeq, const valWithRange<double > &scores,
    const vector< valWithRange<int> > &one_hot){
        vector< vector< double > > nn_in(nData,
            vector< double >(nn_in_var.size(), 0.0));
        int cInd, iVal;
        
        if ((!scores.empty()) &&
            (find_index<string>(nn_in_var, scores.get_name(), cInd))){
                /* We use find_index so that the elements in the nn_in matrix
                    will be put at the correct column (given by the found cInd)
                    depending on the order of the variables that is given by
                    nn_in_var. Note that if the given score name isn't found
                    in the nn_in_var, then this scores won't be added to the
                    nn_in.
                */
                if (scores.size() != nData){
                    throw(string("When encoding nn_in, the scores have ")
                        + to_string(scores.size()) + " elements but nData = "
                        + to_string(nData) + "!");
                }
                for (size_t i(0); i < nData; i++){
                    nn_in[i][cInd] = (scores[i] - scores.min()) /
                        (scores.max() - scores.min());
                    if (nn_in[i][cInd] < 0){
                        nn_in[i][cInd] = 0;
                    } else if (nn_in[i][cInd] > 1){
                        nn_in[i][cInd] = 1;
                    }
                    // We do a min-max scaling of the score and make sure that
                    // they remain between 0 and 1 in case the scores in the
                    // training data range would have been different than in
                    // the current dataset.
                }
        }
        for (size_t j(0); j < one_hot.size(); j++){
            if (one_hot[j].size() != nData){
                throw(string("When encoding nn_in, the one_hot ")
                    + one_hot[j].get_name() + " has "                    
                    + to_string(one_hot[j].size()) + " elements but nData = "
                    + to_string(nData) + "!");
            }
            vector< string > col_names;
            for (int x(one_hot[j].min()); x <= one_hot[j].max(); x++){
                col_names.push_back(one_hot[j].get_name() + "_" + to_string(x));
            }
            vector< vector< double > > one_hot_vals(nData,
                vector< double >(col_names.size(), 0.0));
            for (size_t i(0); i < nData; i++){
                iVal = one_hot[j][i];
                // First put the values to the allowed range in case some aren't
                if (iVal < one_hot[j].min()){
                    iVal = one_hot[j].min();
                } else if (iVal > one_hot[j].max()){
                    iVal = one_hot[j].max();
                }
                cInd = iVal - one_hot[j].min();
                // Tells to which column the current value corresponds
                if (one_hot[j].use_multi1()){
                    for (int cCol(0); cCol <= cInd; cCol++){
                        one_hot_vals[i][cCol] = 1;
                    }
                } else {
                    one_hot_vals[i][cInd] = 1;
                }
                // We set values of 1 either to all columns before cInd or only
                // to cInd depending on the multi1 value.
            }
            // And we now put these one_hot_vals into the nn_in, depending on
            // which columns are present in nn_in_var (and in which order).
            for (size_t k(0); k < col_names.size(); k++){
                if (find_index<string>(nn_in_var, col_names[k], cInd)){
                    for (size_t i(0); i < nData; i++){
                        nn_in[i][cInd] = one_hot_vals[i][k];
                    }
                }
            }
        }

        if (!aaSeq.empty()) {
            if (aaSeq.size() != nData) {
                throw(string("When encoding nn_in, the aaSeq has ") +
                    to_string(aaSeq.size()) + " elements but nData = "
                    + to_string(nData) + "!");
            }
            size_t seqL(aaSeq[0].size());
            vector<string> col_names;
            for (size_t j(1); j <= seqL; j++) {
                // We start at 1 to seqL to use same indexing than in our R data.
                for (auto cAA : pars::blosumColumns){
                    col_names.push_back(string("AA_") + cAA + "_" + to_string(j));
                }
            }
            vector<vector<double>> aa_encoded(nData, vector<double>());
            for (size_t i(0); i < nData; i++) {
                if (aaSeq[i].size() != seqL){
                    throw(string("The sequence " + to_string(i) + " (" + aaSeq[i] +
                        ") doesn't have same size than other sequences of L: "
                        + to_string(seqL)));
                }
                aa_encoded[i].reserve(col_names.size());
                for (size_t j(0); j < seqL; j++) {
                    aa_encoded[i].insert(aa_encoded[i].end(), 
                        pars::blo_plus_1[aaSeq[i][j]].begin(),
                        pars::blo_plus_1[aaSeq[i][j]].end());
                    // Add the various blo_plus_1 elements of the current AA of
                    // the sequence.
                }
            }
            // And we now put this aa_encoded into the nn_in, depending on
            // which columns are present in nn_in_var (and in which order).
            for (size_t k(0); k < col_names.size(); k++) {
                if (find_index<string>(nn_in_var, col_names[k], cInd)) {
                    for (size_t i(0); i < nData; i++) {
                        nn_in[i][cInd] = aa_encoded[i][k];
                    }
                }
            }
        }

        return nn_in;
}


vector< vector< double > > nn_relu(const vector< vector<double> > &x,
    const vector< vector<double> > &w, const vector<double> &b){
  if ((x[0].size() != w.size()) || (w[0].size() != b.size()))
    throw(string("The dimensions of x, w and b aren't coherent."));
  
  vector< vector<double> > y(x.size(), vector<double>(b.size()));
  double val;
  for (size_t k(0); k < x.size(); k++) {
      for (size_t i(0); i < b.size(); i++) {
          val = b[i];
          for (size_t j(0); j < w.size(); j++) {
              val += w[j][i] * x[k][j];
          }
          if (val < 0) {
              val = 0;
          }
          y[k][i] = val;
      }
  }
  return y; 
}

vector<vector<double>> nn_logistic(const vector<vector<double>> &x,
    const vector<vector<double>> &w, const vector<double> &b) {
    if ((x[0].size() != w.size()) || (w[0].size() != b.size()))
        throw(string("The dimensions of x, w and b aren't coherent."));

    vector<vector<double>> y(x.size(), vector<double>(b.size()));
    double val;
    for (size_t k(0); k < x.size(); k++) {
        for (size_t i(0); i < b.size(); i++) {
            val = b[i];
            for (size_t j(0); j < w.size(); j++) {
                val += w[j][i] * x[k][j];
            }
            y[k][i] = 1.0 / (1. + exp(-val));
        }
    }
    return y;
}

vector<vector<double>> nn_softmax(const vector<vector<double>> &x,
    const vector<vector<double>> &w, const vector<double> &b) {
    if ((x[0].size() != w.size()) || (w[0].size() != b.size()))
        throw(string("The dimensions of x, w and b aren't coherent."));

    vector<vector<double>> y(x.size(), vector<double>(b.size()));
    double val, sum_y_kj;
    for (size_t k(0); k < x.size(); k++) {
        sum_y_kj = 0;
        for (size_t i(0); i < b.size(); i++) {
            val = b[i];
            for (size_t j(0); j < w.size(); j++) {
                val += w[j][i] * x[k][j];
            }
            sum_y_kj += (y[k][i] = exp(val));
        }
        for (size_t i(0); i < b.size(); i++) {
            y[k][i] /= sum_y_kj;
        }
    }
    return y;
}

void pred_neural_networks(string nnetsDefFile, vector<string> alleles,
    vector<Peptide> &Peptides, bool noContext) {
    size_t nPep(Peptides.size());
    int cP1;
    double log_min, cRank, pRank_T;
    vector<double> PWM_range, rank_thr;
    vector<int> pepL_range, P1_range;
    vector<string> nn_in_var;
    vector<vector<vector<double>>> nnetsWeights;
    vector<vector<double>> nnetsBiases;
    map<int, vector<double>> sc_thr;
    map<int, double> f_L;

    // Getting the weights implemented in the current nnets.
    get_nnets_def(nnetsDefFile, PWM_range, pepL_range, P1_range, log_min,
        nn_in_var, nnetsWeights, nnetsBiases, rank_thr, sc_thr,
        pRank_T, f_L, noContext);
    
    // Defining the variables that'll be used for the nn-encoding.
    vector<string> pepBords;
    pepBords.resize(nPep);
    valWithRange<double> scores("PWM_log", log(PWM_range[0] + log_min),
                                log(PWM_range[1] + log_min), nPep);
    vector<valWithRange<int>> one_hot;
    one_hot.push_back(
        valWithRange<int>("pepL", pepL_range[0], pepL_range[1], nPep));
    one_hot.push_back(
        valWithRange<int>("P1", P1_range[0], P1_range[1], nPep));
    for (size_t Aind(0); Aind < alleles.size(); Aind++){
        // Repeat the computation for each allele separately.
        for (size_t i(0); i < nPep; i++) {
            scores[i] = log(Peptides[i].get_PWMrank_L(Aind) + log_min);
            pepBords[i] = Peptides[i].getContext();
            one_hot[0][i] = Peptides[i].get_pepL();

            cP1 = Peptides[i].get_P1_centered(Aind);
            if (Peptides[i].get_subSpec_ID()[Aind] < 0){
                /* peptide is a reverse binder for given allele, need to
                    revert its P1 to count from C-terminal of the peptide. */
                cP1 = -cP1;
            }
            one_hot[1][i] = cP1;
        }

        vector<vector<double>> nn_in =
            encode_nn_in(nPep, nn_in_var, pepBords, scores, one_hot);
        /* Note that we could win a little bit of time by keeping the pepBords
            one-hot encoding between the different nnets-repetitions, but this
            shouldn't be a big time-saver. */

        for (size_t i(0); i < nnetsWeights.size() - 1; i++) {
            nn_in = nn_relu(nn_in, nnetsWeights[i], nnetsBiases[i]);
            /* Repeat this possibly multiple times depending on the number of
              hidden layers. And for the last layer, we'll use another
              activation function.
            */
        }
        nn_in = nn_logistic(nn_in, nnetsWeights.back(), nnetsBiases.back());
        if (nn_in[0].size() > 1) {
            throw(string(
                      "There should be a single column for the nnet output. ") +
                  "The nnet architecture is probably different than expected.");
        }

        // And now we'll transform these scores to %Rank scores.
        for (size_t i(0); i < nPep; i++) {
            if (!isnan(nn_in[i][0]) && !Peptides[i].na_score()){
                // The score was nan if peptide is of bad size or contains bad AA.
                cRank = transf_toPrank(1.0 - nn_in[i][0],
                    Peptides[i].get_pepL(), rank_thr,
                    sc_thr, pRank_T, f_L);
                Peptides[i].set_nnets_rank(cRank, nn_in[i][0], Aind);
            } else {
                Peptides[i].set_na_score();
            }
        }
    }
}

double transf_toPrank(double cScore, int pepL,
    const vector<double> &rank_thr,
    const map<int, vector<double>> &sc_thr,
    double pRank_T, const map<int, double> &f_L) {
    
    // We first transform the scores to %Rank_perL
    int cInd;
    double pRank, a_L, d_L, e_L, p_L, ref_L(15);
    // ref_L is the reference length towards which the %ranks are transformed.
    const vector<double> &c_sc_thr(sc_thr.at(pepL));

    if (cScore < c_sc_thr[0]) {
        pRank = rank_thr[0];
    } else if (cScore > c_sc_thr.back()) {
        pRank = rank_thr.back();
        // Put to the %Ranks limits the scores that are better or worst than
        // all the tested values.
    } else {
        cInd = 1;
        while (cScore > c_sc_thr[cInd]) {
            cInd++;
        }

        pRank = (rank_thr[cInd - 1] +
                (rank_thr[cInd] - rank_thr[cInd - 1]) *
                    (cScore - c_sc_thr[cInd - 1]) /
                    (c_sc_thr[cInd] - c_sc_thr[cInd - 1]));
        // For the intermediate scores, we interpolate between the given
        // thresholds
    }

    /* And now we transform the %Ranks_perL to %Ranks that follow
        the expected pepL distribution. */
    if (pepL != ref_L){
        // No need to transform when pepL = ref_L.
        a_L = f_L.at(ref_L) / f_L.at(pepL);
        p_L = pRank_T / a_L;
        if (pRank <= p_L){
            pRank *= a_L;
        } else {
            d_L = (100. - pRank_T) / (100. - p_L);
            e_L = 100. - 100. * d_L;
            pRank = d_L * pRank + e_L;
        }
    }
  return pRank;
}
