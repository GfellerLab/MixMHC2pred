//
//  MixMHC2pred.cpp
//  Predicts the binding of peptides to HLA class II.
//
//  Created by Julien Racle on 08.03.18.
//  Copyright © 2019 CCB. All rights reserved.
//

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include "helper_MixMHC2pred.hpp"
#include "helper_nnets.hpp"
#include "../external/argparse/argparse.hpp"
#include "../external/whereami/src/whereami.h"
using namespace std;

int main(int argc, const char *argv[]) {
    pars::version = "2.0.2";

    // Define the options that can be passed to the program
	ArgumentParser parser;
	parser.addArgument("-i", "--input", 1, false);
	parser.addArgument("-o", "--output", 1, false);
    parser.addArgument("-a", "--alleles", '+', false);
    // parser.addArgument("-s", "--ws_scoring", 1); // Not making this an option anymore
    parser.addArgument("--no_context", 0);
    parser.addArgument("-f", "--allelesFolder", '+');
    parser.addArgument("-e", "--extra_out", 0);
    parser.addHiddenArgument("--ignore_version", '*');
    /* ignore_version is used to avoid checking that the PWMdef file
        corresponds to pars::version. This is an advanced argument used for
        testing different PWM def files and it should not be used in the
        standard case, this is why it is hidden. It can be given with a value
        to replace the version name by this value. */
    // parser.addArgument("--allowPredWithMissingAlleles", 0);

    vector<string> alleles, allelesDefFolders = {"default"};
    string pepFile, oFile, rPepFolder, execFolder,
        nnetsDefRoot("nnetsWeights"), nnetsDefFile, contextType("full");
    bool needAllAlleles, noContext, print_extra_results;
    int ws_scoring(1);

	
    // Parsing the arguments.
    parser.parse(argc, argv);
    pepFile = parser.retrieve<string>("input");
	oFile = parser.retrieve<string>("output");
    alleles = parser.retrieve<vector<string>>("alleles");
    // needAllAlleles = !parser.exists("allowPredWithMissingAlleles");
    noContext = parser.exists("no_context");
    // if (parser.exists("ws_scoring"))
    //     ws_scoring = stoi(parser.retrieve<string>("ws_scoring"));
    if (parser.exists("allelesFolder")){
        allelesDefFolders = parser.retrieve<vector<string>>("allelesFolder");
    }
    print_extra_results = parser.exists("extra_out");

    try {
        pars::ignore_version = parser.exists("ignore_version");
        if (pars::ignore_version && (parser.count("ignore_version") >= 1)){
            pars::version = parser.retrieve<vector<string>>("ignore_version")[0];
            if (parser.count("ignore_version") > 1){
                throw(string("'ignore_version' was given with more than 1 ") +
                    "argument.");
            }
        }

        if (noContext){
            contextType = "none";
        }

        if (ws_scoring != 0 && ws_scoring != 1 && ws_scoring != 2 &&
        ws_scoring != 3 && ws_scoring != 4)
            throw(string("ws_scoring can only take values 0,1,2,3 and 4"));

        std::cout << "Runing MixMHC2pred (v" << pars::version
             << ") for peptide file: " << pepFile << endl;

        vector<Peptide> peptides;
        set<int> pepL_used; // Tells the pepL present in input file.

        char *path = NULL;
        int path_length, dirname_length;
        path_length = wai_getExecutablePath(NULL, 0, &dirname_length);
        try {
            if (path_length > 0) {
                path = (char *)malloc(path_length + 1);
                if (!path)
                    throw(1);
                wai_getExecutablePath(path, path_length, &dirname_length);
                execFolder = string(path, dirname_length);
                // Only keep the first chars corresponding to the directory
                // name.
                free(path);
            } else
                throw(1);
        } catch (...) {
            throw(string("There is an error locating the executable of") +
                  " MixMHC2pred. It is maybe not in a standard location.");
        }

        // Possibly redefine the path to the alleles definition folders:
        for (string &cFolder : allelesDefFolders) {
            if (cFolder == "default") {
                cFolder = execFolder + "/PWMdef/";
            } else if (cFolder.substr(0, 5) == "exec:") {
                cFolder =
                    execFolder + "/" + cFolder.substr(5) + "/";
            } else {
                cFolder += "/";
                // Make sure path ends with a '/'.
            }
        }

        rPepFolder = execFolder + "/rpep/";
        // Loading the peptides and computing their scores.
        pepData_import(pepFile, peptides, pepL_used, contextType, alleles.size());
        std::cout << "Imported " << peptides.size() << " peptides." << flush;

        std::cout << " Computing now the PWM-based scores from each peptide." << endl;

        // First get the PWMs-scores. We redefine some variables based on what
        // is now always used for the nnets (!!! will soon delete these, but
        // still kept for the moment).
        ws_scoring = 1;
        needAllAlleles = true;
        /* !!! Will need to update the way to treat needAllAlleles - currently
            with the nnets, I always force this to true.
            */
        pars::alleles = alleles;
        comp_pepScores(peptides, ws_scoring, needAllAlleles, rPepFolder,
                       pepL_used, allelesDefFolders);
        std::cout << "Computing now the full scores from each peptide." << endl;
        for (size_t c_nnet(1); c_nnet <= pars::n_nnets; c_nnet++){
            nnetsDefFile = execFolder + "/nnetsDef/" + nnetsDefRoot + "_"
                + to_string(c_nnet) + ".txt";
            pred_neural_networks(nnetsDefFile, alleles, peptides, noContext);
        }
        std::cout << "Saving the results in the output file." << endl;

        // Making first the header from the results file
        ofstream outStream(oFile);
        outStream << string(20, '#') << "\n# Output from MixMHC2pred (v"
            << pars::version << ")" << endl
                  << "# Input file: " << pepFile << endl
                  << "# Alleles: ";
        for (string cAllele : alleles) {
            outStream << cAllele << ", ";
        }
        outStream.seekp(-2, ios_base::cur);  // Remove ', ' from the end of outStream.
        outStream << "\n# Alleles folder: ";
        for (string cFolder : allelesDefFolders){
            if (cFolder.substr(0, execFolder.size()) == execFolder){
                cFolder = "exec:" + cFolder.substr(execFolder.size());
                /* When folder based on path to executable, we replace this path
                    to make it shorter. */
            }
            outStream << "'" << cFolder << "'";
            if (cFolder != allelesDefFolders.back()){
                outStream << "; ";
            }
        }
        outStream << "\n# Context type: " << contextType
            << "\n#\n# MixMHC2pred is freely available for academic users." << endl
            << "# Private companies should contact nbulgin@lcr.org "
            << "at the Ludwig Institute\n#  for Cancer Research Ltd "
            << "for commercial licenses.\n#" << endl
            << "# To cite MixMHC2pred, please refer to:" << endl
            << "# Racle, J., et al., Machine learning predictions of MHC-II "
            << "specificities\n#  reveal alternative binding mode of class II "
            << "epitopes. bioRxiv (2022).\n# and\n"
            << "# Racle, J., et al. Robust prediction of HLA class II epitopes "
            << "by deep motif\n#  deconvolution of immunopeptidomes. Nat. "
            << "Biotechnol. 37, 1283-1286 (2019).\n#\n" << string(20, '#')
            << endl;

        printPepResults(outStream, peptides, noContext, print_extra_results);

    
        std::cout << "Finished the computations." << endl;

    } catch(string& errMsg){
        cerr << string(15, '#') << endl << string(15, '#')
            << "\n** Something went wrong!: " << errMsg << endl;
        remove(oFile.c_str());
        /* There was an error and the oFile has possibly already been started
            so we delete it. */
        return 1;
    }

    return 0;
}
