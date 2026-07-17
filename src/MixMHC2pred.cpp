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
#include <chrono>
#include <string>
#include "helper_MixMHC2pred.hpp"
#include "helper_nnets.hpp"
#include "alleles_def.hpp"
#include "../external/argparse/argparse.hpp"
#include "../external/whereami/src/whereami.h"
using namespace std;

int main(int argc, const char *argv[]) {
    // Define the options that can be passed to the program
	ArgumentParser parser;
	parser.addArgument("-i", "--input", 1, false);
	parser.addArgument("-o", "--output", 1, false);
    parser.addArgument("-a", "--alleles", '+', false);
    // parser.addArgument("-s", "--ws_scoring", 1); // Not making this an option anymore
    parser.addArgument("--no_context", 0);
    parser.addArgument("-f", "--allelesFolder", '+');
    parser.addArgument("--nnetsFolder", 1);
    parser.addArgument("-e", "--extra_out", 0);
    parser.addArgument("--batch_size", 1);
    parser.addHiddenArgument("--ignore_version", '*');
    /* ignore_version is used to avoid checking that the PWMdef file
        corresponds to pars::version. This is an advanced argument used for
        testing different PWM def files and it should not be used in the
        standard case, this is why it is hidden. It can be given with a value
        to replace the version name by this value. */
    // parser.addArgument("--allowPredWithMissingAlleles", 0);

    vector<string> alleles, allelesDefFolders = {"default"}, nnetsDefFiles;
    string pepFile, oFile, rPepFolder, execFolder,
        nnetsDefRoot("nnetsWeights"), nnetsFolder("default"), contextType("full");
    bool needAllAlleles, noContext, print_extra_results;
    int ws_scoring(1), batch_size(1e6), nBatches;
    time_t cur_time;

	
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
    if (parser.exists("nnetsFolder")){
        nnetsFolder = parser.retrieve<string>("nnetsFolder");
    }
    print_extra_results = parser.exists("extra_out");
    if (parser.exists("batch_size")){
        batch_size = stoi(parser.retrieve<string>("batch_size"));
    }

    try {
        pars::ignore_version = parser.exists("ignore_version");
        if (pars::ignore_version && (parser.count("ignore_version") >= 1)){
            pars::version = parser.retrieve<vector<string>>("ignore_version")[0];
            if (parser.count("ignore_version") > 1){
                throw(string("'ignore_version' was given with more than 1 ") +
                    "argument.");
            }
        }

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

        if (noContext){
            contextType = "none";
        }

        if (ws_scoring != 0 && ws_scoring != 1 && ws_scoring != 2 &&
        ws_scoring != 3 && ws_scoring != 4){
            throw(string("ws_scoring can only take values 0,1,2,3 and 4"));
        }
        
        /* Checking the nnets def files exist and getting version name if new 
            file format. */
        nnetsDefFiles = nnets_def_files(nnetsFolder, nnetsDefRoot, execFolder);
        
        cur_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
        std::cout << "[" << put_time(localtime(&cur_time), "%Y-%m-%d %H:%M:%S")
			<< "] Runing MixMHC2pred (v" << pars::version
             << ") for peptide file: " << pepFile << endl;

        vector<Peptide> peptides_full, peptides;
        set<int> pepL_used; // Tells the pepL present in input file.

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
        
        cur_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
        std::cout << "[" << put_time(localtime(&cur_time), "%Y-%m-%d %H:%M:%S")
			<< "] Imported " << peptides.size() << " peptides." << endl;
        
        pars::alleles = alleles;
        ofstream outStream(oFile);

        /* To avoid memory issues, we'll split the computations in multiple
            batches when there are a lot of peptides to compute. */
        nBatches = ceil(double(peptides.size()) / batch_size);
        if (nBatches > 1) {
            peptides_full = peptides;
            std::cout << "Computing the scores in " << nBatches
                << " batches of " << batch_size << " peptides." << endl;
        }

        for (int cBatch(0); cBatch < nBatches; cBatch++) {
            if (nBatches > 1) {
                peptides.assign(peptides_full.begin() + cBatch*batch_size,
                    peptides_full.begin() + 
                        min((cBatch+1)*batch_size, (int)peptides_full.size()));
                std::cout << "--- Starting batch " << cBatch+1 << endl;
            }

            cur_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
            std::cout << "[" << put_time(localtime(&cur_time), "%Y-%m-%d %H:%M:%S")
                << "] Computing the PWM-based scores from each peptide." << endl;

            // First get the PWMs-scores. We redefine some variables based on what
            // is now always used for the nnets (!!! will soon delete these, but
            // still kept for the moment).
            ws_scoring = 1;
            needAllAlleles = true;
            /* !!! Will need to update the way to treat needAllAlleles - currently
                with the nnets, I always force this to true.
                */
            comp_pepScores(peptides, ws_scoring, needAllAlleles, rPepFolder,
                        pepL_used, allelesDefFolders);
            cur_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
            std::cout << "[" << put_time(localtime(&cur_time), "%Y-%m-%d %H:%M:%S")
                << "] Computing now the full scores from each peptide." << endl;
            for (size_t c_nnet(0); c_nnet < pars::n_nnets; c_nnet++){
                pred_neural_networks(nnetsDefFiles[c_nnet], alleles, peptides, noContext);
            }
            cur_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
            std::cout << "[" << put_time(localtime(&cur_time), "%Y-%m-%d %H:%M:%S")
                << "] Saving the results in the output file." << endl;

            if (cBatch == 0){
                // Making first the header from the results file
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
                outStream << "\n# nnets definition folder: ";
                string cFolder(nnetsFolder);
                if (cFolder.substr(0, execFolder.size()) == execFolder){
                    cFolder = "exec:" + cFolder.substr(execFolder.size());
                    /* When folder based on path to executable, we replace this path
                        to make it shorter. */
                }
                outStream << "'" << cFolder << "'"
                    << "\n# Context type: " << contextType
                    << "\n#\n# MixMHC2pred is freely available for academic users." << endl
                    << "# Private companies should contact Nadette Bulgin (nbulgin@lcr.org) "
                    << "at the Ludwig Institute\n#  for Cancer Research Ltd "
                    << "for commercial licenses.\n#" << endl
                    << "# To cite MixMHC2pred, please refer to:" << endl
                    << "# Racle, J., et al., Machine learning predictions of MHC-II "
                    << "specificities\n#  reveal alternative binding mode of class II "
                    << "epitopes. Immunity (2023).\n# and\n"
                    << "# Racle, J., et al. Robust prediction of HLA class II epitopes "
                    << "by deep motif\n#  deconvolution of immunopeptidomes. Nat. "
                    << "Biotechnol. 37, 1283-1286 (2019).\n#\n" << string(20, '#')
                    << endl;
            }

            printPepResults(outStream, peptides, noContext, print_extra_results,
                (cBatch == 0));
        }

        cur_time = chrono::system_clock::to_time_t(chrono::system_clock::now());
        std::cout << "[" << put_time(localtime(&cur_time), "%Y-%m-%d %H:%M:%S")
			<< "] Finished the computations." << endl;

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
