//
//  helper_general.hpp
//  MixMHC2pred
//
// File defining some project un-related helper functions/tools.
//
//  Created by Julien Racle on 01.03.19.
//  Copyright © 2019 CCB. All rights reserved.
//

#ifndef newhelper_general_h
#define newhelper_general_h

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

std::istream& safeGetline(std::istream& is, std::string& t);


template< class T> bool find_index(const vector<T> &v, const T &x, int &index)
{
    /* Function inspired by https://www.geeksforgeeks.org/how-to-find-index-of-a-given-element-in-a-vector-in-cpp/
        to get the index at which an element is found in a given vector (or
        returning false if it isn't found) */
	auto it = find(v.begin(), v.end(), x);

	// If element was found
	if (it != v.end()) {
		// calculating the index of x
		index = it - v.begin();
		return true;
	} else {
		return false;
	}
}

#endif /* helper_general_h */
