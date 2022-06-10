//
//  helper_general.hpp
//  MixMHC2pred
//
// File defining some project un-related helper functions/tools.
//
//  Created by Julien Racle on 01.03.19.
//  Copyright © 2019 CCB. All rights reserved.
//

#include "helper_general.hpp"

std::istream& safeGetline(std::istream& is, std::string& t)
{
    /* Function similar to getline function but that takes care of the various
        line endings appearing in different OS (i.e. \n or \r\n or ...).
        This was obtained from https://stackoverflow.com/a/6089413 */
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case std::streambuf::traits_type::eof():
            // Also handle the case when the last line has no line ending
            is.setstate(std::ios::eofbit);
            if(t.empty())
                is.setstate(std::ios::badbit);
            return is;
        default:
            t += (char)c;
        }
    }
}
