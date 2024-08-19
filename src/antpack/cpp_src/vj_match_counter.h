#ifndef VJ_MATCH_COUNTER_HEADER_H
#define VJ_MATCH_COUNTER_HEADER_H

#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <iostream>


namespace py = pybind11;

// The input sequence must have this length in order to be
// considered. Caller should already have extracted the IMGT
// standard positions which are all that is used for VJ
// assignment.
#define REQUIRED_SEQUENCE_LENGTH 128


class VJMatchCounter {
    public:
        VJMatchCounter(std::vector<std::string> geneSeqs,
                std::vector<std::string> geneNames);

        std::tuple<std::string, double> vjMatch(std::string query_sequence);
        std::string findVJSequenceByName(std::string query_name);
        std::tuple<std::vector<std::string>, std::vector<std::string>> getSeqLists();

    protected:
        std::vector<std::string> geneSeqs;
        std::vector<std::string> geneNames;
        std::map<std::string, int> namesToPositions;

};

#endif
