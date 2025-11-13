//
// Created by raphael on 11/13/25.
//

#include "compare.h"

#include <stdexcept>

CompareKMers::CompareKMers(std::vector<char> nuclist, std::vector<size_t> reads, size_t kmersize)
    : nuclist(std::move(nuclist)), reads(std::move(reads)), kmersize(kmersize) {}

CompareKMers::CompareKMers(std::vector<char> nuclist, std::vector<size_t> reads)
    : nuclist(std::move(nuclist)), reads(std::move(reads)), kmersize(31) {}

// Setters
void CompareKMers::set_kmersize(size_t k) {
    kmersize = k;
}

// Getters
size_t CompareKMers::get_kmersize() const {
    return kmersize;
}

std::vector<char>& CompareKMers::get_nuclist() {
    return nuclist;
}

std::vector<size_t>& CompareKMers::get_reads() {
    return reads;
}

size_t CompareKMers::get_read_end_pos(const size_t read_idx) const
{
    if (read_idx < 0 || read_idx >= static_cast<size_t>(reads.size())) {
        throw std::out_of_range("read_idx is out of range");
    }
    
    if (read_idx == 0) {
        return 0;
    }
    
    return reads[read_idx];
}

size_t CompareKMers::get_n_kmers(const size_t ref) const {
    const size_t start = get_read_end_pos(ref-1);
    const size_t end = get_read_end_pos(ref);
    return (end - start) - (kmersize - 1);
}

size_t CompareKMers::get_all_n_kmers() const {
    return nuclist.size() - (kmersize - 1);
}

size_t CompareKMers::get_n_reads() const {
    return reads.size();
}

// Methods
size_t CompareKMers::compare_line(const size_t line1, const size_t line2) const
{
    size_t matches = 0;
    for (size_t i = 0; i < kmersize - 1; ++i)
    {
        if (nuclist[line1 + i + 1] == nuclist[line2 + i])
        {
            matches++;
        }
    }

    return matches == kmersize - 1 ? 1 : 0;
}

std::vector<size_t> CompareKMers::compare_lines(const size_t ref) const
{
    std::vector<size_t> results;
    results.reserve(kmersize);
    for (size_t i = 0; i < reads.size(); ++i)
    {
        if (i == ref) {
            results.push_back(1);  // Identical lines
            continue;
        }

        size_t line_start = reads[i];
        size_t cmp_result = compare_line(ref, line_start);
        results.push_back(cmp_result);
    }
    return results;
}

std::vector<std::vector<size_t>> CompareKMers::compare_all() const
{
    std::vector<std::vector<size_t>> all_results;
    all_results.reserve(reads.size());
    all_results.reserve(reads.size());
    for (size_t i = 0; i < reads.size(); ++i)
    {
        std::vector<size_t> line_results = compare_lines(i);
        all_results.push_back(std::move(line_results));
    }
    return all_results;
}