//
// Created by raphael on 11/13/25.
//

#ifndef ASSEMBLEUR_COMPARE_H
#define ASSEMBLEUR_COMPARE_H
#include <string>
#include <vector>



struct CompareKMers {
    std::vector<char> nuclist;
    std::vector<size_t> reads;
    size_t kmersize = 31;

    CompareKMers(std::vector<char> nuclist, std::vector<size_t> reads, size_t kmersize);
    CompareKMers(std::vector<char> nuclist, std::vector<size_t> reads);

    // Setters
    void set_kmersize(size_t k);

    // Getters
    size_t get_kmersize() const;
    std::vector<char>& get_nuclist();
    std::vector<size_t>& get_reads();
    size_t get_read_end_pos(size_t read_idx) const;
    size_t get_n_kmers(size_t ref) const;
    size_t get_all_n_kmers() const;
    size_t get_n_reads() const;

    // Methods
    size_t                           compare_line(size_t line1, size_t line2) const;
    std::vector<size_t>              compare_lines(size_t ref) const;
    std::vector<std::vector<size_t>> compare_all() const;
};



#endif //ASSEMBLEUR_COMPARE_H