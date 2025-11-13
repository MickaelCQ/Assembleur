#include <gtest/gtest.h>
#include <vector>
#include <stdexcept>
#include "../src/compare.h"

// A test fixture to set up common data for all tests
class CompareKMersTest : public ::testing::Test {
protected:
    std::vector<char> test_nuclist;
    std::vector<size_t> test_reads;
    CompareKMers* kmer_comparer;

    // SetUp is called before each test
    void SetUp() override {
        // Mock data based on prompt.
        // Reads: "ACGT", "AACG", "GTAC"
        // We add extra "GGG" at the end to prevent compare_lines from
        // reading out of bounds with the provided (flawed) logic.
        test_nuclist = {
            'A', 'C', 'G', 'T', // Read 0 (k-mers: 0="ACGT")
            'A', 'A', 'C', 'G', // Read 1 (k-mers: 4="AACG")
            'G', 'T', 'A', 'C', // Read 2 (k-mers: 8="GTAC")
            'G', 'G', 'G'       // Buffer
        }; // Total size 15
        
        // Cumulative end positions
        test_reads = {4, 8, 12}; // 3 reads

        // Initialize with k=4 for predictable logic tests
        kmer_comparer = new CompareKMers(test_nuclist, test_reads, 4);
    }

    // TearDown is called after each test
    void TearDown() override {
        delete kmer_comparer;
    }
};

// --- Test Constructors ---

TEST_F(CompareKMersTest, ConstructorWithKmerSize) {
    EXPECT_EQ(kmer_comparer->get_kmersize(), 4);
    EXPECT_EQ(kmer_comparer->get_nuclist().size(), 15);
    EXPECT_EQ(kmer_comparer->get_reads().size(), 3);
}

TEST_F(CompareKMersTest, ConstructorDefaultKmerSize) {
    // Delete the one from SetUp and create a new one
    delete kmer_comparer;
    kmer_comparer = new CompareKMers(test_nuclist, test_reads);
    
    EXPECT_EQ(kmer_comparer->get_kmersize(), 31); // Test default
}

// --- Test Getters and Setters ---

TEST_F(CompareKMersTest, SetAndGetKmerSize) {
    kmer_comparer->set_kmersize(10);
    EXPECT_EQ(kmer_comparer->get_kmersize(), 10);
}

TEST_F(CompareKMersTest, GetNuclist) {
    EXPECT_EQ(kmer_comparer->get_nuclist().size(), 15);
    EXPECT_EQ(kmer_comparer->get_nuclist()[0], 'A');
}

TEST_F(CompareKMersTest, GetReads) {
    EXPECT_EQ(kmer_comparer->get_reads().size(), 3);
    EXPECT_EQ(kmer_comparer->get_reads()[0], 4);
}

TEST_F(CompareKMersTest, GetNReads) {
    EXPECT_EQ(kmer_comparer->get_n_reads(), 3);
}

// --- Test Position and Count Getters ---
// NOTE: These tests validate the *actual* implementation,
// even if its logic is confusing.

TEST_F(CompareKMersTest, GetReadEndPos) {
    // Per implementation: returns 0 for idx 0, and reads[idx] otherwise
    EXPECT_EQ(kmer_comparer->get_read_end_pos(0), 0);
    EXPECT_EQ(kmer_comparer->get_read_end_pos(1), 8);  // returns reads[1]
    EXPECT_EQ(kmer_comparer->get_read_end_pos(2), 12); // returns reads[2]
}

TEST_F(CompareKMersTest, GetReadEndPosOutOfRange) {
    // Test upper bound
    EXPECT_THROW(kmer_comparer->get_read_end_pos(3), std::out_of_range);
    // Test with a large number
    EXPECT_THROW(kmer_comparer->get_read_end_pos(99), std::out_of_range);
}

TEST_F(CompareKMersTest, GetNKmers) {
    kmer_comparer->set_kmersize(4);
    
    // get_n_kmers(0) will throw because get_read_end_pos(0-1) is called
    EXPECT_THROW(kmer_comparer->get_n_kmers(0), std::out_of_range);

    // get_n_kmers(1): start=pos(0)=0, end=pos(1)=8. len=8. k=4.
    // (8 - 0) - (4 - 1) = 8 - 3 = 5
    EXPECT_EQ(kmer_comparer->get_n_kmers(1), 5);

    // get_n_kmers(2): start=pos(1)=8, end=pos(2)=12. len=4. k=4.
    // (12 - 8) - (4 - 1) = 4 - 3 = 1
    EXPECT_EQ(kmer_comparer->get_n_kmers(2), 1);
}

TEST_F(CompareKMersTest, GetAllNKmers) {
    kmer_comparer->set_kmersize(4);
    // nuclist.size() = 15. k = 4.
    // 15 - (4 - 1) = 15 - 3 = 12
    EXPECT_EQ(kmer_comparer->get_all_n_kmers(), 12);
    
    kmer_comparer->set_kmersize(10);
    // 15 - (10 - 1) = 15 - 9 = 6
    EXPECT_EQ(kmer_comparer->get_all_n_kmers(), 6);
}

// --- Test Core Logic Methods ---

TEST_F(CompareKMersTest, CompareLinePerfectMatch) {
    kmer_comparer->set_kmersize(4);
    // k-mer at 0: "ACGT". k-mer at 1: "CGTA".
    // Overlap: Suffix(k-1) of kmer 0 is "CGT"
    //          Prefix(k-1) of kmer 1 is "CGT"
    // Result: 1 (Match)
    EXPECT_EQ(kmer_comparer->compare_line(0, 1), 1);

    // k-mer at 6: "CGGT". k-mer at 7: "GGTA".
    // Overlap: Suffix "GGT" vs Prefix "GGT"
    // Result: 1 (Match)
    EXPECT_EQ(kmer_comparer->compare_line(6, 7), 1);
}

TEST_F(CompareKMersTest, CompareLineMismatch) {
    kmer_comparer->set_kmersize(4);
    // k-mer at 0: "ACGT". k-mer at 4: "AACG".
    // Overlap: Suffix "CGT" vs Prefix "AAC"
    // Result: 0 (Mismatch)
    EXPECT_EQ(kmer_comparer->compare_line(0, 4), 0);
}

TEST_F(CompareKMersTest, CompareLinePartialMismatch) {
    kmer_comparer->set_kmersize(4);
    // k-mer at 8: "GTAC". k-mer at 9: "TACG".
    // Overlap: Suffix "TAC" vs Prefix "TAC"
    // Result: 1 (Match)
    EXPECT_EQ(kmer_comparer->compare_line(8, 9), 1);

    // k-mer at 7: "GGTA". k-mer at 8: "GTAC".
    // Overlap: Suffix "GTA" vs Prefix "GTA"
    // Result: 1 (Match)
    EXPECT_EQ(kmer_comparer->compare_line(7, 8), 1);
}

TEST_F(CompareKMersTest, CompareLines) {
    // Tests compare_lines(ref). This compares k-mer at `ref`
    // with k-mers starting at `reads[0]`, `reads[1]`, `reads[2]`.
    // `reads` = {4, 8, 12}
    kmer_comparer->set_kmersize(4);

    // --- compare_lines(0) --- (compares k-mer at 0: "ACGT")
    // i=0: ref==i -> 1
    // i=1: compare_line(0, 8) ["ACGT" vs "GTAC"]. Overlap "CGT" vs "GTA" -> 0
    // i=2: compare_line(0, 12) ["ACGT" vs "CGGG"]. Overlap "CGT" vs "CGG" -> 0
    std::vector<size_t> expected0 = {1, 0, 0};
    EXPECT_EQ(kmer_comparer->compare_lines(0), expected0);

    // --- compare_lines(1) --- (compares k-mer at 1: "CGTA")
    // i=0: compare_line(1, 4) ["CGTA" vs "AACG"]. Overlap "GTA" vs "AAC" -> 0
    // i=1: ref==i -> 1
    // i=2: compare_line(1, 12) ["CGTA" vs "CGGG"]. Overlap "GTA" vs "CGG" -> 0
    std::vector<size_t> expected1 = {0, 1, 0};
    EXPECT_EQ(kmer_comparer->compare_lines(1), expected1);

    // --- compare_lines(4) --- (compares k-mer at 4: "AACG")
    // i=0: compare_line(4, 4) ["AACG" vs "AACG"]. ref(4) != i(0) -> compare_line -> 0
    //     (Note: compare_line(4,4) checks suffix "ACG" vs prefix "AAC" -> 0)
    // i=1: compare_line(4, 8) ["AACG" vs "GTAC"]. Overlap "ACG" vs "GTA" -> 0
    // i=2: compare_line(4, 12) ["AACG" vs "CGGG"]. Overlap "ACG" vs "CGG" -> 0
    std::vector<size_t> expected4 = {0, 0, 0};
    EXPECT_EQ(kmer_comparer->compare_lines(4), expected4);
}

TEST_F(CompareKMersTest, CompareAll) {
    kmer_comparer->set_kmersize(4);
    std::vector<std::vector<size_t>> results = kmer_comparer->compare_all();

    // Check dimensions
    ASSERT_EQ(results.size(), 3); // results.size() == reads.size()
    for (const auto& row : results) {
        ASSERT_EQ(row.size(), 3); // row.size() == reads.size()
    }

    // Spot-check results based on the logic from CompareLines test
    
    // Row 0: results of compare_lines(0)
    std::vector<size_t> expected0 = {1, 0, 0};
    EXPECT_EQ(results[0], expected0);

    // Row 1: results of compare_lines(1)
    std::vector<size_t> expected1 = {0, 1, 0};
    EXPECT_EQ(results[1], expected1);
    
    // Row 2: results of compare_lines(2)
    // i=0: compare_line(2, 4) ["GTAA" vs "AACG"]. Overlap "TAA" vs "AAC" -> 0
    // i=1: compare_line(2, 8) ["GTAA" vs "GTAC"]. Overlap "TAA" vs "GTA" -> 0
    // i=2: ref==i -> 1
    std::vector<size_t> expected2 = {0, 0, 1};
    EXPECT_EQ(results[2], expected2);
}