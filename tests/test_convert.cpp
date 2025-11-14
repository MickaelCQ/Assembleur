#include "gtest/gtest.h"
#include "convert.h"
#include <fstream>
#include <string>
#include <vector>
#include <cstdio> // For std::remove

// Fixture for creating a temporary test file
class ConvertTest : public ::testing::Test {
protected:
    std::string test_filename = "temp_test_convert.fa";

    void SetUp() override {
        // Create a temporary FASTA file for testing
        std::ofstream test_file(test_filename);
        ASSERT_TRUE(test_file.is_open());
        test_file << ">read1\n";
        test_file << "ACGT\n";
        test_file << ">read2\n";
        test_file << "TT\n";
        test_file.close();
    }

    void TearDown() override {
        // Clean up the temporary file
        std::remove(test_filename.c_str());
    }
};

// Test case for a simple, known sequence
TEST_F(ConvertTest, ProcessesSimpleSequenceCorrectly) {
    Convert converter;
    converter.process_fasta_file(test_filename);

    const auto& bit_vec = converter.get_bit_vector();
    const auto& end_pos = converter.get_read_end_positions();

    // --- Check Read End Positions ---
    // We expect two reads
    ASSERT_EQ(end_pos.size(), 2);
    // Read 1: "ACGT" -> 4 nucleotides * 2 bits/nuc = 8 bits
    EXPECT_EQ(end_pos[0], 8);
    // Read 2: "TT" -> 2 nucleotides * 2 bits/nuc = 4 bits
    // Total bits = 8 (from read 1) + 4 (from read 2) = 12
    EXPECT_EQ(end_pos[1], 12);

    // --- Check Bit Vector Content ---
    ASSERT_EQ(bit_vec.size(), 12);

    // Expected:
    // Read 1: A=00, C=10, G=01, T=11 -> 00100111
    // Read 2: T=11, T=11             -> 1111
    // Total: 001001111111
    std::vector<bool> expected_vec = {
        // Read 1: ACGT
        false, false, // A (00)
        true,  false, // C (10)
        false, true,  // G (01)
        true,  true,  // T (11)
        // Read 2: TT
        true,  true,  // T (11)
        true,  true   // T (11)
    };

    EXPECT_EQ(bit_vec, expected_vec);
}

// Test case for loading one of the real FASTA files
TEST(ConvertLoadTest, LoadsRealFastaFile) {
    Convert converter;

    // Path relative to the build directory where 'run_tests' is executed
#ifndef TEST_DATA_DIR
    // Fallback if macro not defined
    std::string real_fasta_path = "Tests_Et_Ref/minia.contigs.fa";
#else
    std::string real_fasta_path = std::string(TEST_DATA_DIR) + "/minia.contigs.fa";
#endif

    // Check if file exists before trying to open, to give a better error
    std::ifstream f(real_fasta_path);
    ASSERT_TRUE(f.good()) << "Test file not found at: " << real_fasta_path
                          << ". Make sure you are running tests from the 'build' directory.";
    f.close();


    // Process the file
    ASSERT_NO_THROW(converter.process_fasta_file(real_fasta_path));

    const auto& bit_vec = converter.get_bit_vector();
    const auto& end_pos = converter.get_read_end_positions();

    // Basic sanity checks for a real file
    EXPECT_GT(bit_vec.size(), 0);
    EXPECT_GT(end_pos.size(), 0);

    // The size of the full bit vector should be equal to the last end position
    EXPECT_EQ(bit_vec.size(), end_pos.back());
}

