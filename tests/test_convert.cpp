#include "gtest/gtest.h"
#include "convert.h"
#include <fstream>
#include <string>
#include <vector>
#include <filesystem> // C++17 requis (utilisé dans votre main.cpp)

namespace fs = std::filesystem;

// --- Helper : Gestion des fichiers temporaires ---

class ConvertTest : public ::testing::Test {
protected:
    std::string temp_file = "temp_test.fasta";

    void SetUp() override {
        // S'assure que le fichier n'existe pas avant le test
        if (fs::exists(temp_file)) {
            fs::remove(temp_file);
        }
    }

    void TearDown() override {
        // Nettoyage après le test
        if (fs::exists(temp_file)) {
            fs::remove(temp_file);
        }
    }

    // Crée un fichier FASTA avec le contenu donné
    void createFasta(const std::string& content) {
        std::ofstream out(temp_file);
        out << content;
        out.close();
    }
};

// --- Tests ---

TEST_F(ConvertTest, FileNotFoundError) {
    Convert converter;
    // On essaie de lire un fichier qui n'existe pas (SetUp nettoie tout)
    EXPECT_THROW(converter.processFile("non_existent_ghost_file.fasta"), std::runtime_error);
}

TEST_F(ConvertTest, SingleReadSimple) {
    Convert converter;
    // >header
    // ACGT
    createFasta(">seq1\nACGT");

    converter.processFile(temp_file);

    // Vérification du BitVector
    // ACGT -> 4 nucléotides -> 8 bits
    const BitVector& bv = converter.getBitVector();
    EXPECT_EQ(bv.size(), 8);
    EXPECT_EQ(bv.readBitVector(), "ACGT");

    // Vérification des positions de fin (endPos)
    const std::vector<size_t>& ends = converter.getEndPos();
    ASSERT_EQ(ends.size(), 1);
    EXPECT_EQ(ends[0], 8); // La lecture finit au bit 8
}

TEST_F(ConvertTest, MultiLineSequence) {
    Convert converter;
    // FASTA permet de couper les séquences sur plusieurs lignes
    // >seq1
    // AC
    // GT
    createFasta(">seq1\nAC\nGT");

    converter.processFile(temp_file);

    const BitVector& bv = converter.getBitVector();
    EXPECT_EQ(bv.readBitVector(), "ACGT");

    const std::vector<size_t>& ends = converter.getEndPos();
    ASSERT_EQ(ends.size(), 1);
    EXPECT_EQ(ends[0], 8);
}

TEST_F(ConvertTest, MultipleReads) {
    Convert converter;
    // >r1 (AC -> 4 bits)
    // AC
    // >r2 (TG -> 4 bits)
    // TG
    // >r3 (A -> 2 bits)
    // A
    createFasta(">r1\nAC\n>r2\nTG\n>r3\nA");

    converter.processFile(temp_file);

    const BitVector& bv = converter.getBitVector();
    // Total : AC + TG + A = "ACGTA" (5 nucs -> 10 bits)
    EXPECT_EQ(bv.readBitVector(), "ACTGA"); // AC TG A
    EXPECT_EQ(bv.size(), 10);

    // Vérification endPos (cumulatif)
    const std::vector<size_t>& ends = converter.getEndPos();
    ASSERT_EQ(ends.size(), 3);

    EXPECT_EQ(ends[0], 4);  // Fin de r1 (AC) à 4 bits
    EXPECT_EQ(ends[1], 8);  // Fin de r2 (TG) à 4+4=8 bits
    EXPECT_EQ(ends[2], 10); // Fin de r3 (A) à 8+2=10 bits
}

TEST_F(ConvertTest, IgnoreWhitespaceAndEmptyLines) {
    Convert converter;
    // Test robustesse parsing : espaces dans la séquence, lignes vides
    // >seq1
    // A C
    //
    // G T
    createFasta(">seq1\nA C\n\n\nG T");

    converter.processFile(temp_file);

    EXPECT_EQ(converter.getBitVector().readBitVector(), "ACGT");
}

TEST_F(ConvertTest, ResetBetweenFiles) {
    Convert converter;

    // Passe 1 : Fichier A
    createFasta(">s1\nAA"); // 4 bits
    converter.processFile(temp_file);
    EXPECT_EQ(converter.getBitVector().size(), 4);
    EXPECT_EQ(converter.getEndPos().size(), 1);

    // Passe 2 : Fichier B (réutilisation de l'objet Convert)
    // L'objet doit se vider avant de traiter le nouveau fichier
    createFasta(">s2\nCC\n>s3\nGG"); // 4 + 4 = 8 bits
    converter.processFile(temp_file);

    const BitVector& bv = converter.getBitVector();
    const std::vector<size_t>& ends = converter.getEndPos();

    // Ne doit contenir QUE les données du Fichier B (CCGG)
    EXPECT_EQ(bv.size(), 8);
    EXPECT_EQ(bv.readBitVector(), "CCGG");

    ASSERT_EQ(ends.size(), 2);
    EXPECT_EQ(ends[0], 4);
    EXPECT_EQ(ends[1], 8);
}

TEST_F(ConvertTest, EmptyFile) {
    Convert converter;
    createFasta(""); // Fichier vide

    converter.processFile(temp_file);

    EXPECT_EQ(converter.getBitVector().size(), 0);
    EXPECT_EQ(converter.getEndPos().size(), 0);
}