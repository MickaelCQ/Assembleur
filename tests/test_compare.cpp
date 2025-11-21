#include "gtest/gtest.h"
#include "compare.h"
#include "bitvector.h"
#include <vector>
#include <string>
#include <stdexcept>

// --- Helper : Construction facile ---
// Construit un BitVector et le vecteur endPos à partir d'une liste de chaînes (lectures)
class CompareTest : public ::testing::Test {
protected:
    BitVector bv;
    std::vector<size_t> reads;

    // Prépare les données pour le test
    void setupReads(const std::vector<std::string>& sequences) {
        bv = BitVector();
        reads.clear();

        for (const std::string& seq : sequences) {
            for (char c : seq) {
                bv.addCha(c);
            }
            reads.push_back(bv.size());
        }
    }
};

// --- Tests des Getters / Setters / Calculs de base ---

TEST_F(CompareTest, KmerCountingLogic) {
    // Setup : 3 lectures
    // R1: "ACGT" (4bp)
    // R2: "A"    (1bp)
    // R3: "AAAAA"(5bp)
    setupReads({"ACGT", "A", "AAAAA"});

    CompareKMers comp(bv, reads, 3); // k=3

    // R1 (4bp) -> 4 - 3 + 1 = 2 k-mers
    EXPECT_EQ(comp.get_nKmers(0), 2);

    // R2 (1bp) -> 1 < 3 -> 0 k-mers
    EXPECT_EQ(comp.get_nKmers(1), 0);

    // R3 (5bp) -> 5 - 3 + 1 = 3 k-mers
    EXPECT_EQ(comp.get_nKmers(2), 3);

    // Total : 2 + 0 + 3 = 5
    EXPECT_EQ(comp.get_all_nKmers(), 5);
}

TEST_F(CompareTest, ReadEndPositions) {
    setupReads({"AC", "GT"}); // 2bp (4 bits), 2bp (4 bits)
    CompareKMers comp(bv, reads);

    EXPECT_EQ(comp.get_nReads(), 2);
    EXPECT_EQ(comp.get_read_end_pos(0), 4);       // Fin R1
    EXPECT_EQ(comp.get_read_end_pos(1), 8);       // Fin R2 (4+4)
    EXPECT_THROW(comp.get_read_end_pos(2), std::out_of_range);
}

// --- Tests de la Logique de Comparaison (Overlap) ---

/**
 * IMPORTANT : Analyse de compare_line dans compare.cpp
 * La boucle compare :
 * bit_idx1 + (i + 1) * 2   VS   bit_idx2 + i * 2
 * * Cela signifie qu'on compare le SUFFIXE du k-mer 1 (commençant à l'index 1)
 * avec le PREFIXE du k-mer 2 (commençant à l'index 0).
 * C'est une détection de chevauchement (De Bruijn Edge), pas une égalité stricte.
 */
TEST_F(CompareTest, CompareLine_OverlapDetection) {
    // k=3. On compare k-1 = 2 nucléotides.
    // R0: "ACG" (Suffixe k-1 : CG)
    // R1: "CGT" (Prefixe k-1 : CG)
    // R2: "TGC"
    setupReads({"ACG", "CGT", "TGC"});
    CompareKMers comp(bv, reads, 3);

    size_t posR0 = 0;
    size_t posR1 = 6;  // 3 nucs * 2 bits = 6
    size_t posR2 = 12; // 6 + 6 = 12

    // Cas 1 : Match (ACG -> CGT)
    // Le suffixe de R0 ("CG") == Préfixe de R1 ("CG")
    EXPECT_EQ(comp.compare_line(posR0, posR1), 1);

    // Cas 2 : Mismatch Inverse (CGT -> ACG)
    // Le suffixe de R1 ("GT") != Préfixe de R0 ("AC")
    EXPECT_EQ(comp.compare_line(posR1, posR0), 0);

    // Cas 3 : Chaine (CGT -> TGC)
    // Suffixe R1 ("GT") == Préfixe R2 ("TG") ? NON.
    // Suffixe R1 est "GT". Préfixe R2 est "TG". Ça ne matche pas.
    // "CGT" -> "GT". "TGC" -> "TG".
    EXPECT_EQ(comp.compare_line(posR1, posR2), 0);

    // Cas 4 : Vrai Match pour R1->R2 si on change les données
    // Si R1="CGT" et Rnew="GTA"
    // Suffixe R1 "GT", Prefixe Rnew "GT" -> Match.
}

TEST_F(CompareTest, CompareLine_ExactLogic) {
    // Vérification plus fine bit à bit
    // R1: AAAA (k=3 -> AAA). Suffixe AA
    // R2: AAAT (k=3 -> AAA). Préfixe AA
    setupReads({"AAAA", "AAAT"});
    CompareKMers comp(bv, reads, 3);

    // R1 start=0. R2 start=8 (4 nucs * 2 bits).
    // Overlap : R1(AAA) -> R2(AAA) ?
    // Suffixe R1 (AA) == Prefixe R2 (AA). OUI.
    EXPECT_EQ(comp.compare_line(0, 8), 1);
}

// --- Tests de Compare Lines (1 vs Tous) ---

TEST_F(CompareTest, CompareLines_VectorResult) {
    // k=3 (Overlap k-1 = 2)
    // R0: "ACG" -> Suffix "CG"
    // R1: "CGT" -> Prefix "CG" (Match avec R0)
    // R2: "CGA" -> Prefix "CG" (Match avec R0)
    // R3: "TGC" -> Prefix "TG" (Pas de match avec R0)
    setupReads({"ACG", "CGT", "CGA", "TGC"});
    CompareKMers comp(bv, reads, 3);

    // On compare R0 contre tout le monde
    std::vector<size_t> results = comp.compare_lines(0);

    ASSERT_EQ(results.size(), 4);

    EXPECT_EQ(results[0], 1); // R0 vs R0 (Auto-comparaison est forcée à 1 dans le code)
    EXPECT_EQ(results[1], 1); // R0 -> R1 (CG == CG) : OUI
    EXPECT_EQ(results[2], 1); // R0 -> R2 (CG == CG) : OUI
    EXPECT_EQ(results[3], 0); // R0 -> R3 (CG != TG) : NON
}

// --- Tests de Compare All (Matrice N x N) ---

TEST_F(CompareTest, CompareAll_Matrix) {
    // k=2 (Overlap k-1 = 1 nucléotide)
    // R0: "AC" (Suffixe 'C')
    // R1: "CG" (Préfixe 'C') -> R0 connecte à R1
    // R2: "GT" (Préfixe 'G') -> R1 connecte à R2
    setupReads({"AC", "CG", "GT"});
    CompareKMers comp(bv, reads, 2);

    auto matrix = comp.compare_all();

    /* Matrice attendue :
       Src \ Dest | R0(AC) | R1(CG) | R2(GT)
       -----------+--------+--------+-------
       R0 (Suf C) |   1    |   1    |   0   (C==A?No, C==C?Yes, C==G?No)
       R1 (Suf G) |   0    |   1    |   1   (G==A?No, G==C?No, G==G?Yes)
       R2 (Suf T) |   0    |   0    |   1

       Note: La diagonale est toujours 1 par définition dans le code.
    */

    // Ligne 0 (Source R0)
    EXPECT_EQ(matrix[0][0], 1);
    EXPECT_EQ(matrix[0][1], 1); // Match
    EXPECT_EQ(matrix[0][2], 0);

    // Ligne 1 (Source R1)
    EXPECT_EQ(matrix[1][0], 0);
    EXPECT_EQ(matrix[1][1], 1);
    EXPECT_EQ(matrix[1][2], 1); // Match

    // Ligne 2 (Source R2)
    EXPECT_EQ(matrix[2][0], 0);
    EXPECT_EQ(matrix[2][1], 0);
    EXPECT_EQ(matrix[2][2], 1);
}