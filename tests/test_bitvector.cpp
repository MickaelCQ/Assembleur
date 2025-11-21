#include "gtest/gtest.h"
#include "bitvector.h"
#include <string>
#include <vector>
#include <stdexcept>

// --- Tests de base (Construction et Accesseurs) ---

TEST(BitVectorTest, ConstructorAndDefaults) {
    BitVector bv;
    
    // Par défaut, taille 0 et sizeElement 2 (pour l'ADN)
    EXPECT_EQ(bv.size(), 0);
    EXPECT_EQ(bv.getSizeElement(), 2);
    
    // Constructeur avec paramètres
    BitVector bv2(4); 
    EXPECT_EQ(bv2.getSizeElement(), 4);
}

TEST(BitVectorTest, SetSizeElement) {
    BitVector bv;
    bv.setSizeElement(1);
    EXPECT_EQ(bv.getSizeElement(), 1);

    // Test invalid argument
    EXPECT_THROW(bv.setSizeElement(65), std::invalid_argument); // > 64
    EXPECT_THROW(bv.setSizeElement(0), std::invalid_argument);
}

// --- Tests de manipulation de bits bruts ---

TEST(BitVectorTest, PushBackAndTest) {
    BitVector bv;
    
    bv.push_back(true);  // idx 0
    bv.push_back(false); // idx 1
    bv.push_back(true);  // idx 2

    EXPECT_EQ(bv.size(), 3);
    
    EXPECT_TRUE(bv.test(0));
    EXPECT_FALSE(bv.test(1));
    EXPECT_TRUE(bv.test(2));
    
    // Test opérateur []
    EXPECT_TRUE(bv[0]);
    EXPECT_FALSE(bv[1]);
    
    // Test hors limites (fail-safe renvoie false)
    EXPECT_FALSE(bv.test(999));
}

TEST(BitVectorTest, ToVector) {
    BitVector bv;
    bv.push_back(true);
    bv.push_back(false);
    bv.push_back(true);

    std::vector<bool> vec = bv.to_vector();
    ASSERT_EQ(vec.size(), 3);
    EXPECT_EQ(vec[0], true);
    EXPECT_EQ(vec[1], false);
    EXPECT_EQ(vec[2], true);
}

// --- Tests de gestion mémoire (Block Boundaries) ---

TEST(BitVectorTest, BlockBoundaryCrossing) {
    // BitVector stocke par blocs de 64 bits.
    // On va remplir 64 bits, puis ajouter le 65ème pour forcer la création d'un nouveau bloc.
    BitVector bv;
    size_t boundary = 64;

    // Remplir le premier bloc avec des 1
    for (size_t i = 0; i < boundary; ++i) {
        bv.push_back(true);
    }

    // Ajouter un bit au bloc suivant (idx 64)
    bv.push_back(false); 
    // Ajouter un bit au bloc suivant (idx 65)
    bv.push_back(true);

    EXPECT_EQ(bv.size(), 66);
    
    // Vérifier la fin du bloc 1
    EXPECT_TRUE(bv.test(63));
    // Vérifier le début du bloc 2
    EXPECT_FALSE(bv.test(64));
    EXPECT_TRUE(bv.test(65));
}

TEST(BitVectorTest, ReserveAndClear) {
    BitVector bv;
    bv.reserve(1000); // Devrait allouer ~16 blocs de 64 bits
    
    for(int i=0; i<100; ++i) bv.push_back(true);
    EXPECT_EQ(bv.size(), 100);

    bv.clear();
    EXPECT_EQ(bv.size(), 0);
    
    // Vérifie qu'on peut réutiliser après un clear
    bv.push_back(false);
    EXPECT_EQ(bv.size(), 1);
    EXPECT_FALSE(bv[0]);
}

// --- Tests Spécifiques ADN (Nucléotides) ---

TEST(BitVectorTest, AddChaEncoding) {
    BitVector bv;
    // Encodage défini dans BitVector::addCha :
    // A : 00
    // C : 10 (Attention : push_back(true) puis push_back(false)) -> Bit poids fort stocké en premier ?
    // G : 01
    // T : 11
    
    bv.addCha('A'); 
    bv.addCha('C');
    bv.addCha('G');
    bv.addCha('T');

    // Taille attendue : 4 nucléotides * 2 bits = 8 bits
    EXPECT_EQ(bv.size(), 8);

    // Vérification binaire manuelle
    // A (00)
    EXPECT_FALSE(bv[0]); EXPECT_FALSE(bv[1]);
    // C (10) -> Premier bit ajouté est 1, deuxième est 0
    EXPECT_TRUE(bv[2]);  EXPECT_FALSE(bv[3]);
    // G (01)
    EXPECT_FALSE(bv[4]); EXPECT_TRUE(bv[5]);
    // T (11)
    EXPECT_TRUE(bv[6]);  EXPECT_TRUE(bv[7]);
}

TEST(BitVectorTest, InvalidCharThrow) {
    BitVector bv;
    EXPECT_THROW(bv.addCha('Z'), std::invalid_argument);
    EXPECT_NO_THROW(bv.addCha('a')); // minuscule devrait marcher (toupper utilisé)
}

TEST(BitVectorTest, ReadBitVectorRoundTrip) {
    BitVector bv;
    std::string input = "ACGTACGT";
    for (char c : input) {
        bv.addCha(c);
    }

    std::string output = bv.readBitVector();
    EXPECT_EQ(input, output);
}

TEST(BitVectorTest, GetNucleotideAt) {
    BitVector bv;
    // A=00 (0), C=10 (2), G=01 (1), T=11 (3)
    // Note : getNucleotideAt recrée l'entier : 
    // if (b1) val |= 2; if (b2) val |= 1;
    
    bv.addCha('A');
    bv.addCha('C');
    bv.addCha('G');
    bv.addCha('T');

    EXPECT_EQ(bv.getNucleotideAt(0), 0); // A
    EXPECT_EQ(bv.getNucleotideAt(1), 2); // C (Car encoded 10 -> bit 1 set -> val 2)
    EXPECT_EQ(bv.getNucleotideAt(2), 1); // G
    EXPECT_EQ(bv.getNucleotideAt(3), 3); // T
}

// --- Tests Fonctionnalités Avancées (Append, Resize) ---

TEST(BitVectorTest, Append) {
    BitVector bv1;
    bv1.addCha('A'); // 00

    BitVector bv2;
    bv2.addCha('C'); // 10
    bv2.addCha('T'); // 11

    // Append bv2 à la fin de bv1
    bv1.append(bv2);

    EXPECT_EQ(bv1.size(), 6); // 2 + 4 bits
    EXPECT_EQ(bv1.readBitVector(), "ACT");
}

TEST(BitVectorTest, AppendWithSkip) {
    BitVector bv1;
    bv1.addCha('A');

    BitVector bv2;
    // "CGT"
    bv2.addCha('C'); 
    bv2.addCha('G');
    bv2.addCha('T');

    // On veut ajouter bv2 à bv1, mais en sautant le premier nucléotide de bv2 (C)
    // Donc on ajoute "GT"
    bv1.append(bv2, 1); 

    EXPECT_EQ(bv1.readBitVector(), "AGT");
}

TEST(BitVectorTest, ResizeShrink) {
    BitVector bv;
    bv.addCha('A'); // 0,1
    bv.addCha('C'); // 2,3
    bv.addCha('G'); // 4,5
    bv.addCha('T'); // 6,7
    
    EXPECT_EQ(bv.size(), 8);

    // On garde seulement les 4 premiers bits (AC)
    bv.resize(4);
    
    EXPECT_EQ(bv.size(), 4);
    EXPECT_EQ(bv.readBitVector(), "AC");
    
    // Resize ne doit pas augmenter la taille (sécurité définie dans le .h)
    bv.resize(10);
    EXPECT_EQ(bv.size(), 4);
}