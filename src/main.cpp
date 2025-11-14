//
// Created by raphael on 11/13/25.
//

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream> // NOUVEAU: Pour l'exportation de fichiers (std::ofstream)
#include <bitset>



#include "convert.h"   // Notre classe pour lire le FASTA
#include "compare.h" // Notre classe pour comparer les k-mers

/**
 * @brief Imprime la matrice de comparaison de manière lisible.
 */
void print_matrix(const std::vector<std::vector<size_t>>& matrix) {
    if (matrix.empty()) {
        std::cout << "Matrice de comparaison vide." << std::endl;
        return;
    }

    // Imprimer l'en-tête (index des colonnes)
    std::cout << "      ";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << "[" << j << "]\t";
    }
    std::cout << std::endl << "------";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << "----";
    }
    std::cout << std::endl;

    // Imprimer les lignes
    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << "[" << i << "] | ";
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

// NOUVEAU: Fonction pour exporter la matrice vers un fichier TSV
/**
 * @brief Exporte une matrice 2D vers un fichier TSV (Tab-Separated Values).
 * @param matrix La matrice de données à exporter.
 * @param filename Le nom du fichier de sortie.
 * @throws std::runtime_error si le fichier ne peut pas être ouvert.
 */
void export_matrix_to_tsv(const std::vector<std::vector<size_t>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Erreur: Impossible d'ouvrir le fichier de sortie " + filename);
    }

    // Pas d'en-tête pour une matrice simple, écriture directe des données.
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j];
            // Ajouter une tabulation sauf pour le dernier élément
            if (j < matrix[i].size() - 1) {
                file << "\t";
            }
        }
        // Nouvelle ligne à la fin de chaque ligne de la matrice
        file << "\n";
    }

    file.close();
}


int main(int argc, char* argv[]) {
    // --- 1. Vérifier les arguments ---
    // MODIFIÉ: Attend 3 arguments (exécutable, entrée, sortie)
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fichier_fasta.fa> <fichier_sortie.tsv>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string output_filename = argv[2]; // NOUVEAU: Nom du fichier de sortie
    const size_t KMER_SIZE = 31;

    // --- 2. Utiliser Convert pour traiter le fichier ---
    Convert converter;
    try {
        std::cout << "Traitement du fichier FASTA: " << filename << "..." << std::endl;
        converter.process_fasta_file(filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "Erreur lors du traitement du fichier: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Conversion terminee." << std::endl;

    std::vector<bool> bit_vector = std::move(converter.get_bit_vector());
    std::vector<size_t> read_ends = std::move(converter.get_read_end_positions());

    // --- 3. Utiliser CompareKMers avec les données ---
    std::cout << "Initialisation de la comparaison (k=" << KMER_SIZE << ")..." << std::endl;

    if (read_ends.empty()) {
        std::cout << "Aucune lecture trouvee dans le fichier." << std::endl;
        return 0;
    }

    CompareKMers comparer(std::move(bit_vector), std::move(read_ends), KMER_SIZE);

    // --- 4. Afficher les résultats ---
    size_t total_reads = comparer.get_n_reads();
    size_t total_kmers = comparer.get_all_n_kmers();

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Statistiques :" << std::endl;
    std::cout << "  Nombre total de lectures : " << total_reads << std::endl;
    std::cout << "  Nombre total de k-mers : " << total_kmers << std::endl;
    std::cout << "-----------------------------------" << std::endl;

    // MODIFIÉ: Logique de génération et d'exportation
    if (total_reads > 0) {
        std::cout << "Generation de la matrice de comparaison..." << std::endl;
        std::vector<std::vector<size_t>> comparison_matrix = comparer.compare_all();

        // NOUVEAU: Exporter la matrice vers le fichier TSV
        std::cout << "Exportation de la matrice vers " << output_filename << "..." << std::endl;
        try {
            export_matrix_to_tsv(comparison_matrix, output_filename);
            std::cout << "Exportation terminee." << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << "Erreur lors de l'exportation: " << e.what() << std::endl;
        }

        // NOUVEAU: N'imprimer la matrice dans la console que si elle est petite
        if (total_reads > 50) {
            std::cout << "Matrice trop grande pour l'affichage console (>50 lectures)." << std::endl;
            std::cout << "Les resultats ont ete enregistres dans " << output_filename << std::endl;
        } else {
            std::cout << "Matrice de resultats (Lecture[i] vs Lecture[j]):" << std::endl;
            print_matrix(comparison_matrix);
        }

    } else {
        std::cout << "Aucune lecture a comparer, aucun fichier TSV genere." << std::endl;
    }

    return 0;
}