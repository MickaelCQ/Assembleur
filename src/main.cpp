//
// Created by raphael on 11/13/25.
// Commentaires ajoutés et traduits pour expliquer la logique globale.

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream> // Pour l'exportation de fichiers
#include <bitset>

#include "graphdbj.h"
#include "convert.h"
#include "compare.h"

/**
 * @brief Affiche la matrice de comparaison dans la console.
 * Utilisé principalement pour le débogage sur de petits jeux de données.
 */
void print_matrix(const std::vector<std::vector<size_t>>& matrix) {
    if (matrix.empty()) {
        std::cout << "Matrice de comparaison vide." << std::endl;
        return;
    }

    // En-tête des colonnes
    std::cout << "      ";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << "[" << j << "]\t";
    }
    std::cout << std::endl << "------";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << "----";
    }
    std::cout << std::endl;

    // Lignes
    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << "[" << i << "] | ";
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

/**
 * @brief Exporte une matrice vers un fichier TSV (Tab-Separated Values).
 * Format compatible avec Excel, R, Python (Pandas).
 */
void export_matrix_to_tsv(const std::vector<std::vector<size_t>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Erreur: Impossible d'ouvrir le fichier de sortie " + filename);
    }

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j];
            if (j < matrix[i].size() - 1) {
                file << "\t";
            }
        }
        file << "\n";
    }
    file.close();
}

int main(int argc, char* argv[]) {
    // --- 1. Validation des arguments CLI ---
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <fichier_fasta.fa> <fichier_sortie_prefixe>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::string output_filename = argv[2];
    const size_t KMER_SIZE = 31; // Taille fixe des k-mers pour cet exemple

    // --- 2. Lecture et Conversion du Fichier ---
    Convert converter;
    try {
        std::cout << "Traitement du fichier FASTA: " << filename << "..." << std::endl;
        converter.processFile(filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "Erreur lors du traitement du fichier: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Conversion terminee." << std::endl;

    // --- 3. Analyse via CompareKMers (Statistiques) ---
    const BitVector& bitVector_ref = converter.get_bitVector();
    // std::move transfère la propriété du vecteur (évite une copie coûteuse)
    std::vector<size_t> read_ends = converter.get_read_end_positions();

    if (read_ends.empty()) {
        std::cout << "Aucune lecture trouvee dans le fichier." << std::endl;
        return 0;
    }

    CompareKMers comparator(bitVector_ref, read_ends, KMER_SIZE);

    size_t total_reads = comparator.get_nReads();
    size_t total_kmers = comparator.get_all_nKmers();

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Statistiques :" << std::endl;
    std::cout << "  Nombre total de lectures : " << total_reads << std::endl;
    std::cout << "  Nombre total de k-mers theoriques : " << total_kmers << std::endl;
    std::cout << "-----------------------------------" << std::endl;

    if (total_reads > 0) {
        // --- 4. Construction du Graphe de De Bruijn ---
        std::cout << "Construction du graphe de Bruijn..." << std::endl;

        // On reconstruit le converter car 'read_ends' a été déplacé (move) plus haut.
        // NOTE: Dans un code optimisé, on éviterait de relire le fichier ou de déplacer les vecteurs si tôt.
        // Ici, nous le faisons pour la clarté du flux logique GraphDBJ.
        converter.processFile(filename);
        GraphDBJ graph(converter, KMER_SIZE);

        auto nodes = graph.getNodes();
        std::cout << "Graphe construit. Noeuds uniques (k-1 mers) : " << nodes.size() << std::endl;

        // --- 5. Simplification du Graphe (Boucle Itérative) ---
        std::cout << "--- Debut de la simplification iterative ---" << std::endl;

        bool changed = true;
        int max_iterations = 10; // Sécurité pour éviter boucle infinie
        int iter = 0;

        while (changed && iter < max_iterations) {
            changed = false;
            iter++;
            std::cout << "ITERATION " << iter << " :" << std::endl;

            // 1. On écrase les bulles d'abord (pour recoller les chemins)
            if (graph.resolveBubbles() > 0) {
                changed = true;
            }

            // 2. On coupe les pointes ensuite (qui ont pu apparaitre après fusion de bulles)
            if (graph.clipTips() > 0) {
                changed = true;
            }
        }

        std::cout << "Export du graphe en GFA..." << std::endl;
        graph.exportToGFA(output_filename + ".gfa");

        // --- 6. Génération et Export des Contigs ---
        std::cout << "Generation des Contigs (Extension par couverture)..." << std::endl;

        std::vector<std::string> contigs = graph.generateContigs();
        std::cout << "Nombre de contigs bruts : " << contigs.size() << std::endl;

        // --- ETAPE DE CONSENSUS ---
        // On demande un overlap minimum de 'k' nucléotides pour être sûr
        contigs = GraphDBJ::mergeContigs(contigs, KMER_SIZE);

        std::cout << "Nombre de contigs apres fusion et nettoyage : " << contigs.size() << std::endl;

        std::string contig_filename = output_filename + ".contigs.fasta";
        std::ofstream out_contigs(contig_filename);

        int contigs_exported = 0;
        for (size_t i = 0; i < contigs.size(); ++i) {
            // Filtre : on ne garde que les contigs ayant une taille raisonnable (ex: > 2*k)
            if (contigs[i].length() >= (size_t)KMER_SIZE * 2) {
                out_contigs << ">contig_" << i << "_len_" << contigs[i].length() << "\n";
                out_contigs << contigs[i] << "\n";
                contigs_exported++;
            }
        }
        out_contigs.close();
        std::cout << "Contigs exportes : " << contigs_exported << " (filtres par taille > " << KMER_SIZE*2 << ")" << std::endl;

    } else {
        std::cout << "Aucune lecture a traiter." << std::endl;
    }

    return 0;
}