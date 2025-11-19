//
// Created by raphael on 11/13/25.
// Additions to understanding and reflections on our programming choices 14/11/2025.

#include <iostream> // Indispensable pr communiquer avec l'utilisateur.
#include <string>
#include <stdexcept> // pour avoir des erreurs "propres", les exceptions standards...
#include <vector>   // La structure de base retrouvée dans le stockage de nos données. (bit vector \& matrice ..)
#include <fstream> // NOUVEAU: Pour l'exportation de fichiers (std::ofstream)
#include <bitset>


#include "graphdbj.h"
#include "convert.h"   // Notre classe pour lire le FASTA
#include "compare.h" // Notre classe pour comparator les k-mers (issus du FASTA).

/**
 * @brief Imprime la matrice de comparaison de manière lisible.
 * La fonction sera utilisée pour des matrices de taille raisonnable, considérant que l'objectif étant un affichage de confivialité mais pas d'inonder la console avec des milliers de lignes de code.
 * @param matrix Matrice qui contient nos valeurs de comparaison.
 * 
 * @complexity Temps: O(n*m) affichage de la matrice.
 * @complexity Espace: O(1) a priori pas d'allocation mémoire notable pr l'affichage dans le terminal. 	
 */
void print_matrix(const std::vector<std::vector<size_t>>& matrix) {
    if (matrix.empty()) {
        std::cout << "Matrice de comparaison vide." << std::endl;
        return;
    }

    // Imprimer l'en-tête (index des colonnes) : nous avons choisi cet affichage pour faciliter la lecture.
    std::cout << "      ";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << "[" << j << "]\t"; // Avec des tab on garantit un alignement convenable. 
    }
    std::cout << std::endl << "------";
    for (size_t j = 0; j < matrix[0].size(); ++j) {
        std::cout << "----";
    }
    std::cout << std::endl;

    // Imprimer les lignes
    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << "[" << i << "] | "; // Rappel de l'indice du read.
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            std::cout << matrix[i][j] << "\t"; // La valeur affichée représente notre nombre de kmers paratagés.
        }
        std::cout << std::endl;
    }
}

// NOUVEAU: Fonction pour exporter la matrice vers un fichier TSV
/**
 * @brief Exporte une matrice 2D vers un fichier TSV (Tab-Separated Values). 
 * Nous avons choisi un format tsv car c'est à la fois simple à parser (dans R, Python...), on évite les pièges des CSV. Puis c'est compatible avec une majorité des pipelines courants.
 * Si on s'inscrit dans une logique d'intégration de l'outil (au dela du cadre pédagogique).
 * @param matrix La matrice de données à exporter.
 * @param filename Le path du fichier de sortie.
 * @throws std::runtime_error si l'ouverture du fichier échoue.
 * @complexity Temps: 0(n*m) , chaque cellule écrite une fois.
 * @complexity Espace:O(1), vu qu'on à une écriture en flux, à priori pas d'usage énorme du buffer...(à vérifier par loik).
 */
void export_matrix_to_tsv(const std::vector<std::vector<size_t>>& matrix, const std::string& filename) {
    std::ofstream file(filename);
    //Vérification(s) d'usage(s) pour éviter les comportements silencieux difficiles à diagnostiquer. 
    if (!file.is_open()) {
        throw std::runtime_error("Erreur: Impossible d'ouvrir le fichier de sortie " + filename);
    }
	// MICKAEL : A mon sens faut qu'on en rajoute ...
    
    // Pas d'en-tête pour une matrice simple, écriture directe des données, forme "Val\tVal\tVal.."
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[i].size(); ++j) {
            file << matrix[i][j];
            // Ajouter une tabulation sauf pour le dernier élément
            if (j < matrix[i].size() - 1) {
                file << "\t";
            }
        }
        file << "\n";
    }

    file.close();// Fermeture explicite toujours plus propre pour les gros fichiers. 
}



/**@brief

 *complexity Temps: (A compléter à la fin) 
 *complexity Temps: (A compléter à la fin)

*/
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
        converter.processFile(filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "Erreur lors du traitement du fichier: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Conversion terminee." << std::endl;

    // --- Construction du Graphe de De Bruijn ---
    std::cout << "Construction du graphe de Bruijn..." << std::endl;

    // Création de l'objet graphe (convertit automatiquement les bitvectors en noeuds)
    GraphDBJ graph(converter, KMER_SIZE);

    // Récupération du vecteur de pointeurs de noeuds
    std::vector<Noeud*> tousLesNoeuds = graph.getNodes();

    std::cout << "Graphe construit." << std::endl;
    std::cout << "Nombre de noeuds uniques (k-1 mers) : " << tousLesNoeuds.size() << std::endl;

    // Exemple d'affichage pour un petit graphe
    if (tousLesNoeuds.size() < 20) {
        for (Noeud* n : tousLesNoeuds) {
            std::cout << "Noeud " << n->p << " a " << n->c.size() << " enfants." << std::endl;
        }
    }

    const BitVector& bitVector_ref = converter.get_bitVector();
    std::vector<size_t> read_ends = std::move(converter.get_read_end_positions());

    // --- 3. Utiliser CompareKMers avec les données ---
    std::cout << "Initialisation de la comparaison (k=" << KMER_SIZE << ")..." << std::endl;

    if (read_ends.empty()) {
        std::cout << "Aucune lecture trouvee dans le fichier." << std::endl;
        return 0;
    }

    CompareKMers comparator(bitVector_ref, std::move(read_ends), KMER_SIZE);

    // --- 4. Afficher les résultats ---
    size_t total_reads = comparator.get_nReads();
    size_t total_kmers = comparator.get_all_nKmers();

    std::cout << "-----------------------------------" << std::endl;
    std::cout << "Statistiques :" << std::endl;
    std::cout << "  Nombre total de lectures : " << total_reads << std::endl;
    std::cout << "  Nombre total de k-mers : " << total_kmers << std::endl;
    std::cout << "-----------------------------------" << std::endl;

    // MODIFIÉ: Logique de génération et d'exportation
    if (total_reads > 0) {
        // --- Construction du Graphe de De Bruijn ---
        std::cout << "Construction du graphe de Bruijn..." << std::endl;
        GraphDBJ graph(converter, KMER_SIZE);

        // --- SIMPLIFICATION DU GRAPHE ---
        // 1. Elagage des pointes
        graph.removeTips(KMER_SIZE - 1);

        // 2. Resolution des bulles (utilisez la version avancée si vous l'avez implémentée)
        graph.resolveBubbles();

        // --- GENERATION DE L'ASSEMBLAGE FINAL ---
        std::cout << "Generation des Contigs (Extension par couverture)..." << std::endl;

        // On utilise generateContigs AU LIEU DE generateUnitigs
        std::vector<std::string> contigs = graph.generateContigs();

        std::cout << "Nombre de contigs generes : " << contigs.size() << std::endl;

        // Export
        std::string contig_filename = output_filename + ".contigs.fasta";
        std::ofstream out_contigs(contig_filename);

        // On ne garde souvent que les contigs d'une certaine taille (ex: > 100bp)
        int contigs_exported = 0;
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (contigs[i].length() >= (size_t)KMER_SIZE * 2) { // Filtre simple
                out_contigs << ">contig_" << i << "_len_" << contigs[i].length() << "\n";
                out_contigs << contigs[i] << "\n";
                contigs_exported++;
            }
        }
        out_contigs.close();
        std::cout << "Contigs exportes : " << contigs_exported << " (filtres par taille)" << std::endl;

    } else {
        std::cout << "Aucune lecture a comparator, aucun fichier TSV genere." << std::endl;
    }

    return 0;
}
