#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <iomanip>

struct Alignment {
    std::string aligned_s;
    std::string aligned_t;
    int num_gaps;
    int num_mismatches;
};

int max(int a, int b, int c) {
    return std::max({a, b, c});
}

void traceback(const std::vector<std::vector<int>>& M, const std::string& s, const std::string& t,
               std::string aligned_s, std::string aligned_t, int i, int j, int match, int mismatch, int gap,
               std::vector<Alignment>& alignments) {
    if (i == 0 && j == 0) {
        int num_gaps = std::count(aligned_s.begin(), aligned_s.end(), '-') + std::count(aligned_t.begin(), aligned_t.end(), '-');
        int num_mismatches = 0;
        for (size_t k = 0; k < aligned_s.length(); ++k) {
            if (aligned_s[k] != aligned_t[k] && aligned_s[k] != '-' && aligned_t[k] != '-') {
                num_mismatches++;
            }
        }
        alignments.push_back({aligned_s, aligned_t, num_gaps, num_mismatches});
        return;
    }

    if (i > 0 && j > 0 && M[i][j] == M[i-1][j-1] + (s[i-1] == t[j-1] ? match : mismatch)) {
        traceback(M, s, t, s[i-1] + aligned_s, t[j-1] + aligned_t, i-1, j-1, match, mismatch, gap, alignments);
    }
    if (i > 0 && M[i][j] == M[i-1][j] + gap) {
        traceback(M, s, t, s[i-1] + aligned_s, "-" + aligned_t, i-1, j, match, mismatch, gap, alignments);
    }
    if (j > 0 && M[i][j] == M[i][j-1] + gap) {
        traceback(M, s, t, "-" + aligned_s, t[j-1] + aligned_t, i, j-1, match, mismatch, gap, alignments);
    }
}

void needlemanWunsch(const std::string &s, const std::string &t, int match, int mismatch, int gap, const std::string& output_file, const std::string& optimal_file) {
    auto start = std::chrono::high_resolution_clock::now();  // Iniciar el cronómetro

    int n = s.length();
    int m = t.length();
    std::vector<std::vector<int>> M(n + 1, std::vector<int>(m + 1, 0));

    // Inicialización de la primera fila y columna
    for (int i = 0; i <= n; ++i) M[i][0] = i * gap;
    for (int j = 0; j <= m; ++j) M[0][j] = j * gap;

    // Relleno de la matriz
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int scoreDiagonal = M[i - 1][j - 1] + (s[i - 1] == t[j - 1] ? match : mismatch);
            int scoreUp = M[i - 1][j] + gap;
            int scoreLeft = M[i][j - 1] + gap;
            M[i][j] = max(scoreDiagonal, scoreUp, scoreLeft);
        }
    }

    // Vector para almacenar todos los alineamientos
    std::vector<Alignment> alignments;

    // Llamada a la función recursiva para obtener todas las soluciones óptimas
    traceback(M, s, t, "", "", n, m, match, mismatch, gap, alignments);

    // Encontrar la mejor solución según el criterio de menor número de gaps y mismatches
    Alignment best_alignment = alignments[0];
    for (const auto& alignment : alignments) {
        if (alignment.num_gaps < best_alignment.num_gaps ||
            (alignment.num_gaps == best_alignment.num_gaps && alignment.num_mismatches < best_alignment.num_mismatches)) {
            best_alignment = alignment;
        }
    }

    // Calcular el tiempo total
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Escribir el resultado completo en el archivo output_file
    std::ofstream outfile(output_file);
    outfile << "Matriz de Score:\n";
    for (const auto& row : M) {
        for (int val : row) {
            outfile << val << " ";
        }
        outfile << std::endl;
    }
    outfile << "\nScore optimo: " << M[n][m] << std::endl;

    outfile << "\nAlineamientos optimos:\n";
    for (const auto& alignment : alignments) {
        outfile << alignment.aligned_s << "\n" << alignment.aligned_t << "\n";
        outfile << "Gaps: " << alignment.num_gaps << ", Mismatches: " << alignment.num_mismatches << "\n\n";
    }

    // Escribir el mejor alineamiento en el archivo optimal_file
    std::ofstream optfile(optimal_file);
    optfile << "El alineamiento mas optimo es:\n";
    optfile << best_alignment.aligned_s << "\n" << best_alignment.aligned_t << "\n";
    optfile << "Gaps: " << best_alignment.num_gaps << ", Mismatches: " << best_alignment.num_mismatches << std::endl;
    optfile << std::fixed << std::setprecision(6);  // Fijar la precisión a 6 decimales
    optfile << "Tiempo de ejecución: " << duration.count() << " segundos" << std::endl;

    // Visualización de la Matriz de Puntos en output_file
    outfile << "\nMatriz de Puntos (Coincidencias):\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (s[i] == t[j]) {
                outfile << "* ";
            } else {
                outfile << ". ";
            }
        }
        outfile << std::endl;
    }
}

int main() {
    std::string bacteria = "AAAC";
    std::string sars_cov = "AAAC";
    std::string influenza = "AGC";

    int match = 1;
    int mismatch = -1;
    int gap = -2;

    // Probar combinaciones (2, 3) de las secuencias
    needlemanWunsch(bacteria, sars_cov, match, mismatch, gap, "C:/Users/harry/Documents/Clion_projects/test/bacteria_sars_cov.txt", "C:/Users/harry/Documents/Clion_projects/test/optimal_bac_scars.txt");
    needlemanWunsch(sars_cov, influenza, match, mismatch, gap, "C:/Users/harry/Documents/Clion_projects/test/sars_cov_influenza.txt", "C:/Users/harry/Documents/Clion_projects/test/optimal_scars_infl.txt");
    needlemanWunsch(bacteria, influenza, match, mismatch, gap, "C:/Users/harry/Documents/Clion_projects/test/bacteria_influenza.txt", "C:/Users/harry/Documents/Clion_projects/test/optimal_bac_infl.txt");

    return 0;
}