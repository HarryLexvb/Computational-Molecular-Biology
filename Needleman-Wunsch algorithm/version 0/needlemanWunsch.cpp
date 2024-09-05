#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

using namespace std;

// Función para calcular el máximo de tres valores
int max(int a, int b, int c) {
    return std::max({a, b, c});
}

void needlemanWunsch(const std::string &s, const std::string &t, int match, int mismatch, int gap) {
    int n = s.length();
    int m = t.length();
    vector<vector<int>> M(n + 1, vector<int>(m + 1, 0));

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

    // Mostrar la matriz de scores
    cout << "Matriz de Score:" << endl;
    for (const auto& row : M) {
        for (int val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    // Mostrar la matriz de puntos
    cout << "\nMatriz de Puntos (Coincidencias):" << endl;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (s[i - 1] == t[j - 1]) {
                cout << "* ";
            } else {
                cout << ". ";
            }
        }
        cout << endl;
    }

    // Backtracking para obtener los alineamientos óptimos
    string aligned_s = "", aligned_t = "";
    int i = n, j = m;

    while (i > 0 && j > 0) {
        int score = M[i][j];
        int scoreDiagonal = M[i - 1][j - 1];
        int scoreUp = M[i - 1][j];
        int scoreLeft = M[i][j - 1];

        if (score == scoreDiagonal + (s[i - 1] == t[j - 1] ? match : mismatch)) {
            aligned_s = s[i - 1] + aligned_s;
            aligned_t = t[j - 1] + aligned_t;
            i--; j--;
        } else if (score == scoreUp + gap) {
            aligned_s = s[i - 1] + aligned_s;
            aligned_t = "-" + aligned_t;
            i--;
        } else {
            aligned_s = "-" + aligned_s;
            aligned_t = t[j - 1] + aligned_t;
            j--;
        }
    }

    while (i > 0) {
        aligned_s = s[i - 1] + aligned_s;
        aligned_t = "-" + aligned_t;
        i--;
    }
    while (j > 0) {
        aligned_s = "-" + aligned_s;
        aligned_t = t[j - 1] + aligned_t;
        j--;
    }

    // Mostrar los resultados
    cout << "\nAlineamiento Optimo:" << endl;
    cout << aligned_s << endl;
    cout << aligned_t << endl;
    cout << "Score: " << M[n][m] << endl;

    // Guardar los resultados en un archivo
    ofstream outfile("alignment_output.txt");
    outfile << "Matriz de Score:\n";
    for (const auto& row : M) {
        for (int val : row) {
            outfile << val << " ";
        }
        outfile << endl;
    }

    outfile << "\nMatriz de Puntos (Coincidencias):\n";
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            if (s[i - 1] == t[j - 1]) {
                outfile << "* ";
            } else {
                outfile << ". ";
            }
        }
        outfile << endl;
    }

    outfile << "\nAlineamiento Optimo:\n";
    outfile << aligned_s << endl;
    outfile << aligned_t << endl;
    outfile << "Score: " << M[n][m] << endl;

    outfile.close();
}

int main() {
    std::string bacteria = "AAAC";
    std::string influenza = "AGC";
    int match = 1;
    int mismatch = -1;
    int gap = -2;

    needlemanWunsch(bacteria, influenza, match, mismatch, gap);

    return 0;
}