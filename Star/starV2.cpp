#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <fstream>

using namespace std;

const int GAP_PENALTY = -1;
const int MATCH_SCORE = 2;
const int MISMATCH_SCORE = -1;

// Function to calculate the score of two characters
int score(char a, char b) {
    if (a == b) return MATCH_SCORE;
    if (a == '-' || b == '-') return GAP_PENALTY;
    return MISMATCH_SCORE;
}

// Function to perform pairwise alignment of two sequences
pair<int, pair<string, string>> pairwiseAlignment(const string& s1, const string& s2) {
    int m = s1.length(), n = s2.length();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    // Initialize the first row and column
    for (int i = 0; i <= m; i++) dp[i][0] = i * GAP_PENALTY;
    for (int j = 0; j <= n; j++) dp[0][j] = j * GAP_PENALTY;

    // Fill the dp matrix
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            dp[i][j] = max({dp[i-1][j-1] + score(s1[i-1], s2[j-1]),
                            dp[i-1][j] + GAP_PENALTY,
                            dp[i][j-1] + GAP_PENALTY});
        }
    }

    // Reconstruct the alignment
    string aligned1 = "", aligned2 = "";
    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i-1][j-1] + score(s1[i-1], s2[j-1])) {
            aligned1 = s1[i-1] + aligned1;
            aligned2 = s2[j-1] + aligned2;
            i--; j--;
        } else if (i > 0 && dp[i][j] == dp[i-1][j] + GAP_PENALTY) {
            aligned1 = s1[i-1] + aligned1;
            aligned2 = '-' + aligned2;
            i--;
        } else {
            aligned1 = '-' + aligned1;
            aligned2 = s2[j-1] + aligned2;
            j--;
        }
    }

    return {dp[m][n], {aligned1, aligned2}};
}

// Function to perform multiple sequence alignment using the Star algorithm
vector<string> msaStar(const vector<string>& sequences) {
    int n = sequences.size();
    vector<vector<int>> scores(n, vector<int>(n, 0));
    vector<vector<pair<string, string>>> alignments(n, vector<pair<string, string>>(n));

    // Calculate all pairwise alignments
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            auto [score, alignment] = pairwiseAlignment(sequences[i], sequences[j]);
            scores[i][j] = scores[j][i] = score;
            alignments[i][j] = alignments[j][i] = alignment;
        }
    }

    // Find the central sequence (with the maximum sum of scores)
    int center = 0;
    int maxScore = numeric_limits<int>::min();
    for (int i = 0; i < n; i++) {
        int sum = 0;
        for (int j = 0; j < n; j++) {
            sum += scores[i][j];
        }
        if (sum > maxScore) {
            maxScore = sum;
            center = i;
        }
    }

    // Align all sequences with respect to the central sequence
    vector<string> result(n);
    result[center] = sequences[center];
    for (int i = 0; i < n; i++) {
        if (i != center) {
            auto& [aligned_center, aligned_i] = alignments[center][i];
            result[i] = aligned_i;
        }
    }

    // Adjust lengths
    int maxLen = 0;
    for (const auto& seq : result) {
        maxLen = max(maxLen, (int)seq.length());
    }
    for (auto& seq : result) {
        seq.append(maxLen - seq.length(), '-');
    }

    return result;
}

// Function to calculate the total score of the alignment
int calculateTotalScore(const vector<string>& alignment) {
    int score = 0;
    for (int i = 0; i < alignment[0].length(); i++) {
        for (int j = 0; j < alignment.size(); j++) {
            for (int k = j + 1; k < alignment.size(); k++) {
                score += ::score(alignment[j][i], alignment[k][i]);
            }
        }
    }
    return score;
}

// Function to read sequences from a file
vector<string> readSequencesFromFile(const string& filename, bool forwardOnly) {
    vector<string> sequences;
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        if ((forwardOnly && line.find("F:") != string::npos) || (!forwardOnly && line.find("R:") != string::npos)) {
            size_t start = line.find("5′-") + 3;
            size_t end = line.find("-3′");
            if (start != string::npos && end != string::npos) {
                sequences.push_back(line.substr(start, end - start));
            }
        }
    }
    return sequences;
}

int main() {
    // Read forward and reverse sequences from the BRCA1 file
    vector<string> forwardSequences = readSequencesFromFile("C:/Users/harry/Documents/Clion_projects/test/BRCA1.txt", true);
    vector<string> reverseSequences = readSequencesFromFile("C:/Users/harry/Documents/Clion_projects/test/BRCA1.txt", false);

    // Forward sequences alignment
    cout << "Forward sequences alignment of BRCA1:" << endl;
    auto forwardAlignment = msaStar(forwardSequences);
    for (const auto& seq : forwardAlignment) {
        cout << seq << endl;
    }
    int forwardScore = calculateTotalScore(forwardAlignment);
    cout << "Total score of forward sequences: " << forwardScore << endl << endl;

    // Reverse sequences alignment
    cout << "Reverse sequences alignment of BRCA1:" << endl;
    auto reverseAlignment = msaStar(reverseSequences);
    for (const auto& seq : reverseAlignment) {
        cout << seq << endl;
    }
    int reverseScore = calculateTotalScore(reverseAlignment);
    cout << "Total score of reverse sequences: " << reverseScore << endl;

    return 0;
}
