#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cctype>
#include<limits>

using namespace std;

#define PRIME 131 // A prime number for Rabin-Karp

void naiveApproach(const string &sequence, const string &pattern);
void rabinKarp(const string &sequence, const string &pattern);
void kmpAlgorithm(const string &sequence, const string &pattern);
void printMenu();
void printResult(const string &algorithm, const string &sequence, const string &pattern, const vector<int> &starts, const vector<int> &ends, const string &complexity);
void computeLPSArray(const string &pattern, vector<int> &lps);
void clearScreen();
bool isValidDNASequence(const string &sequence);

int main() {
    string sequence, pattern;
    int choice;


    cout << "=======================================================\n";
    cout << "               DNA Pattern Matcher                    \n";
    cout << "=======================================================\n";
    do {
        cout << "Enter DNA Sequence: ";
        cin >> sequence;
        if (!isValidDNASequence(sequence)) {
            cout << "Invalid DNA sequence. Please enter a sequence containing only A, T, G, and C.\n";
        }
    } while (!isValidDNASequence(sequence));

    do {
        cout << "Enter pattern to search: ";
        cin >> pattern;
        if (!isValidDNASequence(pattern)) {
            cout << "Invalid pattern. Please enter a sequence containing only A, T, G, and C.\n";
        }
    } while (!isValidDNASequence(pattern));

    while (true) {
        clearScreen();
        printMenu();
        cout << "\nChoose an option: ";

        if (!(cin >> choice)) {
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            cout << "Invalid input. Please enter a number between 1 and 4.\n";
            continue;
        }

        switch (choice) {
            case 1:
                naiveApproach(sequence, pattern);
                break;
            case 2:
                rabinKarp(sequence, pattern);
                break;
            case 3:
                kmpAlgorithm(sequence, pattern);
                break;
            case 4:
                cout << "Exiting program. Thank you!\n";
                return 0;
            default:
                cout << "Invalid choice. Please try again.\n";
        }
    }

    return 0;
}

void printMenu() {
    cout << "\n-------------------------------------------------------\n";
    cout << "              Pattern Matching Algorithms             \n";
    cout << "-------------------------------------------------------\n";
    cout << "1. Naive Approach\n";
    cout << "2. Rabin-Karp Algorithm\n";
    cout << "3. Knuth-Morris-Pratt (KMP) Algorithm\n";
    cout << "4. Exit\n";
}

void printResult(const string &algorithm, const string &sequence, const string &pattern, const vector<int> &starts, const vector<int> &ends, const string &complexity) {
    clearScreen();

    cout << "=======================================================\n";
    cout << "               Algorithm: " << algorithm << "\n";
    cout << "=======================================================\n";
    cout << "DNA Sequence:              " << sequence << "\n";
    cout << "Pattern:                   " << pattern << "\n";

    cout << "\nNumber of occurrences: \t   " << starts.size() << "\n";

    for (size_t i = 0; i < starts.size(); i++) {
        cout << "\nInstance " << i + 1 << ":\n";
        cout << "Starting Index:            " << starts[i] << "\n";
        cout << "Ending Index:              " << ends[i] << "\n";
    }

    cout << "\nTime Complexity:           " << complexity << "\n";
    cout << "=======================================================\n";


    cout << "\nPress Enter to continue...";
    cin.ignore();
    cin.get();
}

bool isValidDNASequence(const string &sequence) {
    for (char c : sequence) {
        if (toupper(c) != 'A' && toupper(c) != 'T' && toupper(c) != 'G' && toupper(c) != 'C') {
            return false;
        }
    }
    return true;
}

void naiveApproach(const string &sequence, const string &pattern) {
    int n = sequence.length();
    int m = pattern.length();

    vector<int> starts, ends;

    for (int i = 0; i <= n - m; i++) {
        int j;
        for (j = 0; j < m; j++) {
            if (sequence[i + j] != pattern[j]) {
                break;
            }
        }
        if (j == m) {
            starts.push_back(i);
            ends.push_back(i + m - 1);
        }
    }

    if (!starts.empty()) {
        printResult("Naive Approach", sequence, pattern, starts, ends, "O(n * m)");
    } else {
        cout << "No matching subsequences found.\n";
    }
}

void rabinKarp(const string &sequence, const string &pattern) {
    int n = sequence.length();
    int m = pattern.length();

    vector<int> starts, ends;
    long long hashSeq = 0, hashPat = 0, h = 1;

    for (int i = 0; i < m - 1; i++)
        h = (h * 256) % PRIME;

    for (int i = 0; i < m; i++) {
        hashPat = (256 * hashPat + pattern[i]) % PRIME;
    }

    for (int i = 0; i <= n - m; i++) {
        if (i == 0) {
            for (int j = 0; j < m; j++)
                hashSeq = (256 * hashSeq + sequence[j]) % PRIME;
        } else {
            hashSeq = (256 * (hashSeq - sequence[i - 1] * h) + sequence[i + m - 1]) % PRIME;
            if (hashSeq < 0)
                hashSeq += PRIME;
        }

        if (hashSeq == hashPat) {
            int j;
            for (j = 0; j < m; j++) {
                if (sequence[i + j] != pattern[j]) {
                    break;
                }
            }
            if (j == m) {
                starts.push_back(i);
                ends.push_back(i + m - 1);
            }
        }
    }

    if (!starts.empty()) {
        printResult("Rabin-Karp Algorithm", sequence, pattern, starts, ends, "O(n + m)");
    } else {
        cout << "No matching subsequences found.\n";
    }
}

void kmpAlgorithm(const string &sequence, const string &pattern) {
    int n = sequence.length();
    int m = pattern.length();
    vector<int> starts, ends;
    vector<int> lps(m, 0);
    computeLPSArray(pattern, lps);

    int i = 0, j = 0;
    while (i < n) {
        if (sequence[i] == pattern[j]) {
            i++;
            j++;
        }
        if (j == m) {
            starts.push_back(i - j);
            ends.push_back(i - 1);
            j = lps[j - 1];
        } else if (i < n && sequence[i] != pattern[j]) {
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i++;
            }
        }
    }
    if (!starts.empty()) {
        printResult("Knuth-Morris-Pratt (KMP) Algorithm", sequence, pattern, starts, ends, "O(n + m)");
    } else {
        cout << "No matching subsequences found.\n";
    }
}

void computeLPSArray(const string &pattern, vector<int> &lps) {
    int length = 0;
    lps[0] = 0;

    int i = 1;
    while (i < pattern.length()) {
        if (pattern[i] == pattern[length]) {
            length++;
            lps[i] = length;
            i++;
        } else {
            if (length != 0) {
                length = lps[length - 1];
            } else {
                lps[i] = 0;
                i++;
            }
        }
    }
}

void clearScreen() {
system("cls");
}
