#pragma once
#include <string>
#include <vector>
using namespace std;

// Índices de inicio (base 0) de cada aparición de patron en texto
#pragma once
#include <string>
#include <vector>
using namespace std;

// Índices de inicio (base 0) de cada aparición de patron en texto
// kmp_buscar: algoritmo Knuth-Morris-Pratt (KMP).
// Encuentra todas las apariciones de 'patron' en 'texto' (índices base 0), incluyendo apariciones traslapadas.
inline vector<int> kmp_buscar(const string& texto, const string& patron) {
    vector<int> res;
    int n = texto.size(), m = patron.size();
    if (m == 0 || m > n) return res;

    vector<int> lps(m, 0);
    for (int i = 1, k = 0; i < m; i++) {
        while (k > 0 && patron[i] != patron[k]) k = lps[k - 1];
        if (patron[i] == patron[k]) k++;
        lps[i] = k;
    }

    for (int i = 0, k = 0; i < n; i++) {
        while (k > 0 && texto[i] != patron[k]) k = lps[k - 1];
        if (texto[i] == patron[k]) k++;
        if (k == m) {
            res.push_back(i - m + 1);
            k = lps[k - 1];
        }
    }
    return res;
}