#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include "kmp.h"
#include "codones.h"
#include "hits.h"
using namespace std;
// contar_no_acgt: cuenta caracteres que no son A, C, G ni T. 
//Complejidad: O(n)
inline int contar_no_acgt(const string& s) {
    int c = 0;
    for (char x : s)
        if (x != 'A' && x != 'C' && x != 'G' && x != 'T') c++;
    return c;
}

inline int punto4a(const string& wuhan, const string& texas, const vector<Hit>& hits) {
    cout << "=== Punto 4a: Wuhan 2019 vs Texas 2020 ===\n";
    cout << "Largo Wuhan: " << wuhan.size() << "   Largo Texas: " << texas.size() << "\n";
    cout << "Caracteres no ACGT -> Wuhan: " << contar_no_acgt(wuhan)
         << "   Texas: " << contar_no_acgt(texas) << "\n";

    // Alinear: buscar los primeros 20 caracteres de Wuhan dentro de Texas
    auto t0 = chrono::steady_clock::now();
    vector<int> m = kmp_buscar(texas, wuhan.substr(0, 20));
    auto t1 = chrono::steady_clock::now();
    if (m.empty()) {
        cout << "No se pudo alinear los genomas\n";
        return -1;
    }
    int desfase = m[0];
    cout << "Desfase: " << desfase << "\n";
    cout << "Tiempo de alineacion KMP: "
         << chrono::duration<double, milli>(t1 - t0).count() << " ms\n\n";

    auto tabla = construir_tabla_codones();
    int sinonimas = 0, cambian = 0, no_cod = 0, ambiguas = 0;

    // Comparar solo la parte que existe en ambos
    int n = min((int)wuhan.size(), (int)texas.size() - desfase);

    for (int i = 0; i < n; i++) {
        char bw = wuhan[i], bt = texas[i + desfase];
        if (bw == bt) continue;

                cout << "Pos " << i + 1 << ": " << bw << " -> " << bt << "\n";

        // Codón afectado en cada uno de los 3 frames de lectura del genoma
        for (int m = 0; m < 3; m++) {
            if (i < m) continue;
            int c = m + ((i - m) / 3) * 3;
            if (c + 3 > n) continue;
            string fw = wuhan.substr(c, 3);
            string ft = texas.substr(c + desfase, 3);
            cout << "    frame " << m << " | " << fw << " (" << tabla[fw] << ") -> "
                 << ft << " (" << tabla[ft] << ")\n";
        }
        cout << "   ";

        int idx = hit_en(hits, i);
        if (idx == -1) {
            cout << "   no codificante\n";
            no_cod++;
            continue;
        }
        const Hit& h = hits[idx];

        int ic = inicio_codon(h, i);
        if (ic < 0 || ic + 3 > n) {
            cout << "   ambigua (codon fuera de rango)\n";
            ambiguas++;
            continue;
        }

        string cw = wuhan.substr(ic, 3);
        string ct = texas.substr(ic + desfase, 3);
        cout << "   [" << h.nombre << "] codon " << cw << " -> " << ct;

        if (contar_no_acgt(cw) > 0 || contar_no_acgt(ct) > 0) {
            cout << "   ambigua\n";
            ambiguas++;
            continue;
        }

        auto aw = tabla[cw], at = tabla[ct];
        cout << "   aminoacido " << aw << " -> " << at;
        if (aw == at) { cout << "   sinonima\n";          sinonimas++; }
        else          { cout << "   cambia aminoacido\n"; cambian++;   }
    }

    // Cola de Wuhan sin contraparte en Texas (diferencia de largo, no mutaciones)
    if ((int)wuhan.size() > n) {
        cout << "\nCola de Wuhan sin contraparte en Texas: posiciones "
             << n + 1 << "-" << wuhan.size()
             << " (" << wuhan.size() - n << " bases: " << wuhan.substr(n) << ")\n";
    }

    int total = sinonimas + cambian + no_cod + ambiguas;
    cout << "\n--- Resumen ---\n";
    cout << "Total diferencias:  " << total << "\n";
    cout << "Sinonimas:          " << sinonimas << "\n";
    cout << "Cambian aminoacido: " << cambian << "\n";
    cout << "No codificantes:    " << no_cod << "\n";
    cout << "Ambiguas:           " << ambiguas << "\n";
    return desfase;
}