#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include "kmp.h"
#include "codones.h"
using namespace std;

// inicio y fin en base 0, fin inclusive. marco = 0, 1 o 2
struct Hit { string nombre; int inicio; int fin; int marco; };

// Busca 'sub' (aminoácidos) en las 3 traducciones del genoma.
// Devuelve pares (inicio en nucleótidos base 0, marco).
inline vector<pair<int,int>> buscar_en_marcos(const vector<string>& trad, const string& sub) {
    vector<pair<int,int>> res;
    for (int m = 0; m < 3; m++)
        for (int idx : kmp_buscar(trad[m], sub))
            res.push_back({m + 3 * idx, m});
    return res;
}

inline vector<Hit> buscar_proteinas(const string& genoma,const vector<pair<string,string>>& proteinas) {
    // Se traduce el genoma una sola vez por marco
    vector<string> trad(3);
    for (int m = 0; m < 3; m++) trad[m] = traducir(genoma, m);

    const int MIN_TRAMO = 5;   // largo mínimo (en aminoácidos) de cada tramo del frameshift
    vector<Hit> hits;

    for (const auto& [nombre, prot] : proteinas) {
        int L = prot.size();

        // 1) La proteína completa en un solo marco
        auto completas = buscar_en_marcos(trad, prot);
        if (!completas.empty()) {
            for (auto [ini, m] : completas)
                hits.push_back({nombre, ini, ini + 3 * L - 1, m});
            continue;
        }

        // 2) Frameshift -1: la proteína viene en dos tramos y el segundo
        //    empieza justo en el último nucleótido del primero (traslape de 1 nt)
        bool partida = false;
        for (int k = MIN_TRAMO; k <= L - MIN_TRAMO && !partida; k++) {
            auto h1 = buscar_en_marcos(trad, prot.substr(0, k));
            if (h1.empty()) continue;
            auto h2 = buscar_en_marcos(trad, prot.substr(k));
            for (auto [ini1, m1] : h1) {
                int fin1 = ini1 + 3 * k - 1;
                for (auto [ini2, m2] : h2) {
                    if (ini2 == fin1) {
                        int fin2 = ini2 + 3 * (L - k) - 1;
                        hits.push_back({nombre, ini1, fin1, m1});
                        hits.push_back({nombre, ini2, fin2, m2});
                        partida = true;
                        break;
                    }
                }
                if (partida) break;
            }
        }
        if (!partida) cerr << "No se encontro la proteina " << nombre << endl;
    }
    return hits;
}

// Chequeo: traducir el rango de cada Hit debe dar exactamente la proteína (o su tramo).
// Devuelve cuántos hits pasan. Debe ser igual a hits.size().
inline int verificar_hits(const string& genoma, const vector<pair<string,string>>& proteinas, const vector<Hit>& hits) {
    int bien = 0;
    for (size_t i = 0; i < hits.size(); i++) {
        const Hit& h = hits[i];
        string aa = traducir(genoma.substr(h.inicio, h.fin - h.inicio + 1), 0);
        string full;
        for (const auto& p : proteinas) if (p.first == h.nombre) full = p.second;
        bool tramo1 = (i + 1 < hits.size() && hits[i + 1].nombre == h.nombre);
        bool tramo2 = (i > 0 && hits[i - 1].nombre == h.nombre);
        bool ok = tramo1 ? (aa == full.substr(0, aa.size()))
                : tramo2 ? (aa == full.substr(full.size() - aa.size()))
                         : (aa == full);
        if (ok) bien++;
        else cerr << "FALLA: " << h.nombre << " " << h.inicio + 1 << "-" << h.fin + 1 << endl;
    }
    return bien;
}

// Índice del hit que contiene la posición pos (base 0), o -1 si no cae en ninguna proteína
inline int hit_en(const vector<Hit>& hits, int pos) {
    for (size_t i = 0; i < hits.size(); i++)
        if (pos >= hits[i].inicio && pos <= hits[i].fin) return (int)i;
    return -1;
}

// Inicio (base 0) del codón que contiene pos, dentro de ese hit
inline int inicio_codon(const Hit& h, int pos) {
    return pos - ((pos - h.inicio) % 3);
}