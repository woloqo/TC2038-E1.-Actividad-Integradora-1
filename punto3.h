#pragma once
// Punto 3 (Persona C): impresion de resultados y demostracion del frameshift.
// La busqueda vive en hits.h; la tabla de codones y traducir() en codones.h.
#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <utility>
#include "kmp.h"
#include "codones.h"
#include "hits.h"
using namespace std;

// Nombre comun de cada proteina. Solo se usa para mostrarlo;
// la busqueda trabaja con el ID del archivo seq-proteins.txt.
inline string nombre_comun(const string& id) {
    static const pair<string, string> NOMBRES[] = {
        {"QHD43415_1", "nsp1"},   {"QHD43415_2", "nsp2"},   {"QHD43415_3", "nsp3"},
        {"QHD43415_4", "nsp4"},   {"QHD43415_5", "nsp5"},   {"QHD43415_6", "nsp6"},
        {"QHD43415_7", "nsp7"},   {"QHD43415_8", "nsp8"},   {"QHD43415_9", "nsp9"},
        {"QHD43415_10", "nsp10"}, {"QHD43415_11", "nsp12"}, {"QHD43415_12", "nsp13"},
        {"QHD43415_13", "nsp14"}, {"QHD43415_14", "nsp15"}, {"QHD43415_15", "nsp16"},
        {"QHD43416", "S"},        {"QHD43417", "ORF3a"},    {"QHD43418", "E"},
        {"QHD43419", "M"},        {"QHD43420", "ORF6"},     {"QHD43421", "ORF7a"},
        {"QHD43422", "ORF8"},     {"QHD43423", "N"},        {"QHI42199", "ORF10"}};
    for (const auto& p : NOMBRES)
        if (p.first == id) return p.second;
    return "";
}

// Los 'cuantos' codones que empiezan en 'inicio' (base 0), separados por espacio.
inline string codones_desde(const string& genoma, int inicio, int cuantos) {
    string r;
    for (int i = 0; i < cuantos; i++) {
        if (i > 0) r += ' ';
        r += genoma.substr(inicio + 3 * i, 3);
    }
    return r;
}

// Para cada proteina: nombre, indices (base 1), marco, rango de aminoacidos
// y los primeros 4 aminoacidos con los 4 codones reales del genoma.
// Una proteina partida (frameshift) se imprime en dos tramos.
inline void imprimir_proteinas(const string& genoma,
                               const vector<pair<string, string>>& proteinas,
                               const vector<Hit>& hits) {
    cout << "(indices en base 1; marco = posicion del primer nucleotido modulo 3, contando desde 0)\n\n";
    for (const auto& [nombre, prot] : proteinas) {
        int L = prot.size();
        string comun = nombre_comun(nombre);
        cout << nombre;
        if (!comun.empty()) cout << " (" << comun << ")";
        cout << " - " << L << " aa\n";

        int cubiertos = 0;   // aminoacidos ya cubiertos por tramos anteriores
        bool hay = false;
        int numero = 0;
        for (const Hit& h : hits) {
            if (h.nombre != nombre) continue;
            hay = true;
            int n_aa = (h.fin - h.inicio + 1) / 3;
            bool partido = (n_aa != L);
            int aa_ini = partido ? cubiertos + 1 : 1;
            int aa_fin = aa_ini + n_aa - 1;
            if (partido) {
                numero++;
                cubiertos += n_aa;
                cout << "  tramo " << numero << ": ";
            } else {
                cout << "  ";
            }
            cout << h.inicio + 1 << "-" << h.fin + 1 << "   marco " << h.marco
                 << "   aa " << aa_ini << "-" << aa_fin << "\n";

            int k = n_aa < 4 ? n_aa : 4;
            cout << "    primeros " << k << " aa: " << traducir(genoma.substr(h.inicio, 3 * k), 0)
                 << "   codones: " << codones_desde(genoma, h.inicio, k) << "\n";
        }
        if (!hay) cout << "  NO SE ENCONTRO en el genoma\n";
    }
}

// Demo para el video: la secuencia deslizante TTTAAAC (UUUAAAC en RNA) y
// el traslape de 1 nt entre los dos tramos de la proteina partida (nsp12).
inline void demo_frameshift(const string& genoma, const vector<Hit>& hits) {
    const string DESLIZANTE = "TTTAAAC";
    const int N = DESLIZANTE.size();
    vector<int> pos = kmp_buscar(genoma, DESLIZANTE);

    int i1 = -1;   // primer tramo de una proteina partida
    for (size_t i = 0; i + 1 < hits.size(); i++) {
        if (hits[i].nombre == hits[i + 1].nombre && hits[i + 1].inicio == hits[i].fin) {
            i1 = (int)i;
            break;
        }
    }

    cout << "===== Demo del frameshift =====\n";
    cout << "Secuencia deslizante " << DESLIZANTE << " (UUUAAAC en RNA): "
         << pos.size() << " apariciones en el genoma\n";
    for (int p : pos) {
        cout << "  " << p + 1 << "-" << p + N;
        if (i1 >= 0) {
            int traslape = hits[i1].fin;   // base 0
            if (traslape >= p && traslape < p + N)
                cout << "   <-- aqui se traslapan los dos tramos de " << hits[i1].nombre;
        }
        cout << "\n";
    }
    if (i1 < 0) {
        cout << "No hay proteinas partidas en los resultados.\n";
        return;
    }

    const Hit& a = hits[i1];
    const Hit& b = hits[i1 + 1];
    string antes = genoma.substr(a.fin - 8, 9);    // ultimos 3 codones del tramo 1
    string despues = genoma.substr(b.inicio, 9);   // primeros 3 codones del tramo 2
    cout << "Nucleotido compartido por los dos tramos (base 1): " << a.fin + 1 << "\n";
    cout << "  tramo 1 termina (marco " << a.marco << "): " << codones_desde(genoma, a.fin - 8, 3)
         << "  ->  " << traducir(antes, 0) << "\n";
    cout << "  tramo 2 empieza (marco " << b.marco << "): " << codones_desde(genoma, b.inicio, 3)
         << "  ->  " << traducir(despues, 0) << "\n";
}

// Corre todo el punto 3: busqueda, impresion, verificacion, tiempos y demo.
inline vector<Hit> ejecutar_punto3(const string& genoma, const vector<pair<string, string>>& proteinas) {
    if (genoma.empty() || proteinas.empty()) {
        cout << "Punto 3: faltan datos (genoma o proteinas vacios)\n";
        return {};
    }
    cout << "===== PUNTO 3: proteinas en el genoma de Wuhan =====\n";
    cout << "Genoma: " << genoma.size() << " nt   Proteinas: " << proteinas.size() << "\n";

    auto t0 = chrono::steady_clock::now();
    vector<Hit> hits = buscar_proteinas(genoma, proteinas);
    auto t1 = chrono::steady_clock::now();
    double ms = chrono::duration<double, milli>(t1 - t0).count();

    imprimir_proteinas(genoma, proteinas, hits);

    int ubicadas = 0;
    for (const auto& p : proteinas) {
        for (const Hit& h : hits) {
            if (h.nombre == p.first) { ubicadas++; break; }
        }
    }
    int ok = verificar_hits(genoma, proteinas, hits);

    cout << "\nResumen: " << ubicadas << " de " << proteinas.size() << " proteinas ubicadas, "
         << hits.size() << " hits\n";
    cout << "Verificacion (traducir el rango de cada hit): " << ok << "/" << hits.size() << "\n";
    cout << "Tiempo de buscar_proteinas: " << ms << " ms\n\n";

    demo_frameshift(genoma, hits);
    return hits;
}