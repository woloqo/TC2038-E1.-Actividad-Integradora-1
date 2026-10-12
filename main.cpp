#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include "diferencias.h"
#include "fasta.h"
#include "kmp.h"
#include "manacher.h"
#include "lcs.h"
#include "codones.h"
#include "hits.h"
#include "punto3.h"
using namespace std;

// ---------- PERSONA A: punto 1 (genes) ----------
void punto1(const string& wuhan) {
    cout << "===== PUNTO 1: genes en el genoma de Wuhan =====\n";
    const string genes[] = {"M", "S", "ORF1AB"};
    double tiempo_kmp_ms = 0.0;
    for (const string& gen : genes) {
        auto [encabezado, secuencia] = leer_fasta("archivos/gen-" + gen + ".txt");
        auto t0 = chrono::steady_clock::now();
        vector<int> posiciones = kmp_buscar(wuhan, secuencia);
        auto t1 = chrono::steady_clock::now();
        tiempo_kmp_ms += chrono::duration<double, milli>(t1 - t0).count();
        cout << "Gen " << gen << ": ";
        if (secuencia.empty() || posiciones.empty()) {
            cout << "no encontrado\n";
            continue;
        }
        int inicio = posiciones.front();
        cout << "posiciones " << inicio + 1 << "-" << inicio + (int)secuencia.size()
             << ", primeros 12 nt: " << secuencia.substr(0, 12) << "\n";
    }
    cout << "Tiempo de busqueda KMP de genes: " << tiempo_kmp_ms << " ms\n\n";
}

// ---------- PERSONA B: punto 2 (palíndromos) ----------
pair<int, int> expandir_centros(const string& secuencia) {
    int mejor_inicio = 0;
    int mejor_largo = 0;
    auto expandir = [&](int izquierda, int derecha) {
        while (izquierda >= 0 && derecha < (int)secuencia.size() &&
               secuencia[izquierda] == secuencia[derecha]) {
            int largo = derecha - izquierda + 1;
            if (largo > mejor_largo) {
                mejor_inicio = izquierda;
                mejor_largo = largo;
            }
            --izquierda;
            ++derecha;
        }
    };
    for (int centro = 0; centro < (int)secuencia.size(); ++centro) {
        expandir(centro, centro);
        expandir(centro, centro + 1);
    }
    return {mejor_inicio, mejor_largo};
}

void punto2() {
    cout << "===== PUNTO 2: palindromos =====\n";
    const string genes[] = {"M", "S", "ORF1AB"};
    ofstream salida("palindromos.txt");
    if (!salida) cerr << "No se pudo crear palindromos.txt\n";

    for (const string& gen : genes) {
        auto [encabezado, secuencia] = leer_fasta("archivos/gen-" + gen + ".txt");
        auto t0 = chrono::steady_clock::now();
        auto resultado_lineal = manacher(secuencia);
        auto t1 = chrono::steady_clock::now();
        auto resultado_cuadratico = expandir_centros(secuencia);
        auto t2 = chrono::steady_clock::now();
        int inicio = resultado_lineal.first;
        int largo = resultado_lineal.second;

        cout << "Gen " << gen << ": longitud " << largo << " (posiciones " << inicio + 1
             << "-" << inicio + largo << ")";
        cout << ", expandir-centro " << resultado_cuadratico.second
             << (resultado_cuadratico.second == largo ? " [coincide]" : " [NO coincide]") << "\n";
        cout << "  Manacher O(n): "
             << chrono::duration<double, milli>(t1 - t0).count() << " ms; expandir centros O(n^2): "
             << chrono::duration<double, milli>(t2 - t1).count() << " ms\n";

        if (resultado_cuadratico.second != largo)
            cerr << "Error: los algoritmos no coinciden para " << gen << "\n";
        if (salida) {
            salida << ">" << gen << " inicio=" << inicio + 1 << " longitud=" << largo << "\n";
            salida << secuencia.substr(inicio, largo) << "\n";
        }
    }
    cout << "\n";
}

struct Racha {
    int inicio_wuhan;
    int inicio_texas;
    int largo;
    int diferencia_antes;
    int diferencia_despues;
};

Racha racha_mas_larga(const string& wuhan, const string& texas, int desfase) {
    Racha mejor{0, 0, 0, -1, -1};
    int actual_inicio = 0;
    int actual_largo = 0;
    int diferencia_anterior = -1;
    int limite = min((int)wuhan.size(), (int)texas.size() - desfase);
    for (int i = 0; i < limite; ++i) {
        if (wuhan[i] == texas[i + desfase]) {
            if (actual_largo == 0) actual_inicio = i;
            ++actual_largo;
        } else {
            if (actual_largo > mejor.largo)
                mejor = {actual_inicio, actual_inicio + desfase, actual_largo,
                         diferencia_anterior, i};
            actual_largo = 0;
            diferencia_anterior = i;
        }
    }
    if (actual_largo > mejor.largo)
        mejor = {actual_inicio, actual_inicio + desfase, actual_largo,
                 diferencia_anterior, limite};
    return mejor;
}

void punto4(const string& wuhan, const string& texas, const vector<Hit>& hits) {
    cout << "===== PUNTO 4: comparacion de genomas =====\n";
    int desfase = punto4a(wuhan, texas, hits);
    if (desfase < 0) return;

    cout << "\n--- Punto 4b: subcadena comun mas larga ---\n";
    auto t0 = chrono::steady_clock::now();
    Comun comun = subcadena_comun(wuhan, texas);
    auto t1 = chrono::steady_clock::now();
    cout << "LCS: longitud " << comun.largo;
    if (comun.largo > 0)
        cout << ", Wuhan " << comun.ia + 1 << "-" << comun.ia + comun.largo
             << ", Texas " << comun.ib + 1 << "-" << comun.ib + comun.largo;
    cout << "\nTiempo de subcadena_comun: "
         << chrono::duration<double, milli>(t1 - t0).count() << " ms\n";

    Racha racha = racha_mas_larga(wuhan, texas, desfase);
    cout << "Racha alineada sin diferencias: longitud " << racha.largo
         << ", Wuhan " << racha.inicio_wuhan + 1 << "-" << racha.inicio_wuhan + racha.largo
         << ", Texas " << racha.inicio_texas + 1 << "-" << racha.inicio_texas + racha.largo << "\n";
    if (racha.diferencia_antes >= 0)
        cout << "Diferencia anterior: posicion " << racha.diferencia_antes + 1 << "\n";
    if (racha.diferencia_despues < (int)wuhan.size())
        cout << "Diferencia siguiente: posicion " << racha.diferencia_despues + 1 << "\n";
    bool coincide = comun.largo == racha.largo && comun.ia == racha.inicio_wuhan &&
                    comun.ib == racha.inicio_texas;
    cout << "Cruce LCS/racha: " << (coincide ? "coinciden" : "no coinciden") << "\n\n";
}

int main() {
    auto [encabezado_wuhan, wuhan] = leer_fasta("archivos/SARS-COV-2-MN908947.3.txt");
    auto [encabezado_texas, texas] = leer_fasta("archivos/SARS-COV-2-MT106054.1.txt");
    auto proteinas = leer_proteinas("archivos/seq-proteins.txt");

    if (wuhan.empty() || texas.empty() || proteinas.empty()) {
        cerr << "No se pudieron cargar todos los archivos de entrada desde archivos/\n";
        return 1;
    }

    cout << fixed << setprecision(3);
    punto1(wuhan);
    punto2();
    vector<Hit> hits = ejecutar_punto3(wuhan, proteinas);
    punto4(wuhan, texas, hits);
    return 0;
}