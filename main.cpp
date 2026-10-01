#include <iostream>
#include <string>
#include "fasta.h"
#include "kmp.h"
#include "manacher.h"
#include "codones.h"
using namespace std;

// ---------- PERSONA A: punto 1 (genes) ----------
void punto1() {
    // TODO: leer genoma Wuhan y genes M, S, ORF1AB.
    // Imprimir nombre, índices (base 1) y primeros 12 caracteres.
}

// ---------- PERSONA B: punto 2 (palíndromos) ----------
void punto2() {
    // TODO: manacher() en cada gen, imprimir longitud, guardar palindromos.txt
}

// ---------- PERSONA C: punto 3 (base de proteínas) ----------
void punto3_base() {
    auto t = construir_tabla_codones();
    cout << t["ATG"] << t["GAG"] << t["TAA"] << " " << t.size() << endl;   // ME* 64
    cout << traducir("ATGGAGAGCCTT", 0) << endl;                           // MESL
    cout << traducir("AATGGAGAGCCTT", 1) << endl;                          // MESL, marco 1

    auto prot = leer_proteinas("archivos/seq-proteins.txt");
    cout << prot.size() << " " << prot[0].first << " " << prot[0].second.substr(0, 17) << endl;   // ~24 QHD43415_1 MESLVPGFNEKTHVQLS
}

int main() {
    punto1();
    punto2();
    punto3_base();
    return 0;
}