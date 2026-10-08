#include <iostream>
#include <string>
#include "fasta.h"
#include "kmp.h"
#include "manacher.h"
#include "codones.h"
using namespace std;

// ---------- Gerardo: punto 1 (genes) ----------
void punto1() {
    // TODO: leer genoma Wuhan y genes M, S, ORF1AB.
    // Imprimir nombre, índices (base 1) y primeros 12 caracteres.
}

// ---------- Alejandro: punto 2 (palíndromos) ----------
void punto2() {
    // TODO: manacher() en cada gen, imprimir longitud, guardar palindromos.txt
}

// ---------- Clarisa: punto 3 (base de proteínas) ----------
void punto3_base() {
    // TODO: probar traducir() y leer_proteinas()
}

int main() {
    punto1();
    punto2();
    punto3_base();
    return 0;
}