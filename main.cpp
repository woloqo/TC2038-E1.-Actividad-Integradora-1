#include <iostream>
#include <fstream>
#include <string>
#include "fasta.h"
#include "kmp.h"
#include "manacher.h"
#include "lcs.h"
#include "codones.h"
#include "hits.h"
using namespace std;

// ---------- PERSONA A: punto 1 (genes) ----------
void punto1() {
    // TODO: leer genoma Wuhan y genes M, S, ORF1AB.
    // Imprimir nombre, índices (base 1) y primeros 12 caracteres.
}

// ---------- PERSONA B: punto 2 (palíndromos) ----------
void punto2() {
    const string genes[] = {"M", "S", "ORF1AB"};
    ofstream salida("palindromos.txt");

    for(auto g: genes){
        auto [encabezado, seq] = leer_fasta("archivos/gen-" + g + ".txt");
        auto [ini, largo] = manacher(seq);

        cout<<"Gen "<<g<<": palindromo mas largo de "<<largo<<" (indices "<<ini+1<<"-"<<ini+largo<<")\n";

        salida<<">"<<g<<" inicio="<<ini+1<<" fin="<<ini+largo<<" longitud="<<largo<<"\n";
        for(int i=ini; i<ini+largo; i++){
            salida<<seq[i];
        }
        salida<<"\n";
    }
}

// ---------- PERSONA B: punto 4 (subcadena común), pruebas chicas ----------
void prueba_lcs() {
    string a[] = {"ACGTACGT", "AAAA", "ABCDXYZ", "ACGT", ""};
    string b[] = {"TTACGTAA", "BBBB", "XYZABCD", "TGCA", "ACGT"};
    for(int i=0; i<5; i++){
        auto r = subcadena_comun(a[i], b[i]);
        cout<<"lcs(\""<<a[i]<<"\",\""<<b[i]<<"\") = "<<r.largo<<" ia="<<r.ia<<" ib="<<r.ib<<"\n";
    }
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
    prueba_lcs();
    punto3_base();
    return 0;
}