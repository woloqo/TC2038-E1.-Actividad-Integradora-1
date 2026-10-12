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
   auto [encGenoma, genoma] = leer_fasta("archivos/SARS-COV-2-MN908947.3.txt");
    if (genoma.empty()) {
        cout << "No se pudo leer el genoma\n";
        return;
    }

    const string genes[] = {"M", "S", "ORF1AB"};
    for (const string& g : genes) {
        auto [enc, seq] = leer_fasta("archivos/gen-" + g + ".txt");
        vector<int> pos = kmp_buscar(genoma, seq);

        cout << "Gen " << g << "\n";
        cout << "  Indices (base 1): ";
        if (pos.empty()) cout << "no aparece";
        for (size_t i = 0; i < pos.size(); i++) {
            if (i > 0) cout << ", ";
            cout << pos[i] + 1 << "-" << pos[i] + seq.size();
        }
        cout << "\n  Primeros 12: " << seq.substr(0, 12) << "\n";
    }
}

// ---------- PERSONA B: punto 2 (palíndromos) ----------
void punto2() {
    const string genes[] = {"M", "S", "ORF1AB"};
    ofstream salida("palindromos.txt");

    for(auto g: genes){
        auto [encabezado, seq] = leer_fasta("archivos/gen-"+ g+ ".txt");
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
    auto [enc, genoma] = leer_fasta("archivos/SARS-COV-2-MN908947.3.txt");
    auto prot = leer_proteinas("archivos/seq-proteins.txt");
    auto hits = buscar_proteinas(genoma, prot);
    cout << "hits: " << hits.size() << " (esperado 25)" << endl;
    cout << "verificar: " << verificar_hits(genoma, prot, hits) << "/" << hits.size() << endl;
}

int main() {
    //punto1();
    //punto2();
    //prueba_lcs();
    punto3_base();
    return 0;
}