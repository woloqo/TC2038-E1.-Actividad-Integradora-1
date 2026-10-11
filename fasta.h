#pragma once
#include <string>
#include <utility>
using namespace std;

// Devuelve (encabezado, secuencia). Salta '>' y une líneas, quitando \r \n y espacios
inline pair<string,string> leer_fasta(const string& ruta) {
    ifstream archivo(ruta);
    if (!archivo.is_open()) {
        cerr << "No pude abrir: " << ruta << "\n";
        return {"", ""};
    }

    string linea, encabezado, seq;
    while (getline(archivo, linea)) {
        if (linea.empty()) continue;
        if (linea[0] == '>') {
            encabezado = linea.substr(1);
            if (!encabezado.empty() && encabezado.back() == '\r') encabezado.pop_back();
            continue;
        }
        for (char c : linea) {
            if (!isspace((unsigned char)c)) seq += (char)toupper((unsigned char)c);
        }
    }
    return {encabezado, seq};
}