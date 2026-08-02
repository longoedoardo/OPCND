# PolyCheapCubatureND

**PolyCheapCubatureND** è una libreria numerica progettata per il calcolo efficiente di regole di cubatura su domini poliedrici multidimensionali. Sviluppata in **MATLAB** e **Fortran**, la libreria offre implementazioni a basso costo computazionale ad alte prestazioni per problemi di integrazione.

* **Approccio Tetrahedra-Free:** Entrambi i moduli (`mesh3D` e `mesh4d`) eliminano completamente la necessità di ricorrere alla tradizionale triangolazione o tetraedrizzazione del dominio interno. Sfruttando il teorema della divergenza, l'integrale volumetrico viene ridotto a un'integrazione di superficie (sulle facce poliedriche in 3D o sulle iperfacce prismatiche spazio-temporali in 4D), rendendo il metodo nativamente applicabile a **domini geometrici generici, sia convessi che concavi**.
* **Tecniche di Cubatura Avanzate:** 
  * Uso di **griglie tensoriali di Gauss-Chebyshev** ad alta efficienza per la costruzione della regola di cubatura di riferimento.
  * Utilizzo di regole di **Gauss-Jacobi** per la quadratura di precisione sulle facce triangolari.
  * Sviluppo di momenti di Chebyshev analitici e numerici.
---

## 📂 Struttura del Progetto

La repository è organizzata in modo modulare per separare le implementazioni dimensionali e i rispettivi ambienti di sviluppo:

```text
PolyCheapCubatureND/
│
├── mesh3D/
│   ├── matlab/      # Codice sorgente, funzioni e script MATLAB per l'integrazione 3D
│   └── fortran/     # Codice sorgente, moduli e Makefile per l'integrazione 3D in Fortran
│
└── mesh4d/
    ├── matlab/      # Codice sorgente e script MATLAB per domini spazio-temporali 4D
    └── fortran/     # Codice sorgente, moduli e Makefile (es. mainLALO4D) per 4D
