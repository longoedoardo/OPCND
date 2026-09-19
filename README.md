# OptimalPolyCubatureND

OptimalPolyCubatureND è un progetto di ricerca orientato al calcolo numerico di regole di cubatura su domini poliedrici in 3D e 4D. Il codice presenta implementazioni in MATLAB e Fortran per costruire quadrature su poliedri convessi e concavi tramite una formulazione basata sul teorema della divergenza.

## Funzionalità

- Implementazioni per cubatura su poliedri 3D e domini spazio-temporali 4D
- Supporto a due approcci 4D: Tensor e MoL
- Uso di griglie tensoriali di Gauss-Chebyshev e quadrature di Gauss-Jacobi
- Esempi di esecuzione già inclusi nei file principali del progetto
- Build automatizzata tramite Makefile per le versioni Fortran

## Tecnologie utilizzate

- MATLAB
- Fortran 90
- GNU Make

## Prerequisiti

Per utilizzare il progetto sono richiesti:

- MATLAB (per le versioni MATLAB)
- Un compilatore Fortran, in particolare gfortran
- Un sistema Unix-like con make (il repository include Makefile attivi su macOS/Linux)

## Installazione

Clona il repository:

```bash
git clone https://github.com/your-username/PolyCheapCubatureND.git
cd PolyCheapCubatureND
```

## Configurazione

Il progetto contiene già esempi di input nella cartella dei dati e nei file principali. In particolare:

- Le versioni Fortran leggono file come vertex.dat, tri.dat e vertex_new.dat relativi alla cartella di esecuzione
- Le versioni MATLAB leggono file analoghi nella cartella corrente
- Il grado di esattezza algebrica è impostato direttamente nei file main, ad esempio nel parametro ade

Se si desidera cambiare la geometria o l'integranda, è necessario modificare i file di input o il codice nel main corrispondente.

## Struttura del progetto

```text
PolyCheapCubatureND/
├── LICENSE
├── README.md
├── mesh3D/
│   ├── Fortran3D/
│   ├── Matlab3D/
│   └── Parallel3D/
├── mesh4D/
│   ├── Fortran4D/
│   │   ├── MoL/
│   │   └── Tensor/
│   └── Matlab4D/
│       ├── MoL/
│       └── Tensor/
```

## Riferimenti scientifici

Il riferimento principale alla base metodologica del progetto è l'articolo:

- Sommariva, A., Vianello, M. "Cheap and stable quadrature on polyhedral elements". DOI: https://dx.doi.org/10.1016/j.finel.2025.104409

Per una trattazione teorica completa, inclusi i dettagli matematici, le dimostrazioni e l'analisi delle varie componenti del metodo, si rimanda alla tesi associata al progetto, che raccoglie i riferimenti e le motivazioni teoriche per la formulazione del problema di cubatura su poliedri e l'estensione del metodo ai domini spazio-temporali 4D.

## Licenza

Questo progetto è distribuito con licenza MIT. Vedere il file [LICENSE](LICENSE) per i dettagli.
