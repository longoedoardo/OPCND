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

## Utilizzo

### Versione Fortran 3D

```bash
cd mesh3D/Fortran3D
make
./OptimalPolyCuba3D
```

### Versione MATLAB 3D

Apri il file:

```matlab
mesh3D/Matlab3D/mainOptimalPolyCuba3D.m
```

e eseguilo da MATLAB. Il file aggiunge automaticamente le cartelle ShapeIndependent e ShapeDependent al path.

### Versione Fortran 4D (Tensor)

```bash
cd mesh4D/Fortran4D/Tensor
make
./mainOptimalPolyCuba4D_Tensor
```

### Versione Fortran 4D (MoL)

```bash
cd mesh4D/Fortran4D/MoL
make
./mainOptimalPolyCuba4D_MoL
```

### Versione MATLAB 4D

Apri uno dei file principali:

```matlab
mesh4D/Matlab4D/Tensor/mainOptimalPolyCuba4D_Tensor.m
mesh4D/Matlab4D/MoL/mainOptimalPolyCuba4D_MoL.m
```

## Struttura del progetto

```text
PolyCheapCubatureND/
├── LICENSE
├── README.md
├── mesh3D/
│   ├── Fortran3D/
│   │   ├── Makefile
│   │   ├── data/
│   │   └── src/
│   └── Matlab3D/
│       ├── ShapeDependent/
│       └── ShapeIndependent/
├── mesh4D/
│   ├── Fortran4D/
│   │   ├── MoL/
│   │   └── Tensor/
│   └── Matlab4D/
│       ├── MoL/
│       └── Tensor/
```

## Esempi

Il repository contiene già esempi di esecuzione nei file main:

- [mesh3D/Matlab3D/mainOptimalPolyCuba3D.m](mesh3D/Matlab3D/mainOptimalPolyCuba3D.m)
- [mesh3D/Fortran3D/src/mainOptimalPolyCuba3D.f90](mesh3D/Fortran3D/src/mainOptimalPolyCuba3D.f90)
- [mesh4D/Matlab4D/Tensor/mainOptimalPolyCuba4D_Tensor.m](mesh4D/Matlab4D/Tensor/mainOptimalPolyCuba4D_Tensor.m)
- [mesh4D/Fortran4D/Tensor/src/mainOptimalPolyCuba4D_Tensor.f90](mesh4D/Fortran4D/Tensor/src/mainOptimalPolyCuba4D_Tensor.f90)

Questi file usano dati di esempio già presenti nelle rispettive cartelle data/ o nella directory corrente.

## Riferimenti scientifici

Il riferimento principale alla base metodologica del progetto è l'articolo:

- Sommariva, A., Vianello, M. "Cheap and stable quadrature on polyhedral elements". DOI: https://dx.doi.org/10.1016/j.finel.2025.104409

Per una trattazione teorica completa, inclusi i dettagli matematici, le dimostrazioni e l'analisi delle varie componenti del metodo, si rimanda alla tesi associata al progetto, che raccoglie i riferimenti e le motivazioni teoriche per la formulazione del problema di cubatura su poliedri e l'estensione del metodo ai domini spazio-temporali 4D.

## Test

Al momento non sono presenti test automatizzati o workflow CI nel repository. La verifica del codice avviene attraverso l'esecuzione manuale dei programmi principali e dei target Makefile forniti.

## Note importanti

- Il progetto è principalmente orientato a ricerca e sperimentazione numerica.
- Le implementazioni Fortran richiedono gfortran e i file di input devono essere accessibili relativamente alla directory di esecuzione.
- I main contengono parametri come ade e la funzione integranda, che possono essere modificati per adattare il calcolo al caso di interesse.
- Il progetto è distribuito con licenza MIT.

## Licenza

Questo progetto è distribuito con licenza MIT. Vedere il file [LICENSE](LICENSE) per i dettagli.
