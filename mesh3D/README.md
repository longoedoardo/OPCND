# OptimalPolyCuba3D

Questa directory contiene le implementazioni del metodo **OptimalPolyCuba3D**, basato sulla costruzione di regole di cubatura mediante una base tensoriale di polinomi di Chebyshev e sul calcolo dei momenti geometrici del dominio poliedrale.

Il progetto comprende due implementazioni dello stesso metodo:

- **MATLAB**, orientata principalmente allo sviluppo, alla validazione numerica e alla visualizzazione;
- **Fortran**, orientata a una implementazione ad alte prestazioni.

---

# Struttura della directory

```text
mesh3D/
│
├── Matlab3D/
│   ├── examples/
│   │   ├── bun_zipper.ply
│   │   ├── concave_tri.dat
│   │   ├── concave_vertex.dat
│   │   ├── convex_tri.dat
│   │   ├── convex_vertex.dat
│   │   ├── example_bunny.m
│   │   ├── example_concave.m
│   │   ├── example_convex.m
│   │   └── example_polynomial.m
│   │
│   ├── OptimalPolyCuba3D.m
│   │
│   └── src/
│       ├── chebpolys.m
│       ├── chebyshev_moments_polyhedron.m
│       ├── cub_gausscheb_tens3D.m
│       ├── cubature_tens_chebyshev_facet_V.m
│       ├── dCHEBVAND.m
│       ├── mapTrianglePoints.m
│       ├── mono_next_grlex.m
│       ├── scale_rule.m
│       ├── tenscheb_norm2sq.m
│       └── TriangleQuadraturePoints.m
│
└── Fortran3D/
    ├── Makefile
    │
    ├── examples/
    │   ├── bunny_tri.dat
    │   ├── bunny_vertex.dat
    │   ├── concave_tri.dat
    │   ├── concave_vertex.dat
    │   ├── convex_tri.dat
    │   ├── convex_vertex.dat
    │   ├── example_bunny.f90
    │   ├── example_concave.f90
    │   ├── example_convex.f90
    │   └── example_polynomial.f90
    │
    └── src/
        ├── CubaCheap.f90
        ├── OptimalPolyCuba3D.f90
        ├── PolyhedronMesh.f90
        ├── PrepCheap.f90
        ├── TypesDef.f90
        └── triangleQuadratureGJ.f90
```

La directory `Matlab3D` contiene l'implementazione MATLAB del metodo, mentre `Fortran3D` contiene l'implementazione Fortran.

Le due implementazioni sono organizzate separatamente, ma mantengono la stessa impostazione algoritmica e gli stessi esempi di riferimento.

---

# Implementazione MATLAB

L'implementazione MATLAB si trova nella directory:

```text
Matlab3D/
```

La funzione principale è:

```text
Matlab3D/OptimalPolyCuba3D.m
```

mentre le funzioni ausiliarie utilizzate internamente dall'algoritmo sono contenute nella directory:

```text
Matlab3D/src/
```

Gli esempi e i relativi dati geometrici sono contenuti in:

```text
Matlab3D/examples/
```

## Funzione principale

L'interfaccia principale è:

```matlab
[XYZ, W] = OptimalPolyCuba3D(ade, vertices, facets)
```

dove:

- `ade` è il grado polinomiale totale massimo considerato nella procedura di matching dei momenti;
- `vertices` è una matrice `N x 3` contenente le coordinate cartesiane dei vertici della mesh;
- `facets` è una matrice `M x 3` contenente la connettività delle facce triangolari.

La funzione restituisce:

- `XYZ`, matrice contenente le coordinate dei nodi di cubatura;
- `W`, vettore contenente i relativi pesi.

La funzione integranda viene valutata successivamente sui nodi restituiti dalla procedura:

```matlab
I = W' * f(XYZ(:,1), XYZ(:,2), XYZ(:,3));
```

La separazione tra costruzione della regola di cubatura e valutazione della funzione integranda permette di utilizzare la stessa regola per funzioni diverse.

---

# Implementazione Fortran

L'implementazione Fortran si trova nella directory:

```text
Fortran3D/
```

La funzione principale è implementata nel modulo:

```text
Fortran3D/src/OptimalPolyCuba3D.f90
```

e viene utilizzata tramite:

```fortran
CALL OptimalPolyCuba3D(ade, vertices, facets, XYZ, W)
```

dove:

- `ade` è il grado polinomiale totale massimo considerato;
- `vertices` contiene le coordinate dei vertici della mesh;
- `facets` contiene la connettività delle facce triangolari;
- `XYZ` contiene i nodi di cubatura restituiti;
- `W` contiene i relativi pesi.

La funzione integranda non viene passata alla routine `OptimalPolyCuba3D`. Come nella versione MATLAB, viene valutata successivamente sui punti di cubatura.

## Compilazione

La compilazione viene gestita tramite il `Makefile` presente nella directory:

```text
Fortran3D/
```

Il compilatore utilizzato è:

```text
gfortran
```

Le directory di compilazione vengono create automaticamente:

```text
build/
modules_build/
```

Gli oggetti compilati vengono collocati in:

```text
build/
```

mentre i moduli Fortran `.mod` vengono collocati in:

```text
modules_build/
```

Per compilare tutti gli esempi:

```bash
make
```

oppure:

```bash
make examples
```

Al termine della compilazione vengono generati:

```text
example_convex
example_concave
example_polynomial
example_bunny
```

È possibile compilare ed eseguire direttamente ciascun esempio tramite i relativi target del `Makefile`.

### Convex

```bash
make convex
```

### Concave

```bash
make concave
```

### Polynomial

```bash
make polynomial
```

### Bunny

```bash
make bunny
```

Gli esempi vengono eseguiti dalla directory `examples/`, in modo che i relativi file `.dat` possano essere caricati direttamente.

---

##  Pulizia

Per rimuovere gli oggetti, i moduli compilati e gli eseguibili:

```bash
make clean
```

Per eseguire una pulizia completa:

```bash
make distclean
```

Per visualizzare le principali informazioni relative alla configurazione del progetto:

```bash
make info
```

---


# Rappresentazione del dominio

Il dominio tridimensionale è rappresentato mediante una **mesh superficiale triangolare chiusa**.

La geometria è descritta dalle due matrici:

```text
vertices
facets
```

La matrice `vertices` contiene, per ogni riga, le coordinate di un vertice. La matrice `facets` contiene, per ogni riga, gli indici dei tre vertici che costituiscono una faccia triangolare.

---

# Orientamento delle facce

I vertici di ciascun triangolo devono essere ordinati in modo coerente affinché il prodotto vettoriale associato alla faccia produca una **normale orientata verso l'esterno del dominio**.

La mesh deve essere:

1. chiusa;
2. composta da facce triangolari;
3. orientata coerentemente;
4. rappresentativa del bordo del dominio.

Queste condizioni sono necessarie perché i momenti geometrici vengono calcolati a partire dalla rappresentazione superficiale del dominio mediante una formulazione basata sul teorema della divergenza.

---

# Formato dei dati geometrici

Per l'implementazione Fortran, le mesh vengono fornite mediante due file:

```text
*_vertex.dat
*_tri.dat
```

Il file dei vertici contiene inizialmente il numero di vertici, seguito dalle coordinate:

```text
N
x_1 y_1 z_1
x_2 y_2 z_2
...
x_N y_N z_N
```

Il file delle facce contiene inizialmente il numero di facce, seguito dalla connettività triangolare:

```text
M
i_1 j_1 k_1
i_2 j_2 k_2
...
i_M j_M k_M
```

Questa rappresentazione permette alla routine `MeshReader` di caricare direttamente la geometria e costruire il tipo `t_polyhedron`.

---

# Esempi

La directory

```text
Matlab3D/examples/
Fortran3D/examples/
```

contiene quattro esempi principali, progettati per verificare il funzionamento della regola di cubatura su geometrie e funzioni integrande differenti. Gli esempi sono sostanzialmente corrispondenti nelle due implementazioni, così da permettere un confronto diretto tra MATLAB e Fortran.

---

## Dominio convesso

### MATLAB

File:

```text
Matlab3D/examples/example_convex.m
```

Dati geometrici:

```text
Matlab3D/examples/convex_vertex.dat
Matlab3D/examples/convex_tri.dat
```

L'esempio applica la regola di cubatura a un dominio poliedrale convesso rappresentato mediante una mesh triangolare chiusa e orientata.

### Fortran

File:

```text
Fortran3D/examples/example_convex.f90
```

Dati geometrici:

```text
Fortran3D/examples/convex_vertex.dat
Fortran3D/examples/convex_tri.dat
```

L'esempio riproduce la stessa procedura della versione MATLAB. La mesh viene caricata mediante `MeshReader` e la regola di cubatura viene costruita attraverso:

```fortran
CALL OptimalPolyCuba3D(ade, poly%vertici, poly%facce, XYZ, W)
```

---

## Dominio concavo

### MATLAB

File:

```text
Matlab3D/examples/example_concave.m
```

Dati geometrici:

```text
Matlab3D/examples/concave_vertex.dat
Matlab3D/examples/concave_tri.dat
```

L'esempio considera un dominio poliedrale concavo e applica la stessa procedura di cubatura utilizzata per il caso convesso.

Il test è particolarmente significativo perché verifica l'applicazione del metodo a geometrie non convesse, nelle quali la rappresentazione mediante una mesh triangolare della superficie e il calcolo dei momenti geometrici assumono un ruolo fondamentale.

Anche in questo caso vengono restituiti i nodi `XYZ` e i pesi `W`, dai quali viene calcolato l'integrale della funzione considerata.

### Fortran

File:

```text
Fortran3D/examples/example_concave.f90
```

Dati geometrici:

```text
Fortran3D/examples/concave_vertex.dat
Fortran3D/examples/concave_tri.dat
```

La versione Fortran riproduce il corrispondente test MATLAB utilizzando la stessa mesh triangolare e la stessa procedura di costruzione della regola di cubatura.

La geometria viene caricata tramite `MeshReader` e passata a `OptimalPolyCuba3D`. Il test permette quindi di verificare che l'implementazione Fortran sia in grado di trattare correttamente anche domini non convessi.

Il confronto con l'esempio MATLAB consente inoltre di verificare la coerenza numerica delle due implementazioni a parità di geometria, grado di esattezza e funzione integranda.

---

## Integrazione di funzioni polinomiali


Il dominio considerato è il cubo la cui geometria viene definita direttamente all'interno dell'esempio.

Vengono considerate tre funzioni polinomiali e i relativi integrali.

### MATLAB

File:

```text
Matlab3D/examples/example_polynomial.m
```

Per ciascuna funzione viene costruita la regola di cubatura con un grado `ade` adeguato e l'integrale viene valutato sui nodi restituiti da `OptimalPolyCuba3D`.

Per ogni test vengono riportati:

- grado di esattezza `ade`;
- funzione integranda;
- integrale esatto;
- integrale numerico;
- errore assoluto;
- tempo di calcolo.

L'esempio costituisce quindi un test diretto dell'esattezza polinomiale della regola di cubatura.

### Fortran

File:

```text
Fortran3D/examples/example_polynomial.f90
```

Per ciascun test viene chiamata `OptimalPolyCuba3D` e, successivamente, la funzione integranda viene valutata nei nodi di quadratura `XYZ`. L'integrale numerico viene quindi calcolato mediante il prodotto scalare tra i valori della funzione e i pesi:

```fortran
Integrale = DOT_PRODUCT(W, fXYZ)
```

Per ogni funzione vengono riportati:

- grado `ade`;
- integrale esatto;
- integrale numerico;
- errore assoluto;
- tempo di calcolo.

La corrispondenza con l'esempio MATLAB permette di verificare direttamente l'esattezza polinomiale e la coerenza numerica delle due implementazioni.

---

## Stanford Bunny

### MATLAB

File:

```text
Matlab3D/examples/example_bunny.m
```

Dati geometrici:

```text
Matlab3D/examples/bun_zipper.ply
```

L'esempio utilizza la mesh triangolare del modello **Stanford Bunny**, una geometria tridimensionale complessa costituita da un numero elevato di vertici e facce.

La mesh viene inizialmente caricata dal file `.ply` e successivamente utilizzata per costruire la regola di cubatura mediante `OptimalPolyCuba3D`.

L'esempio rappresenta quindi un test su una geometria realistica e significativamente più complessa rispetto ai domini convesso e concavo utilizzati negli esempi precedenti.

La versione MATLAB può inoltre essere utilizzata per la visualizzazione della mesh e dei punti di quadratura, sfruttando gli strumenti grafici disponibili nell'ambiente MATLAB.

### Fortran

File:

```text
Fortran3D/examples/example_bunny.f90
```

Dati geometrici:

```text
Fortran3D/examples/bunny_vertex.dat
Fortran3D/examples/bunny_tri.dat
```

La versione Fortran utilizza la stessa geometria del Stanford Bunny, precedentemente convertita dal formato `.ply` alla rappresentazione basata su due file `.dat`:

```text
bunny_vertex.dat
bunny_tri.dat
```

Il file `bunny_vertex.dat` contiene le coordinate dei vertici della mesh, mentre `bunny_tri.dat` contiene la connettività delle facce triangolari.

La mesh viene caricata mediante `MeshReader` e successivamente passata a `OptimalPolyCuba3D`. L'esempio permette di valutare il comportamento dell'implementazione Fortran su una mesh superficiale di dimensioni maggiori rispetto ai test elementari.

A differenza della versione MATLAB, l'esempio Fortran è focalizzato sul calcolo numerico della cubatura e sulla valutazione dei tempi di esecuzione; la visualizzazione della geometria e dei punti di quadratura rimane affidata alla versione MATLAB.

---

# Visualizzazione delle mesh

La visualizzazione tridimensionale delle geometrie è attualmente implementata negli esempi MATLAB.

Gli esempi MATLAB permettono di visualizzare direttamente le mesh mediante `patch`.

Per mesh di grandi dimensioni, la visualizzazione può essere semplificata utilizzando una versione ridotta della mesh esclusivamente per il rendering:

```matlab
[facets_plot, vertices_plot] = reducepatch( ...
    facets, vertices, 0.01);
```

La mesh ridotta deve essere utilizzata **solamente per la visualizzazione** e non deve sostituire la mesh completa utilizzata dalla procedura di cubatura.

L'implementazione Fortran è invece focalizzata sulla costruzione ed esecuzione della regola di cubatura e non include attualmente una componente grafica.

---

# Riferimenti

Per la formulazione matematica del metodo, la costruzione dei momenti, la procedura di matching, la costruzione delle regole di cubatura e i risultati numerici si rimanda alla documentazione scientifica e alla tesi associate al progetto **OptimalPolyCuba3D**.

Se il software viene utilizzato in lavori scientifici, si invita a citare il relativo lavoro scientifico o la tesi di riferimento.
