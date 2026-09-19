# OptimalPolyCuba3D

Questa directory contiene le implementazioni del metodo **OptimalPolyCuba3D**, basato sulla costruzione di regole di cubatura mediante una base tensoriale di polinomi di Chebyshev e sul calcolo dei momenti geometrici del dominio poliedrale tramite il teorema della divergenza applicato alle facce triangolari del bordo.

Il progetto comprende tre implementazioni dello stesso metodo:

- **MATLAB** (`Matlab3D/`): orientata allo sviluppo, alla validazione numerica e alla visualizzazione;
- **Fortran** (`Fortran3D/`): implementazione seriale orientata alle prestazioni;
- **Fortran parallela** (`Parallel3D/`): versione della precedente con parallelizzazione OpenMP.

In tutte le versioni la routine principale si chiama `OPC3D`.

---

## Struttura della directory

```text
Mesh3D/
│
├── Matlab3D/
│   ├── OPC3D.m
│   ├── examples/
│   │   ├── bun_zipper.ply
│   │   ├── bunny_tri.dat
│   │   ├── bunny_vertex.dat
│   │   ├── concave_tri.dat
│   │   ├── concave_vertex.dat
│   │   ├── convex_tri.dat
│   │   ├── convex_vertex.dat
│   │   ├── example_bunny.m
│   │   ├── example_concave.m
│   │   ├── example_convex.m
│   │   └── example_polynomial.m
│   └── src/
│       ├── chebpolys.m
│       ├── chebyshev_moments_polyhedron.m
│       ├── cub_gausscheb_tens3D.m
│       ├── cubature_tens_chebyshev_facet_V.m
│       ├── dCHEBVAND.m
│       ├── mapTriangleDunavantPoints.m
│       ├── mapTriangleGJPoints.m
│       ├── mono_next_grlex.m
│       ├── scale_rule.m
│       ├── tenscheb_norm2sq.m
│       ├── TriangleDunavantQuadraturePoints.m
│       ├── TriangleGJQuadraturePoints.m
│       └── Dunavant/
│           └── (routine di supporto per le regole di Dunavant)
│
├── Fortran3D/
│   ├── Makefile
│   ├── examples/
│   │   ├── bunny_tri.dat
│   │   ├── bunny_vertex.dat
│   │   ├── concave_tri.dat
│   │   ├── concave_vertex.dat
│   │   ├── convex_tri.dat
│   │   ├── convex_vertex.dat
│   │   ├── example_bunny.f90
│   │   ├── example_concave.f90
│   │   ├── example_convex.f90
│   │   └── example_polynomial.f90
│   └── src/
│       ├── CubatureFunctions.f90
│       ├── OPC3D.f90
│       ├── PolyhedronMesh.f90
│       ├── ReferenceFunctions.f90
│       ├── TriangleQuadratureDunavant.f90
│       ├── TriangleQuadratureGJ.f90
│       └── TypesDef.f90
│
└── Parallel3D/
    ├── Makefile
    ├── examples/
    │   └── (stessi esempi e dati di Fortran3D)
    └── src/
        ├── CubatureFunctions.f90
        ├── OPC3D_Parallel.f90
        ├── PolyhedronMesh.f90
        ├── ReferenceFunctions.f90
        ├── TriangleQuadratureDunavant.f90
        ├── TriangleQuadratureGJ.f90
        └── TypesDef.f90
```

Le tre implementazioni sono organizzate separatamente, ma mantengono la stessa impostazione algoritmica e gli stessi esempi di riferimento.

---

## Il metodo in breve

Dato il grado polinomiale `ade`, la procedura:

1. costruisce sul cubo di riferimento `[-1,1]^3` una griglia tensoriale di Gauss-Chebyshev con `(ade+1)^3` nodi;
2. costruisce la base tensoriale di Chebyshev di grado totale `<= ade` (ordinamento GRLEX), di dimensione `N_mom = (ade+1)(ade+2)(ade+3)/6`;
3. calcola i momenti della base sul poliedro trasformandoli, con il teorema della divergenza, in integrali sulle facce triangolari della mesh;
4. riscala i nodi sulla bounding box del poliedro e ricombina i momenti per ottenere i pesi.

La regola restituita ha `(ade+1)^3` nodi ed è esatta per polinomi di grado totale fino a `ade` sul poliedro. I pesi possono essere negativi.

### Quadratura sulle facce

Gli integrali sulle facce sono calcolati con una regola di quadratura sul triangolo, scelta tramite il parametro `method`:

| `method` | Regola | Note |
|----------|--------|------|
| `'GJ'`   | Prodotto conico di Gauss-Jacobi sul triangolo di riferimento | `nGP = ceil((ade+2)/2)` punti per direzione, cioè `nGP^2` punti per faccia. Nessun limite su `ade`. |
| `'D'`    | Regole di Dunavant | Si usa la regola di grado `ade+1`; le regole sono disponibili fino al grado 20, quindi `ade <= 19`. |

Il grado di precisione richiesto sulle facce è `ade+1` (e non `ade`) perché la primitiva usata nel teorema della divergenza aumenta di uno il grado dell'integranda.

---

## Implementazione MATLAB

La funzione principale è:

```text
Matlab3D/OPC3D.m
```

Le funzioni ausiliarie sono in `Matlab3D/src/` (le routine di supporto per le regole di Dunavant sono in `Matlab3D/src/Dunavant/`), mentre gli esempi e i dati geometrici sono in `Matlab3D/examples/`.

L'interfaccia è:

```matlab
[XYZ, W] = OPC3D(ade, vertices, facets, method)
```

dove:

- `ade` è il grado polinomiale totale massimo considerato nella procedura di matching dei momenti;
- `vertices` è una matrice `N x 3` con le coordinate cartesiane dei vertici della mesh;
- `facets` è una matrice `M x 3` con la connettività delle facce triangolari;
- `method` seleziona la quadratura sulle facce (`'GJ'` o `'D'`).

La funzione restituisce:

- `XYZ`, matrice `(ade+1)^3 x 3` con le coordinate dei nodi di cubatura;
- `W`, vettore con i relativi pesi.

La funzione integranda viene valutata dopo la costruzione della regola:

```matlab
I = W' * f(XYZ(:,1), XYZ(:,2), XYZ(:,3));
```

La separazione tra costruzione della regola e valutazione dell'integranda permette di riutilizzare la stessa regola per funzioni diverse.

---

## Implementazione Fortran

La routine principale è implementata nel modulo `OPC3D_Module`, nel file:

```text
Fortran3D/src/OPC3D.f90
```

e si utilizza con:

```fortran
USE OPC3D_Module, ONLY: OPC3D

CALL OPC3D(ade, vertices, facets, method, XYZ, W)
```

dove:

| Argomento  | Tipo | Descrizione |
|------------|------|-------------|
| `ade`      | `INTEGER, INTENT(IN)` | Grado polinomiale totale massimo. |
| `vertices` | `REAL(dp), INTENT(IN)`, `(:,:)` | Coordinate dei vertici, `N x 3`. |
| `facets`   | `INTEGER, INTENT(IN)`, `(:,:)` | Connettività delle facce, `M x 3`, con indici che partono da 1. |
| `method`   | `CHARACTER(LEN=*), INTENT(IN)` | `'GJ'` oppure `'D'`. |
| `XYZ`      | `REAL(dp), ALLOCATABLE, INTENT(OUT)`, `(:,:)` | Nodi di cubatura, `(ade+1)^3 x 3`. |
| `W`        | `REAL(dp), ALLOCATABLE, INTENT(OUT)`, `(:)` | Pesi di cubatura. |

Le variabili reali usano il tipo `REAL(dp)` (doppia precisione) definito in `TypesDef.f90`.

La funzione integranda non viene passata a `OPC3D`. Come nella versione MATLAB, viene valutata successivamente sui punti di cubatura:

```fortran
Integrale = 0.0_dp
DO i = 1, SIZE(W)
    Integrale = Integrale + W(i) * f(XYZ(i,1), XYZ(i,2), XYZ(i,3))
END DO
```

### Moduli

| File | Contenuto |
|------|-----------|
| `TypesDef.f90` | Precisione `dp` e costanti. |
| `PolyhedronMesh.f90` | Tipo `t_polyhedron` e lettura delle mesh (`MeshReader`). |
| `ReferenceFunctions.f90` | Griglia di Gauss-Chebyshev, ordinamento GRLEX, matrice di Vandermonde-Chebyshev, norme della base. |
| `TriangleQuadratureGJ.f90` | Quadratura di Gauss-Jacobi sul triangolo. |
| `TriangleQuadratureDunavant.f90` | Regole di Dunavant sul triangolo. |
| `CubatureFunctions.f90` | Momenti di Chebyshev sul poliedro e sulle singole facce, riscalamento della regola. |
| `OPC3D.f90` | Routine principale `OPC3D`. |

### Compilazione

La compilazione è gestita dal `Makefile` presente in `Fortran3D/`. Il compilatore è `gfortran`.

Gli oggetti compilati vengono collocati in `build/`, i moduli `.mod` in `modules_build/`; entrambe le directory sono create automaticamente.

Per compilare tutti gli esempi:

```bash
make
```

oppure:

```bash
make examples
```

Vengono generati gli eseguibili:

```text
example_convex
example_concave
example_polynomial
example_bunny
```

Ciascun esempio può essere compilato ed eseguito direttamente con il relativo target:

```bash
make convex
make concave
make polynomial
make bunny
```

Gli esempi vengono eseguiti dalla directory `examples/`, in modo che i file `.dat` possano essere caricati direttamente.

Per scegliere il metodo di quadratura, modificare la chiamata a `OPC3D` nell'esempio (`'GJ'` oppure `'D'`).

### Verifica rapida

L'esempio `convex` integra la funzione costante `f = 1` con `ade = 1`. Il risultato deve coincidere con il volume della mesh, circa `1.7743437424807`.

### Pulizia e informazioni

Per rimuovere oggetti, moduli compilati ed eseguibili:

```bash
make clean
```

Per una pulizia completa:

```bash
make distclean
```

Per visualizzare le principali informazioni sulla configurazione del progetto:

```bash
make info
```

---

## Implementazione Fortran parallela

La directory `Parallel3D/` contiene la versione del codice Fortran con parallelizzazione OpenMP. La routine principale si trova in:

```text
Parallel3D/src/OPC3D_Parallel.f90
```

Gli altri moduli, gli esempi e i dati geometrici hanno la stessa struttura di `Fortran3D/`. La compilazione avviene tramite il `Makefile` presente in `Parallel3D/`, con gli stessi target della versione seriale (`make convex`, `make bunny`, ...).

Il numero di thread si controlla con la variabile d'ambiente `OMP_NUM_THREADS`:

```bash
export OMP_NUM_THREADS=4
make bunny
```

---

## Rappresentazione del dominio

Il dominio tridimensionale è rappresentato mediante una **mesh superficiale triangolare chiusa**.

La geometria è descritta dalle due matrici:

```text
vertices
facets
```

La matrice `vertices` contiene, per ogni riga, le coordinate di un vertice. La matrice `facets` contiene, per ogni riga, gli indici dei tre vertici che costituiscono una faccia triangolare.

---

## Orientamento delle facce

I vertici di ciascun triangolo devono essere ordinati in modo coerente, affinché il prodotto vettoriale

```text
(V2 - V1) x (V3 - V1)
```

produca una **normale orientata verso l'esterno del dominio**.

La mesh deve essere:

1. chiusa;
2. composta da facce triangolari;
3. orientata coerentemente;
4. rappresentativa del bordo del dominio.

Queste condizioni sono necessarie perché i momenti geometrici vengono calcolati dalla rappresentazione superficiale del dominio con il teorema della divergenza. Il codice non verifica né corregge l'orientamento: una mesh non chiusa o orientata in modo incoerente produce risultati errati.

---

## Formato dei dati geometrici

Per le implementazioni Fortran le mesh vengono fornite mediante due file:

```text
*_vertex.dat
*_tri.dat
```

Il file dei vertici contiene prima il numero di vertici, poi le coordinate:

```text
N
x_1 y_1 z_1
x_2 y_2 z_2
...
x_N y_N z_N
```

Il file delle facce contiene prima il numero di facce, poi la connettività triangolare (indici dei vertici a partire da 1):

```text
M
i_1 j_1 k_1
i_2 j_2 k_2
...
i_M j_M k_M
```

La routine `MeshReader` carica i due file e costruisce il tipo `t_polyhedron`, che contiene i vertici, le facce e la bounding box.

---

## Visualizzazione delle mesh

La visualizzazione tridimensionale è implementata negli esempi MATLAB, che permettono di mostrare le mesh mediante `patch`.

Per mesh di grandi dimensioni il rendering può essere semplificato con una versione ridotta della mesh:

```matlab
[facets_plot, vertices_plot] = reducepatch(facets, vertices, 0.01);
```

La mesh ridotta deve essere utilizzata **solamente per la visualizzazione** e non deve sostituire la mesh completa utilizzata dalla procedura di cubatura.

Le implementazioni Fortran sono focalizzate sulla costruzione della regola di cubatura e non includono una componente grafica.

---

## Autore

Edoardo Longo, Università degli Studi di Verona.

## Riferimenti

Per la formulazione matematica del metodo, la costruzione dei momenti, la costruzione delle regole di cubatura e i risultati numerici si rimanda alla documentazione scientifica e alla tesi associate al progetto **OptimalPolyCuba3D**.