# OptimalPolyCuba4D

Questa directory contiene le implementazioni del metodo **OptimalPolyCuba4D (OPC4D)**, basato sulla costruzione di regole di cubatura mediante una base tensoriale di polinomi di Chebyshev e sul calcolo dei momenti geometrici del dominio poliedrico mobile tramite il teorema della divergenza applicato alle facce triangolari del bordo.

Rispetto al caso 3D, il dominio di integrazione non è un poliedro fisso ma un poliedro che si deforma linearmente, per un tempo adimensionale `tau` in `[0,1]`, dalla configurazione iniziale a quella finale, a parità di connettività della mesh superficiale. L'integrale è quindi calcolato sul politopo spazio-temporale `(x, y, z, tau)` generato da questo movimento.

Il metodo è disponibile in due varianti algoritmiche:

- **MoL** (Method of Lines): discretizza il tempo con nodi di Clenshaw-Curtis e, per ciascun istante `tau_k`, applica la procedura OPC3D al poliedro congelato in quell'istante;
- **Tensor**: costruisce direttamente una base tensoriale di Chebyshev in `(x, y, z, tau)` e calcola i momenti sulle ipersuperfici laterali del politopo (prismi triangolari spazio-temporali generati dal movimento di ciascuna faccia), senza discretizzare il tempo a parte.

Ciascuna variante è disponibile in tre implementazioni:

- **MATLAB** (`Matlab4D/`): orientata allo sviluppo, alla validazione numerica e alla visualizzazione;
- **Fortran** (`Fortran4D/`): implementazione seriale orientata alle prestazioni;
- **Fortran parallela** (`Parallel4D/`): versione della precedente con parallelizzazione OpenMP.

Nella variante MoL la routine principale si chiama `OPC4D_MoL`; nella variante Tensor si chiama `OPC4D_Tensor` (`OPC4D_Parallel_MoL` e `OPC4D_Parallel_Tensor` nelle versioni parallele).

---

## Struttura della directory

```text
Mesh4D/
│
├── README.md
│
├── Matlab4D/
│   ├── MoL/
│   │   ├── examples/
│   │   │   ├── concave_tri.dat
│   │   │   ├── concave_vertex.dat
│   │   │   ├── concave_vertex_new.dat
│   │   │   ├── convex_tri.dat
│   │   │   ├── convex_vertex.dat
│   │   │   ├── convex_vertex_new.dat
│   │   │   ├── example_concave.m
│   │   │   ├── example_convex.m
│   │   │   ├── example_polynomial.m
│   │   │   └── example_rigidmotion.m
│   │   └── src/
│   │       ├── Dunavant/
│   │       ├── ClenshawCurtisTime.m
│   │       ├── OPC4D_MoL.m
│   │       ├── ShiftingTriangleQuadrature.m
│   │       ├── TriangleQuadrature.m
│   │       ├── chebpolys.m
│   │       ├── chebyshev_moments_polyhedron.m
│   │       ├── cub_gausscheb_tens3D.m
│   │       ├── cubature_tens_chebyshev_facet_V.m
│   │       ├── dCHEBVAND.m
│   │       ├── mono_next_grlex.m
│   │       ├── scale_rule.m
│   │       └── tenscheb_norm2sq.m
│   └── Tensor/
│       ├── examples/
│       │   ├── concave_tri.dat
│       │   ├── concave_vertex.dat
│       │   ├── concave_vertex_new.dat
│       │   ├── convex_tri.dat
│       │   ├── convex_vertex.dat
│       │   ├── convex_vertex_new.dat
│       │   ├── example_concave.m
│       │   ├── example_convex.m
│       │   └── example_polynomial.m
│       └── src/
│           ├── Dunavant/
│           ├── OPC4D_Tensor.m
│           ├── PrismQuadrature.m
│           ├── ShiftingPrismQuadrature.m
│           ├── chebpolys.m
│           ├── chebyshev_moments_polyhedron_4D.m
│           ├── clenshaw_curtis.m
│           ├── cub_gausscheb_tens4D.m
│           ├── cubature_tens_chebyshev_facet_4D.m
│           ├── dCHEBVAND.m
│           ├── gaujac.m
│           ├── lgwt.m
│           ├── mono_next_grlex.m
│           ├── scale_rule.m
│           └── tenscheb_norm2sq.m
│
├── Fortran4D/
│   ├── MoL/
│   │   ├── Makefile
│   │   ├── examples/
│   │   │   ├── concave_tri.dat
│   │   │   ├── concave_vertex.dat
│   │   │   ├── concave_vertex_new.dat
│   │   │   ├── convex_tri.dat
│   │   │   ├── convex_vertex.dat
│   │   │   ├── convex_vertex_new.dat
│   │   │   ├── example_concave.f90
│   │   │   ├── example_convex.f90
│   │   │   └── example_polynomial.f90
│   │   └── src/
│   │       ├── CubatureFunctions.f90
│   │       ├── OPC4D_MoL.f90
│   │       ├── PolyhedronMesh.f90
│   │       ├── ReferenceFunctions.f90
│   │       ├── TimeDiscretization.f90
│   │       ├── TriangleQuadratureDunavant.f90
│   │       ├── TriangleQuadratureGJ.f90
│   │       └── TypesDef.f90
│   └── Tensor/
│       ├── Makefile
│       ├── examples/
│       │   ├── concave_tri.dat
│       │   ├── concave_vertex.dat
│       │   ├── concave_vertex_new.dat
│       │   ├── convex_tri.dat
│       │   ├── convex_vertex.dat
│       │   ├── convex_vertex_new.dat
│       │   ├── example_concave.f90
│       │   ├── example_convex.f90
│       │   └── example_polynomial.f90
│       └── src/
│           ├── CubatureTensor.f90
│           ├── OPC4D_Tensor.f90
│           ├── PolyhedronMesh.f90
│           ├── PrismQuadrature.f90
│           ├── ReferenceFunctions.f90
│           ├── TriangleQuadratureDunavant.f90
│           ├── TriangleQuadratureGJ.f90
│           └── TypesDef.f90
│
└── Parallel4D/
    ├── MoL/
    │   ├── Makefile
    │   ├── examples/
    │   │   ├── concave_tri.dat
    │   │   ├── concave_vertex.dat
    │   │   ├── concave_vertex_new.dat
    │   │   ├── convex_tri.dat
    │   │   ├── convex_vertex.dat
    │   │   ├── convex_vertex_new.dat
    │   │   ├── example_concave.f90
    │   │   ├── example_convex.f90
    │   │   └── example_polynomial.f90
    │   └── src/
    │       ├── CubatureFunctions.f90
    │       ├── OPC4D_Parallel_MoL.f90
    │       ├── PolyhedronMesh.f90
    │       ├── ReferenceFunctions.f90
    │       ├── TimeDiscretization.f90
    │       ├── TriangleQuadratureDunavant.f90
    │       ├── TriangleQuadratureGJ.f90
    │       └── TypesDef.f90
    └── Tensor/
        ├── Makefile
        ├── examples/
        │   ├── concave_tri.dat
        │   ├── concave_vertex.dat
        │   ├── concave_vertex_new.dat
        │   ├── convex_tri.dat
        │   ├── convex_vertex.dat
        │   ├── convex_vertex_new.dat
        │   ├── example_concave.f90
        │   ├── example_convex.f90
        │   └── example_polynomial.f90
        └── src/
            ├── CubatureTensor.f90
            ├── OPC4D_Parallel_Tensor.f90
            ├── PolyhedronMesh.f90
            ├── PrismQuadrature.f90
            ├── ReferenceFunctions.f90
            ├── TriangleQuadratureDunavant.f90
            ├── TriangleQuadratureGJ.f90
            └── TypesDef.f90
```

Le due varianti (MoL, Tensor) e le tre implementazioni sono organizzate separatamente, ma mantengono la stessa impostazione algoritmica di fondo e gli stessi esempi di riferimento. Ogni directory è autonoma: contiene il proprio codice sorgente, gli esempi e i dati geometrici (configurazione iniziale e finale della mesh).

### Esempi disponibili

| Esempio | Dominio | Contenuto |
|---------|---------|-----------|
| `example_convex` | Poliedro convesso in movimento | Caso di base per la validazione. |
| `example_concave` | Poliedro concavo in movimento | Verifica su un dominio non convesso. |
| `example_polynomial` | Cubo in movimento | Verifica dell'esattezza su polinomi di grado `<= ade`. |
| `example_rigidmotion` | Cubo in rototraslazione rigida | Solo `Matlab4D/MoL`. Rotazione attorno all'asse z e traslazione, senza variazione di forma né di volume. |

---

## Il metodo in breve

### Metodo delle linee (MoL)

Dato il grado polinomiale `ade` e un numero di istanti temporali `n_tau`, la procedura:

1. genera `n_tau` nodi e pesi di Clenshaw-Curtis sull'intervallo `[0,1]`;
2. per ciascun nodo temporale `tau_k`, interpola linearmente la mesh tra la configurazione iniziale e quella finale, ottenendo il poliedro congelato all'istante `tau_k`;
3. su ciascun poliedro congelato applica la stessa procedura di OPC3D: griglia tensoriale di Gauss-Chebyshev sul cubo `[-1,1]^3`, base di Chebyshev di grado totale `<= ade`, momenti calcolati con il teorema della divergenza sulle facce, riscalamento sulla bounding box istantanea;
4. combina i pesi spaziali così ottenuti con i pesi temporali di Clenshaw-Curtis.

La regola restituita ha `n_tau * (ade+1)^3` nodi.

### Metodo tensoriale (Tensor)

Dato il grado polinomiale `ade`, la procedura:

1. costruisce sull'ipercubo di riferimento `[-1,1]^4` una griglia tensoriale di Gauss-Chebyshev con `(ade+1)^4` nodi `(x,y,z,tau)`;
2. costruisce la base tensoriale di Chebyshev 4D di grado totale `<= ade` (ordinamento GRLEX), di dimensione `N_mom = (ade+1)(ade+2)(ade+3)(ade+4)/24`;
3. costruisce, per ogni faccia triangolare della mesh, il prisma spazio-temporale generato dal suo movimento tra la configurazione iniziale e quella finale, e calcola i momenti della base sul politopo 4D sommando, con il teorema della divergenza, i contributi di tutti i prismi;
4. riscala i nodi sulla bounding box 4D del politopo e ricombina i momenti per ottenere i pesi.

La regola restituita ha `(ade+1)^4` nodi ed è esatta per polinomi di grado totale fino a `ade` sul politopo spazio-temporale. I pesi possono essere negativi.

### Quadratura sulle facce

**MoL** — a ogni istante congelato gli integrali sulle facce triangolari sono calcolati come nel caso 3D, con il parametro `method`:

| `method` | Regola | Note |
|----------|--------|------|
| `'GJ'`   | Prodotto di Gauss-Jacobi sul triangolo di riferimento | Nessun limite su `ade`. |
| `'D'`    | Regole di Dunavant | Disponibili fino al grado 20, quindi `ade <= 19`. |

**Tensor** — gli integrali sono calcolati sul prisma triangolare spazio-temporale di ogni faccia, combinando una regola sul triangolo e una regola nel tempo, scelte tramite il parametro `method`:

| `method` | Triangolo | Tempo | Note |
|----------|-----------|-------|------|
| `'DCC'`  | Dunavant | Clenshaw-Curtis | `ade <= 19` per il vincolo di Dunavant. |
| `'DGL'`  | Dunavant | Gauss-Legendre | `ade <= 19` per il vincolo di Dunavant. |
| `'GJCC'` | Gauss-Jacobi | Clenshaw-Curtis | Nessun limite su `ade`. |
| `'GJL'`  | Gauss-Jacobi | Gauss-Legendre | Nessun limite su `ade`. |

Come nel caso 3D, il grado di precisione richiesto sul prisma è `ade+1` (e non `ade`), perché la primitiva usata nel teorema della divergenza aumenta di uno il grado dell'integranda.

---

## Implementazione MATLAB

### MoL

La funzione principale è:

```text
Matlab4D/MoL/src/OPC4D_MoL.m
```

```matlab
[XYZtau, W] = OPC4D_MoL(ade, n_tau, vertici_iniziali, vertici_finali, facets, method)
```

dove:

- `ade` è il grado polinomiale totale massimo;
- `n_tau` è il numero di nodi temporali di Clenshaw-Curtis;
- `vertici_iniziali`, `vertici_finali` sono matrici `N x 3` con le coordinate dei vertici della mesh, rispettivamente a `tau = 0` e `tau = 1`;
- `facets` è una matrice `M x 3` con la connettività delle facce triangolari;
- `method` seleziona la quadratura sulle facce (`'GJ'` o `'D'`).

La funzione restituisce `XYZtau` (`N_nodi x 4`, con colonne `[x, y, z, tau]`) e i relativi pesi `W`.

Funzioni ausiliarie principali in `Matlab4D/MoL/src/`: `ClenshawCurtisTime.m` (nodi e pesi temporali), `cub_gausscheb_tens3D.m`, `mono_next_grlex.m`, `chebpolys.m`, `dCHEBVAND.m`, `tenscheb_norm2sq.m`, `chebyshev_moments_polyhedron.m`, `cubature_tens_chebyshev_facet_V.m`, `TriangleQuadrature.m`, `ShiftingTriangleQuadrature.m`, `scale_rule.m`, oltre alle regole di Dunavant in `Dunavant/`.

### Tensor

La funzione principale è:

```text
Matlab4D/Tensor/src/OPC4D_Tensor.m
```

```matlab
[XYZT, W] = OPC4D_Tensor(ade, vertici_iniziali, vertici_finali, facets, method)
```

dove gli argomenti sono gli stessi del caso MoL (senza `n_tau`) e `method` è uno tra `'DCC'`, `'DGL'`, `'GJCC'`, `'GJL'`. La funzione restituisce `XYZT` (`(ade+1)^4 x 4`, colonne `[x, y, z, tau]`) e i pesi `W`.

Funzioni ausiliarie principali in `Matlab4D/Tensor/src/`: `cub_gausscheb_tens4D.m`, `mono_next_grlex.m`, `chebpolys.m`, `dCHEBVAND.m`, `tenscheb_norm2sq.m`, `chebyshev_moments_polyhedron_4D.m`, `cubature_tens_chebyshev_facet_4D.m`, `PrismQuadrature.m`, `ShiftingPrismQuadrature.m`, `clenshaw_curtis.m`, `gaujac.m`, `lgwt.m`, `scale_rule.m`, oltre alle regole di Dunavant in `Dunavant/`.

### Esecuzione degli esempi

Dalla directory `examples/` di ciascuna variante (gli script aggiungono `../` e `../src/` al path):

```matlab
example_convex
example_concave
example_polynomial
example_rigidmotion   % solo Matlab4D/MoL
```

In entrambe le varianti la funzione integranda può essere valutata dopo la costruzione della regola, ad esempio `I = W' * f(XYZT(:,1), XYZT(:,2), XYZT(:,3), XYZT(:,4))`.

---

## Implementazione Fortran

### MoL

La routine principale è implementata nel modulo `OPC4D_MoL_Module`, nel file:

```text
Fortran4D/MoL/src/OPC4D_MoL.f90
```

e si utilizza con:

```fortran
USE OPC4D_MoL_Module, ONLY: OPC4D_MoL

CALL OPC4D_MoL(ade, n_tau, vertices_initial, vertices_final, facets, method, XYZT, W)
```

| Argomento           | Tipo | Descrizione |
|---------------------|------|-------------|
| `ade`               | `INTEGER, INTENT(IN)` | Grado polinomiale totale massimo. |
| `n_tau`             | `INTEGER, INTENT(IN)` | Numero di nodi temporali di Clenshaw-Curtis. |
| `vertices_initial`  | `REAL(dp), INTENT(IN)`, `(:,:)` | Vertici a `tau = 0`, `N x 3`. |
| `vertices_final`    | `REAL(dp), INTENT(IN)`, `(:,:)` | Vertici a `tau = 1`, `N x 3`. |
| `facets`            | `INTEGER, INTENT(IN)`, `(:,:)` | Connettività delle facce, `M x 3`, indici da 1. |
| `method`            | `CHARACTER(LEN=*), INTENT(IN)` | `'GJ'` oppure `'D'`. |
| `XYZT`              | `REAL(dp), ALLOCATABLE, INTENT(OUT)`, `(:,:)` | Nodi di cubatura, `n_tau*(ade+1)^3 x 4`. |
| `W`                 | `REAL(dp), ALLOCATABLE, INTENT(OUT)`, `(:)` | Pesi di cubatura. |

#### Moduli

| File | Contenuto |
|------|-----------|
| `TypesDef.f90` | Precisione `dp` e costanti. |
| `PolyhedronMesh.f90` | Tipo `t_polyhedron` e lettura delle mesh (`MeshReader`), con configurazione iniziale e finale. |
| `ReferenceFunctions.f90` | Griglia di Gauss-Chebyshev 3D, ordinamento GRLEX, matrice di Vandermonde-Chebyshev, norme della base. |
| `TimeDiscretization.f90` | Nodi e pesi di Clenshaw-Curtis nel tempo (`ComputeClenshawCurtis`). |
| `TriangleQuadratureGJ.f90` | Quadratura di Gauss-Jacobi sul triangolo. |
| `TriangleQuadratureDunavant.f90` | Regole di Dunavant sul triangolo. |
| `CubatureFunctions.f90` | Momenti di Chebyshev sul poliedro istantaneo e sulle singole facce, riscalamento della regola. |
| `OPC4D_MoL.f90` | Routine principale `OPC4D_MoL`. |

### Tensor

La routine principale è implementata nel modulo `OPC4D_Tensor_Module`, nel file:

```text
Fortran4D/Tensor/src/OPC4D_Tensor.f90
```

e si utilizza con:

```fortran
USE OPC4D_Tensor_Module, ONLY: OPC4D_Tensor

CALL OPC4D_Tensor(ade, vertices_initial, vertices_final, facets, method, XYZT, W)
```

| Argomento           | Tipo | Descrizione |
|---------------------|------|-------------|
| `ade`               | `INTEGER, INTENT(IN)` | Grado polinomiale totale massimo. |
| `vertices_initial`  | `REAL(dp), INTENT(IN)`, `(:,:)` | Vertici a `tau = 0`, `N x 3`. |
| `vertices_final`    | `REAL(dp), INTENT(IN)`, `(:,:)` | Vertici a `tau = 1`, `N x 3`. |
| `facets`            | `INTEGER, INTENT(IN)`, `(:,:)` | Connettività delle facce, `M x 3`, indici da 1. |
| `method`            | `CHARACTER(LEN=*), INTENT(IN)` | `'DCC'`, `'DGL'`, `'GJCC'` o `'GJL'`. |
| `XYZT`              | `REAL(dp), ALLOCATABLE, INTENT(OUT)`, `(:,:)` | Nodi di cubatura, `(ade+1)^4 x 4`. |
| `W`                 | `REAL(dp), ALLOCATABLE, INTENT(OUT)`, `(:)` | Pesi di cubatura. |

#### Moduli

| File | Contenuto |
|------|-----------|
| `TypesDef.f90` | Precisione `dp` e costanti. |
| `PolyhedronMesh.f90` | Tipo `t_polyhedron` e lettura delle mesh (`MeshReader`), con configurazione iniziale e finale. |
| `ReferenceFunctions.f90` | Griglia di Gauss-Chebyshev 4D, ordinamento GRLEX, matrice di Vandermonde-Chebyshev, norme della base. |
| `TriangleQuadratureGJ.f90` | Quadratura di Gauss-Jacobi sul triangolo. |
| `TriangleQuadratureDunavant.f90` | Regole di Dunavant sul triangolo. |
| `PrismQuadrature.f90` | Regole di quadratura sul prisma triangolare spazio-temporale (`'DCC'`, `'DGL'`, `'GJCC'`, `'GJL'`) e trasformazione sui prismi fisici. |
| `CubatureTensor.f90` | Momenti di Chebyshev sul politopo 4D e sulle singole ipersuperfici laterali, riscalamento della regola. |
| `OPC4D_Tensor.f90` | Routine principale `OPC4D_Tensor`. |

In entrambe le varianti le variabili reali usano il tipo `REAL(dp)` definito in `TypesDef.f90`, e la funzione integranda non viene passata alla routine principale: come nel caso 3D, si valuta successivamente sui nodi restituiti.

### Compilazione

La compilazione di ciascuna variante è gestita dal proprio `Makefile` (`Fortran4D/MoL/Makefile`, `Fortran4D/Tensor/Makefile`). Il compilatore è `gfortran`.

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
```

Ciascun esempio può essere compilato ed eseguito direttamente con il relativo target:

```bash
make convex
make concave
make polynomial
```

Gli esempi vengono eseguiti dalla directory `examples/`, in modo che i file `.dat` possano essere caricati direttamente. Per scegliere il metodo di quadratura, modificare la chiamata a `OPC4D_MoL` o `OPC4D_Tensor` nell'esempio.

### Pulizia e informazioni

```bash
make clean       # rimuove oggetti, moduli compilati ed eseguibili
make distclean    # pulizia completa
make info         # informazioni sulla configurazione del progetto
```

---

## Implementazione Fortran parallela (OpenMP)

La directory `Parallel4D/` contiene le versioni del codice Fortran con parallelizzazione **OpenMP**, una per ciascuna variante:

```text
Parallel4D/MoL/src/OPC4D_Parallel_MoL.f90
Parallel4D/Tensor/src/OPC4D_Parallel_Tensor.f90
```

Gli altri moduli, gli esempi e i dati geometrici hanno la stessa struttura delle rispettive versioni seriali in `Fortran4D/`. Le interfacce di `OPC4D_Parallel_MoL` e `OPC4D_Parallel_Tensor` sono identiche a quelle seriali. In entrambe le varianti il ciclo sulle facce della mesh, in cui si accumulano i momenti di Chebyshev, è parallelizzato con una direttiva `!$OMP PARALLEL DO`.

### MoL

La compilazione avviene tramite `Parallel4D/MoL/Makefile`, con gli stessi target della versione seriale (`make convex`, `make concave`, `make polynomial`, `make clean`, `make distclean`, `make info`). Il numero di thread si controlla con la variabile d'ambiente `OMP_NUM_THREADS`:

```bash
export OMP_NUM_THREADS=4
make convex
```

### Tensor

La compilazione avviene tramite `Parallel4D/Tensor/Makefile`, che aggiunge rispetto alla versione MoL:

- la variabile `OMP=1` (default, OpenMP abilitato) / `OMP=0` (build seriale di riferimento), da impostare prima di `make clean && make`;
- la variabile `THREADS` per fissare il numero di thread in esecuzione senza esportare `OMP_NUM_THREADS`;
- il target `make scaling`, che esegue l'esempio indicato in `SCALING_EXAMPLE` (default `convex`) con un numero di thread crescente e ne riporta il tempo di calcolo.

```bash
make convex THREADS=4
make scaling SCALING_EXAMPLE=polynomial
```

---

## Autore

Edoardo Longo, Università degli Studi di Verona.

## Riferimenti

Per la formulazione matematica del metodo, la costruzione dei momenti sul politopo spazio-temporale, la costruzione delle regole di cubatura MoL e Tensor e i risultati numerici si rimanda alla documentazione scientifica e alla tesi associate al progetto **OptimalPolyCuba4D**.