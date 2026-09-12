# OptimalPolyCuba3D

Questa directory contiene l'implementazione tridimensionale del metodo
**OptimalPolyCuba3D**, basato sulla costruzione di regole di cubatura
mediante una base tensoriale di polinomi di Chebyshev e il calcolo dei
momenti geometrici del dominio poliedrale.

---

## Struttura della directory

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
│   ├── src/
│   │   ├── chebpolys.m
│   │   ├── chebyshev_moments_polyhedron.m
│   │   ├── cub_gausscheb_tens3D.m
│   │   ├── cubature_tens_chebyshev_facet_V.m
│   │   ├── dCHEBVAND.m
│   │   ├── mapTrianglePoints.m
│   │   ├── mono_next_grlex.m
│   │   ├── scale_rule.m
│   │   ├── tenscheb_norm2sq.m
│   │   └── TriangleQuadraturePoints.m
│   │
│   ├── tri.dat
│   └── vertex.dat
│
└── Fortran3D/
    └── ...
```

La directory `Matlab3D` contiene l'implementazione MATLAB del metodo,
mentre `Fortran3D` contiene l'implementazione Fortran.

La struttura è organizzata in modo da mantenere separate le diverse
implementazioni dello stesso metodo.

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

mentre le funzioni ausiliarie utilizzate internamente dall'algoritmo
sono contenute nella directory:

```text
Matlab3D/src/
```

Gli esempi e i relativi dati geometrici sono contenuti in:

```text
Matlab3D/examples/
```

---

## Funzione principale

L'interfaccia principale è:

```matlab
I = OptimalPolyCuba3D(ade, vertices, facets, f)
```

dove:

- `ade` è il grado polinomiale totale massimo considerato nella
  procedura di matching dei momenti;
- `vertices` è una matrice `N x 3` contenente le coordinate cartesiane
  dei vertici della mesh;
- `facets` è una matrice `M x 3` contenente la connettività delle facce
  triangolari;
- `f` è un function handle MATLAB che rappresenta la funzione integranda.

La funzione restituisce l'approssimazione numerica dell'integrale. La funzione integranda deve essere vettorializzata.

---

# Rappresentazione del dominio

Il dominio tridimensionale è rappresentato mediante una **mesh
superficiale triangolare chiusa**. La geometria è descritta dalle due matrici:

```matlab
vertices
```

e

```matlab
facets
```

La matrice `vertices` contiene, per ogni riga, le coordinate di un
vertice. La matrice `facets` contiene, per ogni riga, gli indici dei tre vertici
che costituiscono una faccia triangolare.

---

## Orientamento delle facce

L'orientamento delle facce è fondamentale. I vertici di ciascun triangolo devono 
essere ordinati in modo coerente così che il prodotto vettoriale associato 
alla faccia produca una **normale orientata verso l'esterno del dominio**.
La mesh deve inoltre essere:

1. chiusa;
2. composta da facce triangolari;
3. orientata coerentemente;
4. priva di buchi.

Queste condizioni sono necessarie perché i momenti geometrici vengono
calcolati a partire dalla rappresentazione superficiale del dominio,
utilizzando una formulazione basata sul teorema della divergenza.

---

# Esempi

La directory `examples/` contiene esempi progettati per mostrare
differenti caratteristiche del metodo.

## Dominio convesso

```text
example_convex.m
```

Applica la regola di cubatura a un dominio poliedrale convesso.

La mesh utilizzata è:

```text
convex_vertex.dat
convex_tri.dat
```

Questo esempio rappresenta un caso generale di dominio convesso
triangolare.

---

## Dominio concavo

```text
example_concave.m
```

Applica la regola a un dominio poliedrale concavo.

La mesh utilizzata è:

```text
concave_vertex.dat
concave_tri.dat
```

L'esempio è particolarmente significativo perché mostra l'applicazione
del metodo a domini non convessi.

---

## Integrazione di funzioni polinomiali

```text
example_polynomial.m
```

Questo esempio considera funzioni polinomiali per le quali è possibile
calcolare analiticamente l'integrale esatto.

Il confronto tra:

- integrale esatto;
- integrale numerico;
- errore di quadratura;

permette di verificare la proprietà di esattezza polinomiale della
regola e di studiare l'effetto del parametro `ade`.

---

## Stanford Bunny

```text
example_bunny.m
```

Questo esempio utilizza la mesh triangolare del modello **Stanford
Bunny**, contenuta nel file:

```text
bun_zipper.ply
```

Il modello rappresenta una geometria complessa e realistica ed è
utilizzato principalmente per verificare il comportamento del codice
su mesh superficiali di dimensioni significativamente maggiori rispetto
agli esempi elementari.

---

# Esecuzione degli esempi

Gli script presenti nella directory `examples/` sono progettati per
essere indipendenti dalla directory di lavoro corrente di MATLAB.

Gli script determinano automaticamente la posizione della directory
principale `Matlab3D` tramite:

```matlab
projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(projectRoot);
addpath(fullfile(projectRoot, 'src'));
```

In questo modo vengono aggiunte al MATLAB path:

```text
Matlab3D/
Matlab3D/src/
```

senza assumere una specifica directory di lavoro.

---

# Visualizzazione delle mesh

Gli esempi includono la visualizzazione dei domini poliedrali utilizzati.

Per mesh di grandi dimensioni, la visualizzazione contemporanea di tutti
i bordi dei triangoli o l'utilizzo della trasparenza può aumentare
significativamente il costo del rendering.

In questi casi è possibile creare una mesh ridotta esclusivamente per
la visualizzazione, ad esempio:

```matlab
[facets_plot, vertices_plot] = reducepatch( ...
    facets, vertices, 0.01);
```

La mesh ridotta deve essere utilizzata solamente per il `plot` e non
deve sostituire `vertices` e `facets` utilizzati dalla procedura di
cubatura.

---

# Implementazione Fortran

La directory `mesh3D` è organizzata per permettere la presenza di più
implementazioni dello stesso metodo:

```text
mesh3D/
├── Matlab3D/
└── Fortran3D/
```

L'implementazione MATLAB è principalmente orientata a:

- sviluppo e verifica del metodo;
- esperimenti numerici;
- validazione;
- visualizzazione delle geometrie;
- analisi dei risultati.

L'implementazione Fortran è invece pensata per l'esecuzione ad alte
prestazioni e sarà documentata nella directory:

```text
Fortran3D/
```

---

# Riferimenti

Per la formulazione matematica del metodo, la costruzione dei momenti,
la procedura di matching e i risultati numerici si rimanda alla
documentazione scientifica e alla tesi associate al progetto
OptimalPolyCuba3D.

Se il software viene utilizzato in lavori scientifici, si invita a
citare il relativo lavoro scientifico o la tesi di riferimento.
