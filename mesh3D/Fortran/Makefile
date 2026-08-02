# =====================================================================
# Makefile per LALO3D
# Progetto: Cubatura 3D su poliedri concavi e convessi
# Autore: Edoardo Longo
# Data: 02/08/2026
# =====================================================================

# =====================================================================
# SEZIONE 1 — NOME DELL'ESEGUIBILE
# =====================================================================
EXE = LALO3D

# =====================================================================
# SEZIONE 2 — COMPILATORE
# =====================================================================
F90 = gfortran

# =====================================================================
# SEZIONE 3 — DIRECTORY DI BUILD
# =====================================================================
OBJ_DIR = build
MOD_DIR = modules_build
SRC_DIR = src
OUT_DIR = output

# =====================================================================
# SEZIONE 4 — FLAG DI COMPILAZIONE
# =====================================================================
BASE_FLAGS = -O3 \
             -fdefault-real-8 \
             -fdefault-double-8 \
             -ffree-line-length-none \
             -Wall \
             -Wextra \
             -J$(MOD_DIR) \
             -I$(MOD_DIR)

SERIAL_FLAGS = $(BASE_FLAGS)
OMP_FLAGS    = $(BASE_FLAGS) -fopenmp
F90FLAGS     = $(SERIAL_FLAGS)

# =====================================================================
# SEZIONE 5 — NUMERO DI THREAD OpenMP
# =====================================================================
OMP_NUM_THREADS = 4

# =====================================================================
# SEZIONE 6 — LISTA DEI MODULI
# =====================================================================
MOD_OBJ = \
    $(OBJ_DIR)/TypesDef.o \
    $(OBJ_DIR)/PrepCheap.o \
    $(OBJ_DIR)/PolyhedronMesh.o \
    $(OBJ_DIR)/CubaCheap.o \
    $(OBJ_DIR)/triangleQuadratureGJ.o

# Oggetto del programma principale
MAIN_OBJ = $(OBJ_DIR)/mainLALO3D.o

# Lista completa di tutti gli oggetti da linkare
OBJS = $(MOD_OBJ) $(MAIN_OBJ)

# =====================================================================
# SEZIONE 7 — TARGET PHONY
# =====================================================================
.PHONY: all serial openmp setup clean distclean run run_omp info

# =====================================================================
# SEZIONE 8 — TARGET PRINCIPALI
# =====================================================================
all: setup $(EXE)

serial: F90FLAGS = $(SERIAL_FLAGS)
serial: setup $(EXE)
	@echo ""
	@echo ">>> Compilazione SERIALE completata <<<"

openmp: F90FLAGS = $(OMP_FLAGS)
openmp: setup $(EXE)
	@echo ""
	@echo ">>> Compilazione con OpenMP completata <<<"
	@echo "    Thread disponibili: $(OMP_NUM_THREADS)"

# =====================================================================
# SEZIONE 9 — CREAZIONE DELLE DIRECTORY DI BUILD
# =====================================================================
setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)
	@mkdir -p $(OUT_DIR)

# =====================================================================
# SEZIONE 10 — LINKER
# =====================================================================
$(EXE): $(OBJS)
	@echo ""
	@echo "======================================"
	@echo " Linking finale → $(EXE)"
	@echo "======================================"
	$(F90) $(F90FLAGS) $(OBJS) -o $(EXE)
	@echo " Eseguibile pronto: ./$(EXE)"
	@echo "======================================"

# =====================================================================
# SEZIONE 11 — REGOLA GENERICA DI COMPILAZIONE
# =====================================================================
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@echo "  Compiling: $<"
	$(F90) $(F90FLAGS) -c $< -o $@

# =====================================================================
# SEZIONE 12 — DIPENDENZE ESPLICITE
# =====================================================================
$(OBJ_DIR)/TypesDef.o: $(SRC_DIR)/TypesDef.f90

$(OBJ_DIR)/PrepCheap.o: $(SRC_DIR)/PrepCheap.f90 \
    $(OBJ_DIR)/TypesDef.o

$(OBJ_DIR)/PolyhedronMesh.o: $(SRC_DIR)/PolyhedronMesh.f90 \
    $(OBJ_DIR)/TypesDef.o

$(OBJ_DIR)/triangleQuadratureGJ.o: $(SRC_DIR)/triangleQuadratureGJ.f90 \
    $(OBJ_DIR)/TypesDef.o

$(OBJ_DIR)/CubaCheap.o: $(SRC_DIR)/CubaCheap.f90 \
    $(OBJ_DIR)/TypesDef.o \
    $(OBJ_DIR)/triangleQuadratureGJ.o

$(OBJ_DIR)/mainLALO3D.o: $(SRC_DIR)/mainLALO3D.f90 \
    $(OBJ_DIR)/TypesDef.o \
    $(OBJ_DIR)/PrepCheap.o \
    $(OBJ_DIR)/PolyhedronMesh.o \
	$(OBJ_DIR)/triangleQuadratureGJ.o \
	$(OBJ_DIR)/CubaCheap.o

# =====================================================================
# SEZIONE 13 — TARGET DI ESECUZIONE
# =====================================================================
run: all
	@echo ""
	@echo ">>> Avvio $(EXE) in modalità SERIALE <<<"
	./$(EXE)

run_omp: openmp
	@echo ""
	@echo ">>> Avvio $(EXE) con OpenMP — Thread: $(OMP_NUM_THREADS) <<<"
	OMP_NUM_THREADS=$(OMP_NUM_THREADS) ./$(EXE)

# =====================================================================
# SEZIONE 14 — TARGET DI PULIZIA
# =====================================================================
clean:
	@echo "Pulizia cartelle di build..."
	rm -rf $(OBJ_DIR) $(MOD_DIR)
	@echo "Rimosso: $(OBJ_DIR)/ $(MOD_DIR)/"

distclean: clean
	@echo "Pulizia completa (distclean)..."
	rm -f $(EXE)
	rm -f $(OUT_DIR)/*.dat $(OUT_DIR)/*.out $(OUT_DIR)/*.txt
	rm -f punti.txt *.dat *.out
	@echo "Progetto riportato allo stato sorgente."

# =====================================================================
# SEZIONE 15 — INFO
# =====================================================================
info:
	@echo ""
	@echo "==============================="
	@echo " LALO3D — Configurazione Build (Geometria Solo)"
	@echo "==============================="
	@echo " Eseguibile  : $(EXE)"
	@echo " Compilatore : $(F90)"
	@echo " Flag        : $(F90FLAGS)"
	@echo " Sorgenti    : $(SRC_DIR)/"
	@echo " Oggetti     : $(OBJ_DIR)/"
	@echo " Moduli .mod : $(MOD_DIR)/"
	@echo " Output      : $(OUT_DIR)/"
	@echo " OMP threads : $(OMP_NUM_THREADS)"
	@echo "==============================="
	@echo ""