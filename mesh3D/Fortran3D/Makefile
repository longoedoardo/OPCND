# =====================================================================
# Makefile per OptimalPolyCuba3D
# Progetto: Cubatura 3D su poliedri concavi e convessi
# Autore: Edoardo Longo
# Data: 06/08/2026
# =====================================================================

# =====================================================================
# SEZIONE 1 — NOME DELL'ESEGUIBILE
# =====================================================================
EXE = OptimalPolyCuba3D

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

# =====================================================================
# SEZIONE 4 — FLAG DI COMPILAZIONE
# =====================================================================
F90FLAGS = -O3 \
           -fdefault-real-8 \
           -fdefault-double-8 \
           -ffree-line-length-none \
           -Wall \
           -Wextra \
           -J$(MOD_DIR) \
           -I$(MOD_DIR)

# =====================================================================
# SEZIONE 5 — LISTA DEI MODULI
# =====================================================================
MOD_OBJ = \
    $(OBJ_DIR)/TypesDef.o \
    $(OBJ_DIR)/PrepCheap.o \
    $(OBJ_DIR)/PolyhedronMesh.o \
    $(OBJ_DIR)/CubaCheap.o \
    $(OBJ_DIR)/triangleQuadratureGJ.o

# Oggetto del programma principale
MAIN_OBJ = $(OBJ_DIR)/mainOptimalPolyCuba3D.o

# Lista completa di tutti gli oggetti da linkare
OBJS = $(MOD_OBJ) $(MAIN_OBJ)

# =====================================================================
# SEZIONE 6 — TARGET PHONY
# =====================================================================
.PHONY: all setup clean distclean run info

# =====================================================================
# SEZIONE 7 — TARGET PRINCIPALI
# =====================================================================
all: setup $(EXE)
	@echo ""
	@echo ">>> Compilazione completata <<<"

# =====================================================================
# SEZIONE 8 — CREAZIONE DELLE DIRECTORY DI BUILD
# =====================================================================
setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)

# =====================================================================
# SEZIONE 9 — LINKER
# =====================================================================
$(EXE): $(OBJS)
	@echo ""
	@echo "======================================"
	@echo " Linking finale : $(EXE)"
	@echo "======================================"
	$(F90) $(F90FLAGS) $(OBJS) -o $(EXE)
	@echo " Eseguibile pronto: ./$(EXE)"
	@echo "======================================"

# =====================================================================
# SEZIONE 10 — REGOLA GENERICA DI COMPILAZIONE
# =====================================================================
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@echo "  Compiling: $<"
	$(F90) $(F90FLAGS) -c $< -o $@

# =====================================================================
# SEZIONE 11 — DIPENDENZE ESPLICITE
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

$(OBJ_DIR)/mainOptimalPolyCuba3D.o: $(SRC_DIR)/mainOptimalPolyCuba3D.f90 \
    $(OBJ_DIR)/TypesDef.o \
    $(OBJ_DIR)/PrepCheap.o \
    $(OBJ_DIR)/PolyhedronMesh.o \
    $(OBJ_DIR)/triangleQuadratureGJ.o \
    $(OBJ_DIR)/CubaCheap.o

# =====================================================================
# SEZIONE 12 — TARGET DI ESECUZIONE
# =====================================================================
run: all
	@echo ""
	@echo ">>> Avvio $(EXE) <<<"
	./$(EXE)

# =====================================================================
# SEZIONE 13 — TARGET DI PULIZIA
# =====================================================================
clean:
	@echo "Pulizia cartelle di build..."
	rm -rf $(OBJ_DIR) $(MOD_DIR)
	@echo "Rimosso: $(OBJ_DIR)/ $(MOD_DIR)/"

distclean: clean
	@echo "Pulizia completa (distclean)..."
	rm -f $(EXE)
	rm -f punti.txt *.dat *.out
	@echo "Progetto riportato allo stato sorgente."

# =====================================================================
# SEZIONE 14 — INFO
# =====================================================================
info:
	@echo ""
	@echo "==============================="
	@echo "	    OptimalPolyCuba3D"
	@echo "==============================="
	@echo " Eseguibile  : $(EXE)"
	@echo " Compilatore : $(F90)"
	@echo " Flag        : $(F90FLAGS)"
	@echo " Sorgenti    : $(SRC_DIR)/"
	@echo " Oggetti     : $(OBJ_DIR)/"
	@echo " Moduli .mod : $(MOD_DIR)/"
	@echo "==============================="
	@echo ""