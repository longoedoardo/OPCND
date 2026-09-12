# =====================================================================
#
# Makefile per OptimalPolyCuba3D
#
# Progetto: Cubatura 3D su poliedri concavi e convessi
#
# Autore: Edoardo Longo
#
# Data: 12/09/2026
#
# =====================================================================


# =====================================================================
# SEZIONE 1 — NOME DEGLI ESEMPI
# =====================================================================

EXAMPLES = \
	example_convex \
	example_concave \
	example_polynomial \
	example_bunny


# =====================================================================
# SEZIONE 2 — COMPILATORE
# =====================================================================

F90 = gfortran


# =====================================================================
# SEZIONE 3 — DIRECTORY DI BUILD
# =====================================================================

OBJ_DIR     = build
MOD_DIR     = modules_build
SRC_DIR     = src
EXAMPLE_DIR = examples


# =====================================================================
# SEZIONE 4 — FLAG DI COMPILAZIONE
# =====================================================================

F90FLAGS = -O0 \
	-fdefault-real-8 \
	-fdefault-double-8 \
	-ffree-line-length-none \
	-Wall \
	-Wextra \
	-fcheck=all \
	-fbacktrace \
	-J$(MOD_DIR) \
	-I$(MOD_DIR)


# =====================================================================
# SEZIONE 5 — LISTA DEI MODULI E DELLE SUBROUTINE
# =====================================================================

MOD_OBJ = \
	$(OBJ_DIR)/TypesDef.o \
	$(OBJ_DIR)/PrepCheap.o \
	$(OBJ_DIR)/PolyhedronMesh.o \
	$(OBJ_DIR)/triangleQuadratureGJ.o \
	$(OBJ_DIR)/CubaCheap.o \
	$(OBJ_DIR)/OptimalPolyCuba3D.o


# =====================================================================
# SEZIONE 6 — OGGETTI DEGLI ESEMPI
# =====================================================================

EXAMPLE_OBJ = \
	$(OBJ_DIR)/example_convex.o \
	$(OBJ_DIR)/example_concave.o \
	$(OBJ_DIR)/example_polynomial.o \
	$(OBJ_DIR)/example_bunny.o


# =====================================================================
# SEZIONE 7 — TARGET PHONY
# =====================================================================

.PHONY: all setup \
	examples \
	convex \
	concave \
	polynomial \
	bunny \
	run \
	clean \
	distclean \
	info


# =====================================================================
# SEZIONE 8 — TARGET PRINCIPALE
# =====================================================================

all: examples


# =====================================================================
# SEZIONE 9 — CREAZIONE DELLE DIRECTORY DI BUILD
# =====================================================================

setup:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(MOD_DIR)


# =====================================================================
# SEZIONE 10 — REGOLA GENERICA DI COMPILAZIONE DEI MODULI
#
# Compila i sorgenti presenti nella directory src/ generando gli
# oggetti nella directory build/ e i moduli Fortran nella directory
# modules_build/.
#
# =====================================================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@echo ""
	@echo "  Compiling: $<"
	$(F90) $(F90FLAGS) -c $< -o $@


# =====================================================================
# SEZIONE 11 — DIPENDENZE ESPLICITE DEI MODULI
#
# L'ordine delle dipendenze è necessario perché i file Fortran
# utilizzano moduli definiti in altri sorgenti.
#
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


$(OBJ_DIR)/OptimalPolyCuba3D.o: $(SRC_DIR)/OptimalPolyCuba3D.f90 \
	$(OBJ_DIR)/TypesDef.o \
	$(OBJ_DIR)/PrepCheap.o \
	$(OBJ_DIR)/CubaCheap.o


# =====================================================================
# SEZIONE 12 — COMPILAZIONE DEGLI ESEMPI
# =====================================================================

$(OBJ_DIR)/example_convex.o: $(EXAMPLE_DIR)/example_convex.f90 \
	$(OBJ_DIR)/TypesDef.o \
	$(OBJ_DIR)/PolyhedronMesh.o \
	$(OBJ_DIR)/OptimalPolyCuba3D.o
	@echo ""
	@echo "  Compiling example: $<"
	$(F90) $(F90FLAGS) -c $< -o $@


$(OBJ_DIR)/example_concave.o: $(EXAMPLE_DIR)/example_concave.f90 \
	$(OBJ_DIR)/TypesDef.o \
	$(OBJ_DIR)/PolyhedronMesh.o \
	$(OBJ_DIR)/OptimalPolyCuba3D.o
	@echo ""
	@echo "  Compiling example: $<"
	$(F90) $(F90FLAGS) -c $< -o $@


$(OBJ_DIR)/example_polynomial.o: $(EXAMPLE_DIR)/example_polynomial.f90 \
	$(OBJ_DIR)/TypesDef.o \
	$(OBJ_DIR)/PolyhedronMesh.o \
	$(OBJ_DIR)/OptimalPolyCuba3D.o
	@echo ""
	@echo "  Compiling example: $<"
	$(F90) $(F90FLAGS) -c $< -o $@


$(OBJ_DIR)/example_bunny.o: $(EXAMPLE_DIR)/example_bunny.f90 \
	$(OBJ_DIR)/TypesDef.o \
	$(OBJ_DIR)/PolyhedronMesh.o \
	$(OBJ_DIR)/OptimalPolyCuba3D.o
	@echo ""
	@echo "  Compiling example: $<"
	$(F90) $(F90FLAGS) -c $< -o $@


# =====================================================================
# SEZIONE 13 — LINKING DEGLI ESEMPI
# =====================================================================

example_convex: $(MOD_OBJ) $(OBJ_DIR)/example_convex.o
	@echo ""
	@echo "======================================"
	@echo " Linking example : example_convex"
	@echo "======================================"
	$(F90) $(F90FLAGS) $(MOD_OBJ) $(OBJ_DIR)/example_convex.o -o $@
	@echo " Esempio pronto: ./example_convex"
	@echo "======================================"


example_concave: $(MOD_OBJ) $(OBJ_DIR)/example_concave.o
	@echo ""
	@echo "======================================"
	@echo " Linking example : example_concave"
	@echo "======================================"
	$(F90) $(F90FLAGS) $(MOD_OBJ) $(OBJ_DIR)/example_concave.o -o $@
	@echo " Esempio pronto: ./example_concave"
	@echo "======================================"


example_polynomial: $(MOD_OBJ) $(OBJ_DIR)/example_polynomial.o
	@echo ""
	@echo "======================================"
	@echo " Linking example : example_polynomial"
	@echo "======================================"
	$(F90) $(F90FLAGS) $(MOD_OBJ) $(OBJ_DIR)/example_polynomial.o -o $@
	@echo " Esempio pronto: ./example_polynomial"
	@echo "======================================"


example_bunny: $(MOD_OBJ) $(OBJ_DIR)/example_bunny.o
	@echo ""
	@echo "======================================"
	@echo " Linking example : example_bunny"
	@echo "======================================"
	$(F90) $(F90FLAGS) $(MOD_OBJ) $(OBJ_DIR)/example_bunny.o -o $@
	@echo " Esempio pronto: ./example_bunny"
	@echo "======================================"


# =====================================================================
# SEZIONE 14 — COMPILAZIONE DI TUTTI GLI ESEMPI
# =====================================================================

examples: setup $(EXAMPLES)
	@echo ""
	@echo "======================================"
	@echo " Tutti gli esempi sono stati compilati"
	@echo "======================================"
	@echo ""
	@echo "Eseguibili disponibili:"
	@echo "  ./example_convex"
	@echo "  ./example_concave"
	@echo "  ./example_polynomial"
	@echo "  ./example_bunny"
	@echo ""


# =====================================================================
# SEZIONE 15 — ESECUZIONE DEI SINGOLI ESEMPI
# =====================================================================

convex: setup example_convex
	@echo ""
	@echo ">>> Avvio esempio CONVEX <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_convex


concave: setup example_concave
	@echo ""
	@echo ">>> Avvio esempio CONCAVE <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_concave


polynomial: setup example_polynomial
	@echo ""
	@echo ">>> Avvio esempio POLYNOMIAL <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_polynomial


bunny: setup example_bunny
	@echo ""
	@echo ">>> Avvio esempio BUNNY <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_bunny


# =====================================================================
# SEZIONE 16 — ESECUZIONE
#
# Alias per la compilazione di tutti gli esempi.
#
# =====================================================================

run: examples


# =====================================================================
# SEZIONE 17 — PULIZIA
# =====================================================================

clean:
	@echo "Pulizia cartelle di build..."
	rm -rf $(OBJ_DIR) $(MOD_DIR)
	rm -f $(EXAMPLES)
	@echo "Rimosso: $(OBJ_DIR)/"
	@echo "Rimosso: $(MOD_DIR)/"
	@echo "Rimossi tutti gli eseguibili."


# =====================================================================
# SEZIONE 18 — PULIZIA COMPLETA
# =====================================================================

distclean: clean
	@echo "Pulizia completa (distclean)..."
	rm -f punti.txt *.dat *.out
	@echo "Progetto riportato allo stato sorgente."


# =====================================================================
# SEZIONE 19 — INFORMAZIONI
# =====================================================================

info:
	@echo ""
	@echo "======================================"
	@echo "          OptimalPolyCuba3D"
	@echo "======================================"
	@echo " Esempi:"
	@echo "   example_convex"
	@echo "   example_concave"
	@echo "   example_polynomial"
	@echo "   example_bunny"
	@echo ""
	@echo " Compilatore : $(F90)"
	@echo " Flag        : $(F90FLAGS)"
	@echo " Sorgenti    : $(SRC_DIR)/"
	@echo " Esempi      : $(EXAMPLE_DIR)/"
	@echo " Oggetti     : $(OBJ_DIR)/"
	@echo " Moduli .mod : $(MOD_DIR)/"
	@echo "======================================"
	@echo ""