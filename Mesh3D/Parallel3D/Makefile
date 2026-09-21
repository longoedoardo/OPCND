# ==============================================================================
#
# Makefile
#
# OptimalPolyCuba3D
#
# Autore: Edoardo Longo
# Data: 19/09/2026
#
# ==============================================================================


# ==============================================================================
# SEZIONE 1 — CONFIGURAZIONE
# ==============================================================================

FC = gfortran

BUILD_DIR   = build
MODULE_DIR  = modules_build
SRC_DIR     = src
EXAMPLE_DIR = examples

TARGETS = \
$(EXAMPLE_DIR)/example_convex \
$(EXAMPLE_DIR)/example_concave \
$(EXAMPLE_DIR)/example_polynomial \
$(EXAMPLE_DIR)/example_bunny

# Numero di core logici disponibili su macOS.

NCPU := $(shell sysctl -n hw.ncpu)


# ==============================================================================
# SEZIONE 2 — FLAG DI COMPILAZIONE
# ==============================================================================

FFLAGS = \
-O3 \
-march=native \
-flto \
-funroll-loops \
-ffree-line-length-none \
-Wall \
-Wextra \
-fopenmp \
-J$(MODULE_DIR) \
-I$(MODULE_DIR)

LDFLAGS = \
-O3 \
-march=native \
-flto \
-fopenmp


# ==============================================================================
# SEZIONE 3 — FILE OGGETTO
# ==============================================================================

MODULE_OBJECTS = \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/PolyhedronMesh.o \
$(BUILD_DIR)/TriangleQuadratureGJ.o \
$(BUILD_DIR)/TriangleQuadratureDunavant.o \
$(BUILD_DIR)/CubatureFunctions.o \
$(BUILD_DIR)/OPC3D_Parallel.o


# ==============================================================================
# SEZIONE 4 — TARGET PRINCIPALI
# ==============================================================================

.PHONY: all examples convex concave polynomial bunny run parallel clean distclean info

all: examples

examples: $(TARGETS)
	@echo ""
	@echo "=============================================="
	@echo "  Tutti gli esempi sono stati compilati."
	@echo "=============================================="
	@echo ""
	@echo "  Eseguibili:"
	@echo "    ./examples/example_convex"
	@echo "    ./examples/example_concave"
	@echo "    ./examples/example_polynomial"
	@echo "    ./examples/example_bunny"
	@echo ""


# ==============================================================================
# SEZIONE 5 — CREAZIONE DIRECTORY
# ==============================================================================

$(BUILD_DIR):
	@mkdir -p $@

$(MODULE_DIR):
	@mkdir -p $@


# ==============================================================================
# SEZIONE 6 — COMPILAZIONE DEI MODULI FORTRAN
# ==============================================================================

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR) $(MODULE_DIR)
	@echo "  Compiling  $<"
	$(FC) $(FFLAGS) -c $< -o $@


# ==============================================================================
# SEZIONE 7 — DIPENDENZE DEI MODULI
# ==============================================================================

$(BUILD_DIR)/TypesDef.o: \
$(SRC_DIR)/TypesDef.f90

$(BUILD_DIR)/ReferenceFunctions.o: \
$(SRC_DIR)/ReferenceFunctions.f90 \
$(BUILD_DIR)/TypesDef.o

$(BUILD_DIR)/PolyhedronMesh.o: \
$(SRC_DIR)/PolyhedronMesh.f90 \
$(BUILD_DIR)/TypesDef.o

$(BUILD_DIR)/TriangleQuadratureGJ.o: \
$(SRC_DIR)/TriangleQuadratureGJ.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o

$(BUILD_DIR)/TriangleQuadratureDunavant.o: \
$(SRC_DIR)/TriangleQuadratureDunavant.f90 \
$(BUILD_DIR)/TypesDef.o

$(BUILD_DIR)/CubatureFunctions.o: \
$(SRC_DIR)/CubatureFunctions.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/TriangleQuadratureGJ.o \
$(BUILD_DIR)/TriangleQuadratureDunavant.o

$(BUILD_DIR)/OPC3D_Parallel.o: \
$(SRC_DIR)/OPC3D_Parallel.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/CubatureFunctions.o


# ==============================================================================
# SEZIONE 8 — COMPILAZIONE DEGLI ESEMPI
# ==============================================================================

$(BUILD_DIR)/example_%.o: \
$(EXAMPLE_DIR)/example_%.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/PolyhedronMesh.o \
$(BUILD_DIR)/OPC3D_Parallel.o | $(BUILD_DIR) $(MODULE_DIR)
	@echo "  Compiling  $<"
	$(FC) $(FFLAGS) -c $< -o $@


# ==============================================================================
# SEZIONE 9 — LINKING
# ==============================================================================

$(EXAMPLE_DIR)/example_convex: \
$(MODULE_OBJECTS) \
$(BUILD_DIR)/example_convex.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@

$(EXAMPLE_DIR)/example_concave: \
$(MODULE_OBJECTS) \
$(BUILD_DIR)/example_concave.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@

$(EXAMPLE_DIR)/example_polynomial: \
$(MODULE_OBJECTS) \
$(BUILD_DIR)/example_polynomial.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@

$(EXAMPLE_DIR)/example_bunny: \
$(MODULE_OBJECTS) \
$(BUILD_DIR)/example_bunny.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@


# ==============================================================================
# SEZIONE 10 — ESECUZIONE DEGLI ESEMPI
# ==============================================================================

convex: $(EXAMPLE_DIR)/example_convex
	@echo ""
	@echo ">>> Avvio esempio CONVEX <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ./example_convex

concave: $(EXAMPLE_DIR)/example_concave
	@echo ""
	@echo ">>> Avvio esempio CONCAVE <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ./example_concave

polynomial: $(EXAMPLE_DIR)/example_polynomial
	@echo ""
	@echo ">>> Avvio esempio POLYNOMIAL <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ./example_polynomial

bunny: $(EXAMPLE_DIR)/example_bunny
	@echo ""
	@echo ">>> Avvio esempio BUNNY <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ./example_bunny

run: examples


# ==============================================================================
# SEZIONE 11 — BUILD PARALLELA
# ==============================================================================

parallel:
	@echo ""
	@echo "Compilazione parallela con $(NCPU) core..."
	@echo ""
	$(MAKE) -j$(NCPU) examples


# ==============================================================================
# SEZIONE 12 — PULIZIA
# ==============================================================================

clean:
	@echo "Pulizia della build..."
	rm -rf $(BUILD_DIR) $(MODULE_DIR)
	rm -f $(TARGETS)

distclean: clean
	@echo "Pulizia completa..."
	rm -f punti.txt *.dat *.out
	@echo "Progetto riportato allo stato sorgente."


# ==============================================================================
# SEZIONE 13 — INFORMAZIONI
# ==============================================================================

info:
	@echo ""
	@echo "=============================================="
	@echo "             OptimalPolyCuba3D"
	@echo "=============================================="
	@echo ""
	@echo " Compilatore : $(FC)"
	@echo " Ottimizzaz. : -O3 -flto -funroll-loops"
	@echo " OpenMP      : abilitato"
	@echo " Core CPU    : $(NCPU)"
	@echo " Sorgenti    : $(SRC_DIR)/"
	@echo " Esempi      : $(EXAMPLE_DIR)/"
	@echo " Oggetti     : $(BUILD_DIR)/"
	@echo " Moduli      : $(MODULE_DIR)/"
	@echo ""
	@echo " Esempi:"
	@echo "   make convex"
	@echo "   make concave"
	@echo "   make polynomial"
	@echo "   make bunny"
	@echo ""
	@echo " Build completa:"
	@echo "   make"
	@echo ""
	@echo " Build parallela:"
	@echo "   make parallel"
	@echo ""
	@echo " Pulizia:"
	@echo "   make clean"
	@echo "   make distclean"
	@echo ""
	@echo "=============================================="
	@echo ""