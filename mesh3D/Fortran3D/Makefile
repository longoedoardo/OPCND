# ==============================================================================
#
# Makefile
#
# OptimalPolyCuba3D
#
# Autore: Edoardo Longo
# Data: 18/09/2026
#
# ==============================================================================


# ==============================================================================
# SEZIONE 1 — CONFIGURAZIONE
# ==============================================================================

FC = gfortran

TARGETS = \
	example_convex \
	example_concave \
	example_polynomial \
	example_bunny

BUILD_DIR   = build
MODULE_DIR  = modules_build
SRC_DIR     = src
EXAMPLE_DIR = examples


# ==============================================================================
# SEZIONE 2 — FLAG DI COMPILAZIONE
# ==============================================================================

FFLAGS = \
	-O3 \
	-ffree-line-length-none \
	-Wall \
	-Wextra \
	-J$(MODULE_DIR) \
	-I$(MODULE_DIR)

LDFLAGS = -O3


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
	$(BUILD_DIR)/OPC3D.o

EXAMPLE_OBJECTS = \
	$(BUILD_DIR)/example_convex.o \
	$(BUILD_DIR)/example_concave.o \
	$(BUILD_DIR)/example_polynomial.o \
	$(BUILD_DIR)/example_bunny.o


# ==============================================================================
# SEZIONE 4 — TARGET PRINCIPALI
# ==============================================================================

.PHONY: all examples \
	convex concave polynomial bunny \
	run clean distclean info


all: examples


examples: $(TARGETS)
	@echo ""
	@echo "=============================================="
	@echo "  Tutti gli esempi sono stati compilati."
	@echo "=============================================="
	@echo ""
	@echo "  Eseguibili:"
	@echo "    ./example_convex"
	@echo "    ./example_concave"
	@echo "    ./example_polynomial"
	@echo "    ./example_bunny"
	@echo ""


# ==============================================================================
# SEZIONE 5 — CREAZIONE DIRECTORY
# ==============================================================================

$(BUILD_DIR) $(MODULE_DIR):
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


$(BUILD_DIR)/OPC3D.o: \
	$(SRC_DIR)/OPC3D.f90 \
	$(BUILD_DIR)/TypesDef.o \
	$(BUILD_DIR)/ReferenceFunctions.o \
	$(BUILD_DIR)/CubatureFunctions.o


# ==============================================================================
# SEZIONE 8 — COMPILAZIONE DEGLI ESEMPI
# ==============================================================================

$(BUILD_DIR)/example_convex.o: \
	$(EXAMPLE_DIR)/example_convex.f90 \
	$(BUILD_DIR)/TypesDef.o \
	$(BUILD_DIR)/PolyhedronMesh.o \
	$(BUILD_DIR)/OPC3D.o

$(BUILD_DIR)/example_concave.o: \
	$(EXAMPLE_DIR)/example_concave.f90 \
	$(BUILD_DIR)/TypesDef.o \
	$(BUILD_DIR)/PolyhedronMesh.o \
	$(BUILD_DIR)/OPC3D.o

$(BUILD_DIR)/example_polynomial.o: \
	$(EXAMPLE_DIR)/example_polynomial.f90 \
	$(BUILD_DIR)/TypesDef.o \
	$(BUILD_DIR)/PolyhedronMesh.o \
	$(BUILD_DIR)/OPC3D.o

$(BUILD_DIR)/example_bunny.o: \
	$(EXAMPLE_DIR)/example_bunny.f90 \
	$(BUILD_DIR)/TypesDef.o \
	$(BUILD_DIR)/PolyhedronMesh.o \
	$(BUILD_DIR)/OPC3D.o


$(BUILD_DIR)/example_%.o: $(EXAMPLE_DIR)/example_%.f90
	@echo "  Compiling  $<"
	$(FC) $(FFLAGS) -c $< -o $@


# ==============================================================================
# SEZIONE 9 — LINKING
# ==============================================================================

example_convex: \
	$(MODULE_OBJECTS) \
	$(BUILD_DIR)/example_convex.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@


example_concave: \
	$(MODULE_OBJECTS) \
	$(BUILD_DIR)/example_concave.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@


example_polynomial: \
	$(MODULE_OBJECTS) \
	$(BUILD_DIR)/example_polynomial.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@


example_bunny: \
	$(MODULE_OBJECTS) \
	$(BUILD_DIR)/example_bunny.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@


# ==============================================================================
# SEZIONE 10 — ESECUZIONE DEGLI ESEMPI
# ==============================================================================

convex: example_convex
	@echo ""
	@echo ">>> Avvio esempio CONVEX <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_convex


concave: example_concave
	@echo ""
	@echo ">>> Avvio esempio CONCAVE <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_concave


polynomial: example_polynomial
	@echo ""
	@echo ">>> Avvio esempio POLYNOMIAL <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_polynomial


bunny: example_bunny
	@echo ""
	@echo ">>> Avvio esempio BUNNY <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && ../example_bunny


run: examples


# ==============================================================================
# SEZIONE 11 — PULIZIA
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
# SEZIONE 12 — INFORMAZIONI
# ==============================================================================

info:
	@echo ""
	@echo "=============================================="
	@echo "             OptimalPolyCuba3D"
	@echo "=============================================="
	@echo ""
	@echo " Compilatore : $(FC)"
	@echo " Ottimizzaz. : -O3"
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
	@echo "   make -j"
	@echo ""
	@echo " Pulizia:"
	@echo "   make clean"
	@echo "   make distclean"
	@echo ""
	@echo "=============================================="