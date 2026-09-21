# ==============================================================================
#
# Makefile
#
# OptimalPolyCuba4D Tensor
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

EXAMPLES = $(patsubst $(EXAMPLE_DIR)/example_%.f90,%,$(wildcard $(EXAMPLE_DIR)/example_*.f90))
TARGETS  = $(addprefix $(EXAMPLE_DIR)/example_,$(EXAMPLES))

# Numero di core logici.
# macOS: sysctl
# Linux: nproc

NCPU := $(shell sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)

# Abilitazione OpenMP:
#   make OMP=1   -> OpenMP abilitato
#   make OMP=0   -> OpenMP disabilitato

OMP ?= 1

THREADS   ?= $(NCPU)
STACKSIZE ?= 256M

# Esempio usato da "make scaling".

SCALING_EXAMPLE ?= convex


# ==============================================================================
# SEZIONE 2 — FLAG DI COMPILAZIONE
# ==============================================================================

ifeq ($(OMP),1)

OMPFLAGS = -fopenmp

else

OMPFLAGS =

endif


FFLAGS = \
-O3 \
-march=native \
-flto \
-funroll-loops \
-ffree-line-length-none \
-Wall \
-Wextra \
$(OMPFLAGS) \
-J$(MODULE_DIR) \
-I$(MODULE_DIR)


LDFLAGS = \
-O3 \
-march=native \
-flto \
$(OMPFLAGS)


# ==============================================================================
# SEZIONE 3 — FILE OGGETTO
# ==============================================================================

MODULE_OBJECTS = \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/TriangleQuadratureGJ.o \
$(BUILD_DIR)/TriangleQuadratureDunavant.o \
$(BUILD_DIR)/PrismQuadrature.o \
$(BUILD_DIR)/CubatureTensor.o \
$(BUILD_DIR)/OPC4D_Parallel_Tensor.o \
$(BUILD_DIR)/PolyhedronMesh.o


# ==============================================================================
# SEZIONE 4 — TARGET PRINCIPALI
# ==============================================================================

.PHONY: all examples $(EXAMPLES) run parallel scaling clean distclean info


all: examples


examples: $(TARGETS)
	@echo ""
	@echo "=============================================="
	@echo "  Tutti gli esempi Tensor 4D sono compilati."
	@echo "=============================================="
	@echo ""
	@echo "  Eseguibili:"
	@for e in $(EXAMPLES); do echo "    ./$(EXAMPLE_DIR)/example_$$e"; done
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


$(BUILD_DIR)/TriangleQuadratureGJ.o: \
$(SRC_DIR)/TriangleQuadratureGJ.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o


$(BUILD_DIR)/TriangleQuadratureDunavant.o: \
$(SRC_DIR)/TriangleQuadratureDunavant.f90 \
$(BUILD_DIR)/TypesDef.o


$(BUILD_DIR)/PrismQuadrature.o: \
$(SRC_DIR)/PrismQuadrature.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/TriangleQuadratureGJ.o \
$(BUILD_DIR)/TriangleQuadratureDunavant.o


$(BUILD_DIR)/CubatureTensor.o: \
$(SRC_DIR)/CubatureTensor.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/PrismQuadrature.o


$(BUILD_DIR)/OPC4D_Parallel_Tensor.o: \
$(SRC_DIR)/OPC4D_Parallel_Tensor.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/ReferenceFunctions.o \
$(BUILD_DIR)/CubatureTensor.o \
$(BUILD_DIR)/PrismQuadrature.o


$(BUILD_DIR)/PolyhedronMesh.o: \
$(SRC_DIR)/PolyhedronMesh.f90 \
$(BUILD_DIR)/TypesDef.o


# ==============================================================================
# SEZIONE 8 — COMPILAZIONE DEGLI ESEMPI
# ==============================================================================

$(BUILD_DIR)/example_%.o: \
$(EXAMPLE_DIR)/example_%.f90 \
$(BUILD_DIR)/TypesDef.o \
$(BUILD_DIR)/PolyhedronMesh.o \
$(BUILD_DIR)/OPC4D_Parallel_Tensor.o
	@echo "  Compiling  $<"
	$(FC) $(FFLAGS) -c $< -o $@


# ==============================================================================
# SEZIONE 9 — LINKING
# ==============================================================================

$(TARGETS): $(EXAMPLE_DIR)/example_%: \
$(MODULE_OBJECTS) \
$(BUILD_DIR)/example_%.o
	@echo "  Linking    $@"
	$(FC) $(LDFLAGS) $^ -o $@


# ==============================================================================
# SEZIONE 10 — ESECUZIONE DEGLI ESEMPI
# ==============================================================================

# make convex | concave | polynomial | ... THREADS=n
#
# Gli esempi vengono eseguiti dalla directory examples/
# in modo che eventuali file .dat relativi vengano trovati correttamente.


$(EXAMPLES): %: $(EXAMPLE_DIR)/example_%
	@echo ""
	@echo ">>> Avvio esempio $(shell echo $@ | tr a-z A-Z) (OMP=$(OMP), thread: $(THREADS)) <<<"
	@echo ""
	cd $(EXAMPLE_DIR) && \
	OMP_NUM_THREADS=$(THREADS) \
	OMP_STACKSIZE=$(STACKSIZE) \
	./example_$@


run: examples
	@for e in $(EXAMPLES); do \
		$(MAKE) --no-print-directory $$e THREADS=$(THREADS); \
	done


# ==============================================================================
# SEZIONE 11 — BUILD PARALLELA E SCALABILITA'
# ==============================================================================

# Compilazione dei sorgenti su piu' processi.

parallel:
	@echo ""
	@echo "Compilazione parallela con $(NCPU) core..."
	@echo ""
	$(MAKE) -j$(NCPU) examples


# ------------------------------------------------------------------------------
# Scalabilita' OpenMP
# ------------------------------------------------------------------------------
#
# Esegue l'esempio SCALING_EXAMPLE con:
#
#   1, 2, 4, 8, ..., NCPU
#
# thread.
#
# Richiede OMP=1 e costrutti !$OMP nei sorgenti.
#

scaling: $(EXAMPLE_DIR)/example_$(SCALING_EXAMPLE)
	@echo ""
	@echo "Scalabilita' OpenMP - esempio: $(SCALING_EXAMPLE)"
	@echo "OpenMP=$(OMP), core disponibili=$(NCPU)"
	@echo ""

	@t=1; \
	while [ $$t -le $(NCPU) ]; do \
		printf "  thread = %3d   " $$t; \
		( cd $(EXAMPLE_DIR) && \
		OMP_NUM_THREADS=$$t \
		OMP_STACKSIZE=$(STACKSIZE) \
		./example_$(SCALING_EXAMPLE) | grep "Tempo di calcolo" ); \
		t=$$((t * 2)); \
	done

	@echo ""


# ==============================================================================
# SEZIONE 12 — PULIZIA
# ==============================================================================

clean:
	@echo "Pulizia della build..."
	rm -rf $(BUILD_DIR)
	rm -rf $(MODULE_DIR)
	rm -f $(TARGETS)


distclean: clean
	@echo "Pulizia completa..."
	rm -f punti.txt
	rm -f *.out
	@echo "Progetto riportato allo stato sorgente."


# ==============================================================================
# SEZIONE 13 — INFORMAZIONI
# ==============================================================================

info:
	@echo ""
	@echo "=============================================="
	@echo "          OptimalPolyCuba4D Tensor"
	@echo "=============================================="
	@echo ""
	@echo " Compilatore : $(FC)"
	@echo " Ottimizzaz. : -O3 -flto -funroll-loops"
	@echo " OpenMP      : $(if $(filter 1,$(OMP)),abilitato,disabilitato)"
	@echo " Core CPU    : $(NCPU)"
	@echo " Thread run  : $(THREADS)"
	@echo " Sorgenti    : $(SRC_DIR)/"
	@echo " Esempi      : $(EXAMPLE_DIR)/"
	@echo " Oggetti     : $(BUILD_DIR)/"
	@echo " Moduli      : $(MODULE_DIR)/"
	@echo ""
	@echo " Esempi rilevati:"
	@for e in $(EXAMPLES); do echo "   make $$e"; done
	@echo ""
	@echo " Build completa:"
	@echo "   make"
	@echo ""
	@echo " Build parallela:"
	@echo "   make parallel"
	@echo ""
	@echo " Esecuzione con un dato numero di thread:"
	@echo "   make convex THREADS=4"
	@echo ""
	@echo " Scalabilita' OpenMP:"
	@echo "   make scaling"
	@echo ""
	@echo " Build seriale di riferimento:"
	@echo "   make clean && make OMP=0"
	@echo ""
	@echo " Pulizia:"
	@echo "   make clean"
	@echo "   make distclean"
	@echo ""
	@echo "=============================================="