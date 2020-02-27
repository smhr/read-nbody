#!/usr/bin/make

#main building variables
DSRC    = .
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = gfortran
DEBUGFLAG = -g -fbounds-check -ffree-line-length-none -Wall
OPTSC   = -c -J mod -O3 -ffast-math
OPTSC   += ${DEBUGFLAG}
OPTSL   =  -J mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)READ-NBODY-TAIL: $(MKDIRS) $(DOBJ)read-nbody-tail.o
	@rm -f $(filter-out $(DOBJ)read-nbody-tail.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) READ-NBODY-TAIL
$(DEXE)READ-NBODY-RT: $(MKDIRS) $(DOBJ)read-nbody-rt.o
	@rm -f $(filter-out $(DOBJ)read-nbody-rt.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) READ-NBODY-RT
$(DEXE)READ-NBODY-2RT: $(MKDIRS) $(DOBJ)read-nbody-2rt.o
	@rm -f $(filter-out $(DOBJ)read-nbody-2rt.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) READ-NBODY-2RT

#compiling rules
$(DOBJ)read-nbody-tail.o: ./read-nbody-tail.f90 \
	$(DOBJ)kdtree2.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)read-nbody-rt.o: ./read-nbody-rt.f90 \
	$(DOBJ)kdtree2.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)read-nbody-2rt.o: ./read-nbody-2rt.f90 \
	$(DOBJ)kdtree2.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)kdtree2.o: ./kdtree2.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
