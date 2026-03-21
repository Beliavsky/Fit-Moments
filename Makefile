MAKEFLAGS = --no-print-directory
SHELL = cmd.exe
FC = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private

# Executables
EXES = xrandom.exe xdist_moments.exe xfit_moments_csv_table.exe

all: $(EXES)

# Object files for modules
kind.o: kind.f90
	$(FC) $(FFLAGS) -c kind.f90

random.o: random.f90 kind.o
	$(FC) $(FFLAGS) -c random.f90

dist_moments.o: dist_moments.f90 kind.o
	$(FC) $(FFLAGS) -c dist_moments.f90

# Object files for programs
xrandom.o: xrandom.f90 random.o
	$(FC) $(FFLAGS) -c xrandom.f90

xdist_moments.o: xdist_moments.f90 dist_moments.o kind.o
	$(FC) $(FFLAGS) -c xdist_moments.f90

xfit_moments_csv_table.o: xfit_moments_csv_table.f90 dist_moments.o kind.o
	$(FC) $(FFLAGS) -c xfit_moments_csv_table.f90

# Final executables
xrandom.exe: xrandom.o random.o kind.o
	$(FC) $(FFLAGS) xrandom.o random.o kind.o -o xrandom.exe

xdist_moments.exe: xdist_moments.o dist_moments.o kind.o
	$(FC) $(FFLAGS) xdist_moments.o dist_moments.o kind.o -o xdist_moments.exe

xfit_moments_csv_table.exe: xfit_moments_csv_table.o dist_moments.o kind.o
	$(FC) $(FFLAGS) xfit_moments_csv_table.o dist_moments.o kind.o -o xfit_moments_csv_table.exe

.PHONY: clean
clean:
	del /f /q *.o *.mod $(EXES) 2>nul
