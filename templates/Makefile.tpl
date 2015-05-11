scf:
		gfortran scf.f potential.f -o scf

clean:
		rm scf
		rm SNAP*
		rm SCFCPU SCFCEN SCFORB SCFOUT SCFLOG
