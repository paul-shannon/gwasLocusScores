default:
	@echo targets:  all [roxy install test]

all: roxy install test

roxy:
	R -e "devtools::document()"

install:
	R CMD INSTALL . --no-test-load

test:
	(cd inst/unitTests; make)

