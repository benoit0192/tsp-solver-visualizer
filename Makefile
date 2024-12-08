
all: solver map-renderer

solver:
	$(MAKE) -f Makefile.tsp_solver

solver-run:
	$(MAKE) -f Makefile.tsp_solver run \
					$(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
map-renderer:
	$(MAKE) -f ./Makefile.tsp_map_renderer

map-renderer-run:
	$(MAKE) -f ./Makefile.tsp_map_renderer run \
					$(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))

clean:
	$(MAKE) -f Makefile.tsp_solver       clean
	$(MAKE) -f Makefile.tsp_map_renderer clean

.PHONY: all clean solver solver-run map-renderer map-renderer-run
