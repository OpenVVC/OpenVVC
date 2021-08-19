%:
	@[ -d build/$@/ ] && (cd build/$@/; $(MAKE) --no-print-directory $(filter-out $@,$(MAKECMDGOALS))) || :
