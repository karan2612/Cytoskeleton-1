%.out: %.C
	g++ $^ -lgr -lgr3 \
	-L/Users/q/miniconda2/lib/python2.7/site-packages/gr/ \
	-L/Users/q/miniconda2/lib/python2.7/site-packages/gr3/ \
	-o $@

%.run: %.out
	DYLD_LIBRARY_PATH=/Users/q/miniconda2/lib/python2.7/site-packages/gr/:/Users/q/miniconda2/lib/python2.7/site-packages/gr3/ \
	./$^

