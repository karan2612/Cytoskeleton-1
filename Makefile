%.out: %.cpp
	g++ $^ -lgr -lgr3 \
	-L/usr/local/gr/lib/ \
	-o $@

%.run: %.out 
	DYLD_LIBRARY_PATH=/usr/local/gr/lib/ \
	./$^
