Compiler = g++
CPPFLAGS = -O3 -I /n/home/bstefan/m/boost_1_57_0/boost_1_57_0/ `root-config --cflags --glibs`


#-I /net/opt/root/root-5.34.10-wheezy/include/root/

#this doesn't get used for some reason...
#%.o: %.cc 
#	$(Compiler)  -o -c $@ $< $(decay_flags)

create_root: create_root.o
	$(Compiler) $(CPPFLAGS) -o $@ $^ 

clean:
	rm -f *.o
	touch *.cc



