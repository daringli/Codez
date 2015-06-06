Compiler = g++

ROOT_includes =-I /net/opt/root/root-5.34.10-wheezy/include/root
ROOT_libs =-L/net/opt/root/root-5.34.10-wheezy/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -Wl,-rpath,/net/opt/root/root-5.34.10-wheezy/lib/root -lm -ldl -rdynamic



eg_objects = globals.o model.o mass_model.o codex_gamma_model.o codex_particle_model.o integrate.o deexcite.o  deexcite_table.o decay_chain.o lorentz_boost.o print_gunfile.o elements.o quasi_elastic.o can_be_run.o

chart_objects = globals.o model.o mass_model.o codex_gamma_model.o codex_particle_model.o integrate.o deexcite.o  deexcite_table.o decay_chain.o separation_energy.o decay_name.o particle_probability.o print_gunfile.o lorentz_boost.o elements.o run_talys.o file_names.o talys_probability.o can_be_run.o

nuclist_objects = globals.o model.o mass_model.o codex_gamma_model.o codex_particle_model.o integrate.o deexcite.o  deexcite_table.o decay_chain.o separation_energy.o decay_name.o particle_probability.o print_gunfile.o lorentz_boost.o elements.o run_talys.o file_names.o talys_probability.o can_be_run.o

spectra_objects = globals.o model.o mass_model.o codex_gamma_model.o codex_particle_model.o integrate.o deexcite.o  deexcite_table.o decay_chain.o separation_energy.o decay_name.o particle_probability.o print_gunfile.o lorentz_boost.o elements.o run_talys.o file_names.o talys_probability.o can_be_run.o

rho_objects = mass_model.o globals.o model.o codex_particle_model.o can_be_run.o integrate.o

pot_objects = mass_model.o globals.o model.o codex_particle_model.o can_be_run.o integrate.o

trans_objects = mass_model.o globals.o model.o codex_particle_model.o can_be_run.o integrate.o

sep_objects = mass_model.o globals.o

out2gun_objects = elements.o

out2root_objects=

out-spectra2root_objects=


deexcite_objects = globals.o model.o mass_model.o codex_gamma_model.o codex_particle_model.o integrate.o deexcite.o  deexcite_table.o decay_chain.o lorentz_boost.o print_gunfile.o elements.o quasi_elastic.o can_be_run.o print_products.o

quasi_objects = globals.o model.o mass_model.o codex_gamma_model.o codex_particle_model.o integrate.o deexcite.o  deexcite_table.o decay_chain.o lorentz_boost.o print_gunfile.o elements.o quasi_elastic.o can_be_run.o

#since it for some reason uses implicit rule over the one defined below
CXXFLAGS = -O3 -I /net/opt/root/root-5.34.10-wheezy/include/root -I /n/home/bstefan/m/boost_1_57_0/boost_1_57_0/  -L /n/home/bstefan/m/build-boost/boost/bin.v2/libs/program_options/build/gcc-4.7/release/link-static/threading-multi/ 
Boost_po=-l boost_program_options

eg: $(eg_objects) event_generator.o
	$(Compiler) -o $@ $^ 

chart: $(chart_objects) map_chart.o
	$(Compiler) -o $@ $^ 

rho: $(rho_objects) level_density.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

pot: $(pot_objects) pot.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

trans: $(trans_objects) trans.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

sep: $(sep_objects) sep_energy.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)


all: eg chart

nuclist: $(nuclist_objects) nuclist.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

spectra: $(spectra_objects) spectra.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

quasi: $(quasi_objects) quasi.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

deexcite: $(deexcite_objects) deexcite_list.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

out2gun: $(out2gun_objects) out2gun.o
	$(Compiler) $(CXXFLAGS) -o $@ $^ $(Boost_po)

out2root: $(out2root_objects) out2root.o
	$(Compiler)  $(CXXFLAGS) -g -o out2root out2root.cc $(ROOT_libs)  $(Boost_po)

out-spectra2root: $(out2root_objects) out-spectra2root.o
	$(Compiler)  $(CXXFLAGS) -g -o out-spectra2root out-spectra2root.cc $(ROOT_libs)  $(Boost_po)

clean:
	rm -f *.o


