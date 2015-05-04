
CC=/usr/local/bin/gcc
FUNGIMOL_INCLUDE=/home/tim/fungimol/include
FUNGIMOL_MEMORY_DIR=/home/tim/fungimol/Util/debug
# .DELETE_ON_ERROR means that make will delete files it was trying to
# generate, if the command that was generating them returns an error
# status.  This is always a good thing but it is incompatible with
# some historical behavior that apparently somebody wanted.  See the
# gnu-make info file and search for ".DELETE_ON_ERROR".
.PHONY: .DELETE_ON_ERROR
.DELETE_ON_ERROR:
# $(BRENNERMD) is the Fortran version, so we can compare outputs when testing.
BRENNERMD:=./brennermd
ifeq ($(wildcard $(BRENNERMD)),)
   BRENNERMD:=brennermd
endif
RASAPIDIR = ../rasmolapi
INCLUDES = -I$(RASAPIDIR)
# -funroll-loops makes exactly no difference, so leave it out.
# 5.24 5.22 5.22 5.22 elapsed without -funroll-loops, 
# 5.63 5.22 5.22 5.22 elapsed with -funroll-loops.  
# 
# -mcpu=pentiumpro doesn't help either.
# 
# On my battery-powered portable, there are large variations in the
# time things take, so don't use it for timing runs.  The same run can
# take 7 or 8 or 9 or 10 seconds.  In contrast, times from a 600 Mhz
# Gateway Pentium plugged into the wall are very reproducible after
# the first run.  The portable still has large variations when it's
# plugged into the wall. 
#
# Tim Freeman 24 Jul 2000
#
C_FLAGS = -Wall -Wcast-qual -Wno-unused -W -ffast-math -fno-exceptions -Werror -Winline
ifeq ($(DEBUG),)
   C_FLAGS += -fomit-frame-pointer -DNDEBUG -O3 
else
   ifeq ($(wildcard /usr/include/mcheck.h),)
      C_FLAGS += -g
   else
      C_FLAGS += -g -DMCHECK
   endif
# Comment in next line when I'm using the debugger to get information
# about performance.
#   C_FLAGS += -DNDEBUG -O3
endif
ifneq ($(DEBUGPRINT),)
   C_FLAGS += -DDEBUGPRINT
endif
ifneq ($(FUNGIMOL),)
   C_FLAGS += -DFUNGIMOL -DDOUBLE_PRECISION
   # The next line is required so we get fungimol's include directory
   # on the path.
   FUNGIMOL_MEMORY:=yes
   INFINITE_CUBE:=yes
   # Fungimol is only one precision, and Brennerc needs double
   # precision sometimes, so if we're using both together we need both
   # to always be double precision.  Fungimol needs to change to make
   # this better.
   DOUBLE:=yes
endif
# Whether to use Fungimol's memory debugging.
ifneq ($(FUNGIMOL_MEMORY),)
   C_FLAGS += -DFUNGIMOL_MEMORY -I$(FUNGIMOL_INCLUDE)
endif
ifneq ($(INFINITE_CUBE),)
   C_FLAGS += -DINFINITE_CUBE
endif
ifneq ($(RUNGEKUTTA),)
   C_FLAGS += -DRUNGEKUTTA
endif
ifneq ($(N2_ASYM),)
   C_FLAGS += -DN2_ASYM
endif

# With comparable code, if SINGLE was defined a timing run took 4.38
# seconds and otherwise it took 4.88 seconds, about a 10% speedup.

# The convention is that a "Double" in the source is a double
# precision floating point number that has to be double precision so
# things are accurate enough.  In constrast, a "double" in the source
# has to be a double to make the compiler happy.  So the code should
# compile, run, and give wobbly answers if we have -DDouble=float.
# Once I'm sure what type I want everything to have, we can forget
# about some of these distinctions.  

# Compile with DOUBLE=yes to get double precision.  The default is
# single precision, since that seems to get good results now.
ifeq ($(DOUBLE),)
   SINGLE:=yes
endif
ifeq ($(SINGLE),)
   C_FLAGS += -DFloat=double -DDouble=double
else
# We really do get wrong answers if we do all the math in single
# precision.  It's a chaotic system, so I can't say that I'm sure the
# wrong answers are wrong in an important way, but it does seem
# worthwhile to avoid questionable results.  
   C_FLAGS += -DFloat=float -DDouble=double
#   C_FLAGS += -DFloat=float -DDouble=float
endif
LIBS = -L/usr/X11R6/lib -L$(RASAPIDIR) -lrasmol -lXi -lXext -lX11 -lm

# Have = here, not :=, since we may change OBJDIR below.
ifeq ($(FUNGIMOL),)
  OBJDIR = .
else
   OBJDIR = fungimol
   JUNK:=$(shell mkdir $(OBJDIR) >& /dev/null)
endif
OBJS = $(OBJDIR)/AtomPairInfoState.o $(OBJDIR)/api.o \
       $(OBJDIR)/caguts.o $(OBJDIR)/rb2.o \
       $(OBJDIR)/ljguts.o $(OBJDIR)/pibond.o \
       $(OBJDIR)/sili_germ.o $(OBJDIR)/bcuint.o $(OBJDIR)/radic.o \
       $(OBJDIR)/mtable.o $(OBJDIR)/vector.o \
       $(OBJDIR)/spgch.o $(OBJDIR)/myassert.o \
       $(OBJDIR)/expand.o $(OBJDIR)/xalloc.o $(OBJDIR)/rmax.o $(OBJDIR)/kt.o \
       $(OBJDIR)/rangetree.o $(OBJDIR)/xmass.o \
       $(OBJDIR)/bren2.o $(OBJDIR)/rlist.o $(OBJDIR)/xmms.o \
       $(OBJDIR)/CCNeighborState.o $(OBJDIR)/rslj.o

ifeq ($(FUNGIMOL),)
   # Sources that don't need to be in the Fungimol plugin.
   OBJS += $(OBJDIR)/ccapi.o $(OBJDIR)/safe_fgets.o \
     $(OBJDIR)/minimize.o $(OBJDIR)/vscale.o
   ifneq ($(FUNGIMOL_MEMORY),) 
      # Pick up just the objects required for Fungimol leak detection
      # and memory corruption checking
      OBJS += $(FUNGIMOL_MEMORY_DIR)/MemoryUtil.o \
              $(FUNGIMOL_MEMORY_DIR)/Refcount.o
   endif
   # Where to find the inter*.d files.
   DATA_HOME=.
else
   # Fungimol-only sources.
   OBJS += $(OBJDIR)/BrennerAtom.o $(OBJDIR)/State.o
   # Where to find the inter*.d files.
   DATA_HOME=/usr/share/brennermd
endif
# The quotes are C syntax, since we need DATA_HOME to be a string.
C_FLAGS += -DDATA_HOME=\"$(DATA_HOME)\"

$(OBJDIR)/%.o: %.c
	$(CC) $(INCLUDES) $(C_FLAGS) -c $*.c -o $(OBJDIR)/$*.o
# We use the same compiler whether it's c or cpp so I don't have to
# bother with extern (C) {...}.  We should only have cpp on the list
# when were compiling for Fungimol.
$(OBJDIR)/%.o: %.cpp
	$(CC) $(INCLUDES) $(C_FLAGS) -c $*.cpp -o $(OBJDIR)/$*.o

.PHONY: all

# We need libbren.a to build brenfort.  If we don't have
# libbren.a, then don't try to build brenfort.
ifeq ($(wildcard libbren.a),)
   ifeq ($(FUNGIMOL),)
      all: brennerfile TAGS
   else
      all: $(OBJDIR)/BrennerAtom.so
      $(OBJDIR)/BrennerAtom.so: $(OBJS) Makefile
	g++ $(OBJS) -fpic -shared -o $(OBJDIR)/BrennerAtom.so
   endif
else
   all: brenfort libbrenner.a brennerfile brennerrasmol TAGS
endif

ifneq ($(wildcard /usr/share/brennermd/inter*.d),)
   # If we can get the inter*.d files from the brennermd package, and we
   # don't already have them, get them.
   INTERSTARDOTD = inter2d_iv.d inter3d_ch.d inter3d_h.d inter3d_iv_new.d \
      inter3dtors.d
   $(INTERSTARDOTD): /usr/share/brennermd/inter*.d
	ln -s /usr/share/brennermd/inter*.d .
else
   INTERSTARDOTD :=
endif

$(OBJDIR)/CCNeighborState.o: CCNeighborState.c CCNeighborState.h

old.d: input.d coord_safe.d
	cp coord_safe.d coord.d
	-rm -f xmol.d output.d old.d
	time $(BRENNERMD)
	perl ./normalize_d.pl < xmol.d > old.d

.PHONY: test
new.d coord.d test: $(INTERSTARDOTD) brennerfile old.d
	cp coord_safe.d coord.d
	-rm -f xmol.d new.d
	@#If you change the number of timesteps (that is, the number
	@#immediately after the -s), be sure to change the
	@#first number in input.d to be the same, otherwise the movie
	@#lengths will be different and the comparison will fail.
	time ./brennerfile -s200 -t300
	perl ./normalize_d.pl < xmol.d > new.d
	@# We need the larger error bound if we're doing single precision.
	perl ./diff_d.pl old.d new.d 0.004

# Run the entire regression test, including compilations from scratch.
# This contrasts with test, which just checks the current build.
.PHONY: fulltest
fulltest:
	$(MAKE) clean
	$(MAKE) CC=/usr/local/bin/gcc
	$(MAKE) test
	$(MAKE) clean
	$(MAKE) CC=/usr/local/bin/g++
	$(MAKE) test
	$(MAKE) clean
	@# Next only verifies that it compiles.
	$(MAKE) CC=/usr/local/bin/g++ FUNGIMOL=yes
	$(MAKE) clean

# To pass the diff check at the end of timerun, you have to specify
# thermostat 4 and 10 steps between xmol writes in input.d.  That is,
# the first line of input.d can look like this:
# 20 1 4 10  / # steps, # steps between data, thermostat, xmol writes 

# For this test data, stopping the program with the debugger and
# looking at the stack 8 times gave 4 stops in the neighboring code,
# so we spend somewhere around half our time there.  Tim Freeman 29 Jul 2000.

EXE := brennerfile
.PHONY: timerunNoDiff
timerunNoDiff: $(INTERSTARDOTD)
	cp coord_safe.d coord.d
	-rm -f xmol.d
	time ./$(EXE) -s200 -t300

.PHONY: timerun
timerun:
	-rm -f new.d
	$(MAKE) timerunNoDiff
	perl ./normalize_d.pl < xmol.d > new.d
	perl ./diff_d.pl new.d correct_timerun.d 0.003

# The diff target generates code in both single precision mode and
# double precision mode, runs the two executables on the same inputs,
# and compares the output xmol.d files and the standard output for
# significant differences.  In the text saved from standard output, if
# there's the string "CMP", we expect a floating point number
# immediately afterward and we'll compare them for approximate
# equality in diff_run.pl.
.PHONY: diff
diff:
	-rm -f new.d xmol.d
	$(MAKE) brennerfile.single DEBUG=yes
	$(MAKE) timerunNoDiff EXE=brennerfile.single> out.single
	perl ./normalize_d.pl < xmol.d > new.d
	@#Can't use make clean here, since that would delete out.single.
	-rm -f old.d xmol.d
	$(MAKE) brennerfile.double DEBUG=yes
	$(MAKE) timerunNoDiff EXE=brennerfile.double> out.double
	perl ./normalize_d.pl < xmol.d > old.d
	perl ./diff_run.pl out.single out.double
	perl ./diff_d.pl old.d new.d 

.PHONY: clean
clean:
	-rm -f libbrenner.a *.o new.d old.d xmol.d output.d coord.d \
	   eigenenergies.d eigenvectors.d dos.d brennerfile ldos.d \
	   load.d max_ke.d overwrite.d pair_energy.d TAGS \
	   out.single out.double .'#'* inter*.d
	-rm -rf single double brennerfile.single brennerfile.double fungimol

brennerfile: $(OBJS) bren_file.c
	$(CC) -o brennerfile $(INCLUDES) $(C_FLAGS) bren_file.c $(OBJS) -lm

# We always have to rebuild brennerfile.single and ....double so
# declare them PHONY. 
.PHONY: brennerfile.single
brennerfile.single:
	-rm -f brennerfile
	-mkdir single
	$(MAKE) brennerfile SINGLE=yes OBJDIR=single 
	mv brennerfile brennerfile.single

.PHONY: brennerfile.double
brennerfile.double:
	-rm -f brennerfile
	-mkdir double
	$(MAKE) OBJDIR=double brennerfile
	mv brennerfile brennerfile.double

brennerrasmol: $(OBJS) bren_ras.c libbrenner.a
	$(CC) -o brennerrasmol $(INCLUDES) $(C_FLAGS) bren_ras.c -L. -lbrenner $(LIBS)

libbrenner.a : $(OBJS)
	rm -f $@
	ar cq $@ $(OBJS)
	ranlib $@

brenfort.o: brenfort.c brenfort.h brenner.h

MISC_OBJS = kt.o xmass.o myassert.o
brenfort: libbren.a brenfort.o brenffile.c $(MISC_OBJS)
	$(CC) -o brenfort $(C_FLAGS)  -Wno-missing-prototypes brenfort.o  $(MISC_OBJS) brenffile.c -L. -lbren -lf2c -lm

brenfras: libbren.a brenfort.o brenfras.c $(MISC_OBJS)
	$(CC) -o brenfras $(C_FLAGS) $(INCLUDES) -Wno-missing-prototypes brenfort.o $(MISC_OBJS) brenfras.c -L. $(LIBS) -lbren -lf2c -lm

tar:
	cd .. ; tar cf brennerc.tar `cat brennerc/TARGETS`; gzip -9 brennerc.tar
TAGS:
	etags *.c *.h
