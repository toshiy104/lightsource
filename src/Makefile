SRCS=light.cpp
OBJS=$(SRCS:.cpp=.o)
PROG=light

#CCOPTIONS= -O -DDEBUG
CCOPTIONS= -std=c++11 -Wall -O3

CPP = g++

.cpp.o:
	$(CPP) -c $(ALLCPPFLAGS) $(CCOPTIONS) $(CDEBUGFLAGS) $(INCLUDE) -o $@ -c $<

#-----------------------------------------------------------------------------
# crop
#-----------------------------------------------------------------------------
all:: $(PROG)

$(PROG): $(OBJS)
	$(CPP) -o $@ $(CDEBUGFLAGS) $(OBJS) -lm

clean::
	rm -f $(OBJS) $(PROG) $(PROG).exe
