TOPSRCDIR := .
include $(TOPSRCDIR)/make.inc.gnu  # make.inc.intel, make.inc.llvm

TARGET := drsyl

# Set global defines
DEFINES += -DNDEBUG -DALIGNMENT=64 #-DINTSCALING

# Select all C source files
SRCS := $(wildcard *.c)
OBJS := $(SRCS:.c=.o)

.SUFFIXES: .c .o

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

%.o : %.c
	$(CC) $(CFLAGS) $(DEFINES) -c $< -o $@

clean: 
	rm -f $(TARGET) *.o ipo_out.optrpt


.PHONY: clean
