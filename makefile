INCFLAGS=-I./lib/ -I.
CFLAGS=-O3 $(INCFLAGS)

OBJDIR = bin
SRCDIR = src

SRCS    = $(shell find $(SRCDIR) -name '*.c')
OBJS    = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(SRCS))
CC		= gcc
LIBS	= -lm -lGL -lGLU -lglut

MKDIR_P = mkdir -p

all: molas

molas: directories $(OBJS)
	$(CC) $(OBJS) $(LIBS) $(CFLAGS) -o $(OBJDIR)/molas

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(OPTS) -O3 -Wall -c $< -o $@
	
directories:
	$(MKDIR_P) $(OBJDIR)

clean:
	rm -f $(OBJS)
	rm -rf $(OBJDIR) 

depend:
	makedepend $(INCFLAGS) $(SRCS)


# DO NOT DELETE THIS LINE