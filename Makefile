CFLAGS	      =  -g -DXWINDOWS -I/usr/X11/include

DEST	      = .

EXTHDRS	      =

HDRS	      = 

INSTALL	      = /etc/install

LD	      = CC

LDFLAGS	      =

LDLIBS	      = /usr/X11/lib

MAKEFILE      = Makefile

OBJS	      = main.o \
		lgar_routines.o \
		linked_list.o \
		mem_funcs.o \
		soil_funcs.o \
		util_funcs.o \
                xcategorize.o \
		xinitialize.o \
		xplot.o \
		xshow.o

PRINT	      = pr

PROGRAM       = xlgar 

SHELL	      = /bin/sh

SRCS	      = main.c \
		lgar_routines.c \
		linked_list.c \
		mem_funcs.c \
		soil_funcs.c \
		util_funcs.c \
		xcategorize.c \
		xinitialize.c \
		xplot.c \
		xshow.c


all:		$(PROGRAM)

$(PROGRAM):     $(OBJS) $(LDLIBS)
		@echo "Linking $(PROGRAM)"
		@$(LD) $(LDFLAGS) $(OBJS) -lc -lm -L/usr/X11R6/lib -lXt \
		-L/usr/X11R6/lib -lX11 -L/$(LDLIBS) -o $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS) core

clobber:;	@rm -f $(OBJS) $(PROGRAM) core tags

depend:;	@mkmf -f $(MAKEFILE) ROOT=$(ROOT)

echo:;		@echo $(HDRS) $(SRCS)

index:;		@ctags -wx $(HDRS) $(SRCS)

