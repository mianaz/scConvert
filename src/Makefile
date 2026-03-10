# Makefile for scConvert CLI
#
# Builds standalone scconvert binary for streaming h5ad ↔ h5seurat conversion.
# Requires HDF5 C library (brew install hdf5, or system package).

CC       ?= cc
CFLAGS   ?= -O2 -Wall -Wextra -Wpedantic -std=c11
LDFLAGS  ?=

# Auto-detect HDF5 via pkg-config, then brew, then fallback
HDF5_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null)
HDF5_LIBS   := $(shell pkg-config --libs hdf5 2>/dev/null)

ifeq ($(HDF5_LIBS),)
  # Try Homebrew
  BREW_PREFIX := $(shell brew --prefix hdf5 2>/dev/null)
  ifneq ($(BREW_PREFIX),)
    HDF5_CFLAGS := -I$(BREW_PREFIX)/include
    HDF5_LIBS   := -L$(BREW_PREFIX)/lib -lhdf5
  else
    # Fallback
    HDF5_CFLAGS := -I/usr/local/include
    HDF5_LIBS   := -L/usr/local/lib -lhdf5
  endif
endif

SRCS     = main.c sc_h5ad.c sc_h5mu.c sc_modality.c sc_sparse.c sc_dataframe.c sc_groups.c sc_util.c
OBJS     = $(SRCS:.c=.o)
TARGET   = scconvert

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(HDF5_LIBS)

%.o: %.c sc_convert.h
	$(CC) $(CFLAGS) $(HDF5_CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

install: $(TARGET)
	install -m 755 $(TARGET) /usr/local/bin/
