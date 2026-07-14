.PHONY: all clean

OUTPUT=pde-core
VISUALIZER=pde-visualizer
SHARED_LIB=libpdecore.so

CXX = clang++ -std=c++20

OPTIMIZE = -O3 -march=native -flto -fno-fast-math

OMP_PREFIX = /opt/homebrew/opt/libomp
OMP_INC = -I$(OMP_PREFIX)/include
OMP_LIB = -L$(OMP_PREFIX)/lib -lomp
OMP_FLAGS = -Xpreprocessor -fopenmp $(OMP_INC) $(OMP_LIB)

GLUI_PREFIX = /opt/homebrew/opt/glui
GLUI_INC = -I$(GLUI_PREFIX)/include/GL
GLUI_LIB = -L$(GLUI_PREFIX)/lib -lglui

AARCH = $(shell uname -m)
ifeq ($(AARCH), arm64)
	ARCH = -arch arm64 -mcpu=apple-m1 
else
	ARCH = -arch x86_64
endif	

all: core-debug

core-debug:
	$(CXX) -fsanitize=address $(ARCH) \
		$(GLUI_INC) \
		$(GLUI_LIB) \
		-framework GLUT \
		-framework OpenGL \
		-DGL_SILENCE_DEPRECATION \
		-O0 \
		-g main.cpp -o ${OUTPUT}

core-release:
	clang++ $(ARCH) \
		$(GLUI_INC) \
		$(GLUI_LIB) \
		-framework GLUT \
		-framework OpenGL \
		-DGL_SILENCE_DEPRECATION \
		$(OPTIMIZE) \
		main.cpp -o ${OUTPUT}
