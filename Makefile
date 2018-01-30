#
# problem_* build script for macOS
#

# external control variables:
#
# HOSTTYPE                              -- same as shell's HOSTTYPE
# EXTERN_COMMON_PATH                    -- external dependence: common path
# EXTERN_PROB_PATH                      -- external dependence: problem path

TARGET = problem_7
UNAME := $(shell uname)

SRCS = \
	$(EXTERN_COMMON_PATH)/pthread_barrier.cpp \
	$(EXTERN_COMMON_PATH)/prim_mono_view.cpp \
	$(EXTERN_COMMON_PATH)/get_file_size.cpp \
	$(EXTERN_COMMON_PATH)/util_gl.cpp \
	$(EXTERN_PROB_PATH)/problem_6.cpp \
	$(EXTERN_PROB_PATH)/cl_util.cpp \
	$(EXTERN_PROB_PATH)/cl_wrap.cpp \
	GLEssentials/Source/param.cpp \
	GLEssentials/Source/Classes/OSX/AppDelegate.m \
	GLEssentials/Source/Classes/OSX/GLEssentialsWindowController.m \
	GLEssentials/Source/Classes/OSX/GLEssentialsGLView.m \
	GLEssentials/Source/Classes/OSX/GLEssentialsFullscreenWindow.m \
	GLEssentials/Source/Classes/OpenGLRenderer.m \
	GLEssentials/Source/main.mm
OBJS0 = $(SRCS:.cpp=.o)
OBJS1 = $(OBJS0:.mm=.o)
OBJS = $(OBJS1:.m=.o)

# for now stick with /usr/bin/clang, as this recent gcc fails to build us:
#
# gcc (Ubuntu/Linaro 4.6.1-9ubuntu3) 4.6.1

CC = /usr/bin/clang
LD = /usr/bin/clang++
# CFLAGS += -pipe -fmessage-length=0 -fno-rtti -O3 -funroll-loops -ffast-math -fstrict-aliasing -Wtrigraphs -Wreturn-type -Wunused-variable -Wunused-value -DGLEW_NO_GLU
CFLAGS += \
	-Wno-logical-op-parentheses \
	-Wno-bitwise-op-parentheses \
	-Wno-parentheses
CFLAGS += \
	-fmessage-length=0 \
	-fdiagnostics-show-note-include-stack \
	-fmacro-backtrace-limit=0 \
	-fmodules \
	-gmodules \
	-fmodules-prune-interval=86400 \
	-fmodules-prune-after=345600 \
	-Wnon-modular-include-in-framework-module \
	-Werror=non-modular-include-in-framework-module \
	-fvisibility=hidden \
	-fvisibility-inlines-hidden  \
	-fasm-blocks \
	-fstrict-aliasing \
	-fno-exceptions \
	-fno-rtti \
	-Ofast \
	-Wno-missing-field-initializers \
	-Wno-missing-prototypes \
	-Werror=return-type \
	-Wunreachable-code \
	-Werror=deprecated-objc-isa-usage \
	-Werror=objc-root-class \
	-Wno-non-virtual-dtor \
	-Wno-overloaded-virtual \
	-Wno-exit-time-destructors \
	-Wno-missing-braces \
	-Wunused-function \
	-Wno-unused-label \
	-Wno-unused-parameter \
	-Wunused-variable \
	-Wunused-value \
	-Wempty-body \
	-Wuninitialized \
	-Wconditional-uninitialized \
	-Wno-shadow \
	-Wno-four-char-constants \
	-Wno-conversion \
	-Wconstant-conversion \
	-Wint-conversion \
	-Wbool-conversion \
	-Wenum-conversion \
	-Wno-newline-eof \
	-Wno-c++11-extensions \
	-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk \
	-Winvalid-offsetof \
	-mmacosx-version-min=10.10 \
	-I$(EXTERN_COMMON_PATH) \
	-I$(EXTERN_PROB_PATH) \
	-I./GLEssentials/Source \
	-I./GLEssentials/Source/Classes \
	-I./GLEssentials/Source/Classes/OSX \
	-DMINIMAL_TREE=1 \
	-DFB_RES_FIXED_W=512 \
	-DFB_RES_FIXED_H=512 \
	-DRAY_HIGH_PRECISION_RCP_DIR=1 \
	-DWORKFORCE_NUM_THREADS=2 \
	-DDIVISION_OF_LABOR_VER=2 \
	-DBOUNCE_COMPUTE_VER=1 \
	-DAO_NUM_RAYS=16 \
	-DCLANG_QUIRK_0001=1 \
	-DCLANG_QUIRK_0002=1 \
	-DNDEBUG

ifeq ($(UNAME), Darwin)

	ifeq ($(HOSTTYPE), powerpc)
		CFLAGS += -arch ppc -mtune=7450 -faltivec
		LINKFLAGS += -framework OpenGL -framework CoreServices -arch ppc -Wl,-Y,1455
	else
		LINKFLAGS += \
			-mmacosx-version-min=10.10 \
			-Xlinker -object_path_lto \
			-Xlinker build/GLEssentials.build/Release/GLEssentials-OSX.build/Objects-normal/x86_64/GLEssentials_lto.o \
			-fobjc-arc \
			-fobjc-link-runtime \
			-stdlib=libc++

		ifeq ($(HOSTTYPE), x86_64)
			CFLAGS += -arch x86_64 -march=native -mtune=native
			LINKFLAGS += -framework OpenGL -framework CoreServices -arch x86_64
		else
			CFLAGS += -arch i386 -march=native -mtune=native
			LINKFLAGS += -framework OpenGL -framework CoreServices -arch i386
		endif
	endif

else
endif

CXX = $(CC)
CXXFLAGS = $(CFLAGS)

BINARY := $(TARGET).app/Contents/MacOS/$(TARGET)

all: $(BINARY)

%.o : %.cpp
	$(CC) -x c++ -std=gnu++11 -stdlib=libc++ $(CXXFLAGS) -c $< -o $@

%.o : %.mm
	$(CC) -x objective-c++ -std=gnu++11 -fobjc-arc -stdlib=libc++ $(CXXFLAGS) -c $< -o $@

%.o : %.m
	$(CC) -x objective-c -std=gnu99 -fobjc-arc $(CFLAGS) -c $< -o $@

$(BINARY) : $(OBJS)
	$(LD) -o $@ $(OBJS) $(LINKFLAGS)

clean:
	$(RM) $(BINARY) $(OBJS)

clobber: clean
	$(RM) *.bak *~

.PHONY: all clean clobber

