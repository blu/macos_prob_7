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

CC = /usr/bin/clang
LD = /usr/bin/clang++
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
	-fstrict-aliasing \
	-fno-exceptions \
	-fno-rtti \
	-Ofast \
	-flto \
	-funroll-loops \
	-Wtrigraphs \
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
	-Wreturn-type \
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
	-Winvalid-offsetof \
	-isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.12.sdk \
	-mmacosx-version-min=10.10 \
	-I$(EXTERN_COMMON_PATH)/.. \
	-I$(EXTERN_COMMON_PATH) \
	-I$(EXTERN_PROB_PATH) \
	-I./GLEssentials/Source \
	-I./GLEssentials/Source/Classes \
	-I./GLEssentials/Source/Classes/OSX \
	-DMINIMAL_TREE=1 \
	-DCLANG_QUIRK_0001=1 \
	-DCLANG_QUIRK_0002=1 \
	-DOCL_QUIRK_0001=1 \
	-DOCL_KERNEL_BUILD_VERBOSE=0 \
	-DNDEBUG

ifeq ($(UNAME), Darwin)

	LINKFLAGS += -framework OpenCL -framework OpenGL -framework CoreServices

	LINKFLAGS += \
		-mmacosx-version-min=10.10 \
		-fobjc-arc \
		-fobjc-link-runtime \
		-stdlib=libc++

	CFLAGS += -arch $(HOSTTYPE) -march=native -mtune=native
	LINKFLAGS += -arch $(HOSTTYPE)

else
endif

CXX = $(CC)
CXXFLAGS = $(CFLAGS)

BINARY := $(TARGET).app/Contents/MacOS/$(TARGET)

all: $(BINARY)

%.o : %.cpp
	$(CC) -x c++ -std=c++11 -stdlib=libc++ $(CXXFLAGS) -c $< -o $@

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

