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
	$(EXTERN_COMMON_PATH)/get_file_size.cpp \
	$(EXTERN_PROB_PATH)/problem_6.cpp \
	$(EXTERN_PROB_PATH)/cl_util.cpp \
	$(EXTERN_PROB_PATH)/cl_wrap.cpp \
	Content/param.cpp \
	Renderer/MetalRenderer.m \
	Application/AppDelegate.m \
	Application/main.mm

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
	-fvisibility=hidden \
	-fvisibility-inlines-hidden  \
	-fstrict-aliasing \
	-fno-exceptions \
	-fno-rtti \
	-O3 \
	-ffast-math \
	-fno-unsafe-math-optimizations \
	-fhonor-infinities \
	-fomit-frame-pointer \
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
	-Wno-missing-declarations \
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
	-Wc++11-extensions \
	-Winvalid-offsetof \
	-mmacosx-version-min=10.13 \
	-I$(EXTERN_COMMON_PATH)/.. \
	-I$(EXTERN_COMMON_PATH) \
	-I$(EXTERN_PROB_PATH) \
	-I./Content \
	-I./Renderer \
	-DMINIMAL_TREE=1 \
	-DCLANG_QUIRK_0001=1 \
	-DOCL_QUIRK_0001=1 \
	-DOCL_QUIRK_0005=1 \
	-DOCL_KERNEL_BUILD_VERBOSE=0 \
	-DOCL_BUFFER_COPY=0 \
	-DFRAME_RATE=0 \
	-DNDEBUG

ifeq ($(UNAME), Darwin)

	LINKFLAGS += -framework OpenCL \
				 -framework Metal \
				 -framework MetalKit \
				 -framework CoreServices \
				 -framework CoreVideo \
				 -framework AppKit

	LINKFLAGS += \
		-mmacosx-version-min=10.13 \
		-fobjc-arc \
		-fobjc-link-runtime \
		-stdlib=libc++ \
		-arch $(HOSTTYPE)

	CFLAGS += -arch $(HOSTTYPE)

	ifeq ($(HOSTTYPE), arm64)
		CFLAGS += -march=armv8.4-a -mtune=native

	else ifeq ($(HOSTTYPE), x86_64)
		CFLAGS += -march=native -mtune=native

	endif
else
	# Nothing to do at this stage.
endif

CXX = $(CC)
CXXFLAGS = $(CFLAGS)

BINARY := $(TARGET)

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

