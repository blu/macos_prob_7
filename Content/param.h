#ifndef param_H__
#define param_H__

enum KernelParamType { // type of kernel parameters
	PARAM_TYPE_BUFFER = 0,
	PARAM_TYPE_IMAGE
};

struct cli_param {
	enum Flags { // flags comprising cli_param::flags
		BIT_REPORT_CAPS              = 1,
		BIT_DISCARD_PLATFORM_VERSION = 2,
		BIT_DISCARD_DEVICE_VERSION   = 4,
		BIT_REPORT_KERNEL_TIME       = 8
	};

	unsigned platform_idx;  // cl platform index
	unsigned device_idx;    // cl device index
	unsigned flags;         // custom control flags
	enum KernelParamType kern_param_type;

	unsigned image_w;       // frame width
	unsigned image_h;       // frame height
	unsigned image_hz;      // frame rate target Hz
	unsigned frames;        // frames to run
	unsigned workgroup_size;// cl workgroup size
};


#ifdef __cplusplus
extern "C" {
#endif

extern struct cli_param param;

int parseCLI(int, const char **, struct cli_param *);
int content_init(void);
int content_deinit(void);
int content_frame(void *target);

#ifdef __cplusplus
}
#endif

#endif // param_H__
