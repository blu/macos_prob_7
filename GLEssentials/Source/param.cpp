#include <CL/cl.h>
#if CL_VERSION_1_2 == 0
	#error required CL_VERSION_1_2 to compile
#endif
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "param.h"
#include "timer.h"
#include "native_gl.h"
#include "vectnative.hpp"
#include "pure_macro.hpp"
#include "cl_util.hpp"
#include "cl_wrap.hpp"
#include "platform.hpp"
#include "prim_mono_view.hpp"
#include "scoped.hpp"
#include "stream.hpp"
#include "array.hpp"
#include "problem_6.hpp"

// verify iostream-free status
#if _GLIBCXX_IOSTREAM
#error rogue iostream acquired
#endif

// verify tree minimalism
#if MINIMAL_TREE == 0
#error opencl kernels used expect a minimal tree
#endif

namespace stream {

// deferred initialization by main()
in cin;
out cout;
out cerr;

} // namespace stream

static const char arg_prefix[]                   = "-";
static const char arg_report_caps[]              = "report_caps";
static const char arg_discard_platform_version[] = "discard_platform_version";
static const char arg_discard_device_version[]   = "discard_device_version";
static const char arg_platform[]                 = "platform";
static const char arg_device[]                   = "device";
static const char arg_use_images[]               = "use_images";
static const char arg_report_kernel_time[]       = "report_kernel_time";
static const char arg_screen[]                   = "screen";
static const char arg_bitness[]                  = "bitness";
static const char arg_fsaa[]                     = "fsaa";
static const char arg_frames[]                   = "frames";

static const size_t n_buffering = 2;

namespace testbed {

template <>
class scoped_functor< cl_context > {
public:
	void operator()(cl_context* arg) {
		assert(0 != arg);
		if (cl_context(0) != *arg)
			clReleaseContext(*arg);
	}
};

template <>
class scoped_functor< cl_command_queue > {
public:
	void operator()(cl_command_queue* arg) {
		assert(0 != arg);
		if (cl_command_queue(0) != *arg)
			clReleaseCommandQueue(*arg);
	}
};

template <>
class scoped_functor< cl_mem > {
public:
	void operator()(cl_mem* arg) {
		assert(0 != arg);
		for (size_t i = 0; i < n_buffering; ++i)
			if (cl_mem(0) != arg[i])
				clReleaseMemObject(arg[i]);
	}
};

template <>
class scoped_functor< cl_program > {
public:
	void operator()(cl_program* arg) {
		assert(0 != arg);
		if (cl_program(0) != *arg)
			clReleaseProgram(*arg);
	}
};

template <>
class scoped_functor< cl_kernel > {
public:
	void operator ()(cl_kernel* arg) {
		assert(0 != arg);
		if (cl_kernel(0) != *arg)
			clReleaseKernel(*arg);
	}
};

template <>
class scoped_functor< GLuint > {
public:
	void operator ()(GLuint* arg) {
		assert(0 != arg);
		glDeleteTextures(n_buffering, arg);
	}
};

template <>
class scoped_functor< FILE > {
public:
	void operator()(FILE* arg) {
		assert(0 != arg);
		fclose(arg);
	}
};

template < typename T >
class generic_free {
public:
	void operator()(T* arg) {
		assert(0 != arg);
		std::free(arg);
	}
};

} // namespace testbed


static bool
validate_fullscreen(
	const char* const string,
	unsigned& screen_w,
	unsigned& screen_h)
{
	if (0 == string)
		return false;

	unsigned x, y, hz;

	if (3 != sscanf(string, "%u %u %u", &x, &y, &hz))
		return false;

	if (!x || !y || !hz)
		return false;

	screen_w = x;
	screen_h = y;

	return true;
}


static bool
validate_bitness(
	const char* const string,
	unsigned (& screen_bitness)[4])
{
	if (0 == string)
		return false;

	unsigned bitness[4];

	if (4 != sscanf(string, "%u %u %u %u",
			&bitness[0],
			&bitness[1],
			&bitness[2],
			&bitness[3]))
	{
		return false;
	}

	if (!bitness[0] || 16 < bitness[0] ||
		!bitness[1] || 16 < bitness[1] ||
		!bitness[2] || 16 < bitness[2] ||
		16 < bitness[3])
	{
		return false;
	}

	screen_bitness[0] = bitness[0];
	screen_bitness[1] = bitness[1];
	screen_bitness[2] = bitness[2];
	screen_bitness[3] = bitness[3];

	return true;
}

int
parseCLI(
	const int argc,
	char** const argv,
	cli_param* param) {

	const size_t prefix_len = std::strlen(arg_prefix);
	bool success = true;

	for (int i = 1; i < argc && success; ++i) {
		if (std::strncmp(argv[i], arg_prefix, prefix_len)) {
			success = false;
			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_report_caps)) {
			param->flags |= cli_param::BIT_REPORT_CAPS;
			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_discard_platform_version)) {
			param->flags |= cli_param::BIT_DISCARD_PLATFORM_VERSION;
			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_discard_device_version)) {
			param->flags |= cli_param::BIT_DISCARD_DEVICE_VERSION;
			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_platform)) {
			if (++i == argc || 1 != sscanf(argv[i], "%u", &param->platform_idx))
				success = false;

			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_device)) {
			if (++i == argc || 1 != sscanf(argv[i], "%u", &param->device_idx))
				success = false;

			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_use_images) && PARAM_TYPE_BUFFER == param->kern_param_type) {
			param->kern_param_type = PARAM_TYPE_IMAGE;
			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_report_kernel_time)) {
			param->flags |= cli_param::BIT_REPORT_KERNEL_TIME;
			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_screen)) {
			if (++i == argc || !validate_fullscreen(argv[i], param->image_w, param->image_h))
				success = false;

			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_bitness)) {
			if (++i == argc || !validate_bitness(argv[i], param->bitness))
				success = false;

			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_fsaa)) {
			if (++i == argc || 1 != sscanf(argv[i], "%u", &param->fsaa))
				success = false;

			continue;
		}

		if (!std::strcmp(argv[i] + prefix_len, arg_frames)) {
			if (++i == argc || 1 != sscanf(argv[i], "%u", &param->frames))
				success = false;

			continue;
		}

		success = false;
	}

	if (!success) {
		stream::cerr << "usage: " << argv[0] << " [<option> ...]\n"
			"options (multiple args to an option must constitute a single string, eg. -foo \"a b c\"):\n"
			"\t" << arg_prefix << arg_report_caps << "\t\t\t: report CL capabilities\n"
			"\t" << arg_prefix << arg_discard_platform_version << "\t: discard advertised platform version when producing platform report\n"
			"\t" << arg_prefix << arg_discard_device_version << "\t\t: discard advertised device version when producing device report\n"
			"\t" << arg_prefix << arg_platform << " <index>\t\t: use platform of specified index\n"
			"\t" << arg_prefix << arg_device << " <index>\t\t\t: use device of specified index\n"
			"\t" << arg_prefix << arg_use_images << "\t\t\t: use images instead of buffers as source arguments\n"
			"\t" << arg_prefix << arg_report_kernel_time << "\t\t: report CL kernel time\n"
			"\t" << arg_prefix << arg_screen << " <width> <height> <Hz>\t: set fullscreen output of specified geometry and refresh\n"
			"\t" << arg_prefix << arg_bitness << " <r> <g> <b> <a>\t: set GLX config of specified RGBA bitness; default is screen's bitness\n"
			"\t" << arg_prefix << arg_fsaa << " <positive_integer>\t: set GL fullscreen antialiasing; default is none\n"
			"\t" << arg_prefix << arg_frames << " <unsigned_integer>\t: set number of frames to run; default is max unsigned int\n";

		return 1;
	}

	return 0;
}

static matx4 transpose(const matx4& src) {
	simd::f32x4 r0, r1, r2, r3;
	simd::transpose4x4(src[0], src[1], src[2], src[3], r0, r1, r2, r3);
	return matx4(r0, r1, r2, r3);
}

class matx3_rotate : public matx3 {
	matx3_rotate();

public:
	matx3_rotate(
		const float a,
		const float x,
		const float y,
		const float z) {

		simd::f32x4 sin_ang;
		simd::f32x4 cos_ang;
		simd::sincos(simd::f32x4(a, simd::flag_zero()), sin_ang, cos_ang);

		const float sin_a = sin_ang[0];
		const float cos_a = cos_ang[0];

		m[0] = simd::f32x4(x * x + cos_a * (1 - x * x),         x * y - cos_a * (x * y) + sin_a * z, x * z - cos_a * (x * z) - sin_a * y, simd::flag_zero());
		m[1] = simd::f32x4(y * x - cos_a * (y * x) - sin_a * z, y * y + cos_a * (1 - y * y),         y * z - cos_a * (y * z) + sin_a * x, simd::flag_zero());
		m[2] = simd::f32x4(z * x - cos_a * (z * x) + sin_a * y, z * y - cos_a * (z * y) - sin_a * x, z * z + cos_a * (1 - z * z),         simd::flag_zero());
	}
};


class matx4_rotate : public matx4 {
	matx4_rotate();

public:
	matx4_rotate(
		const float a,
		const float x,
		const float y,
		const float z) {

		simd::f32x4 sin_ang;
		simd::f32x4 cos_ang;
		simd::sincos(simd::f32x4(a, simd::flag_zero()), sin_ang, cos_ang);

		const float sin_a = sin_ang[0];
		const float cos_a = cos_ang[0];

		m[0] = simd::f32x4(x * x + cos_a * (1 - x * x),         x * y - cos_a * (x * y) + sin_a * z, x * z - cos_a * (x * z) - sin_a * y, 0.f);
		m[1] = simd::f32x4(y * x - cos_a * (y * x) - sin_a * z, y * y + cos_a * (1 - y * y),         y * z - cos_a * (y * z) + sin_a * x, 0.f);
		m[2] = simd::f32x4(z * x - cos_a * (z * x) + sin_a * y, z * y - cos_a * (z * y) - sin_a * x, z * z + cos_a * (1 - z * z),         0.f);
		m[3] = simd::f32x4(0.f,                                 0.f,                                 0.f,                                 1.f);
	}
};


static inline float
wrap_at_period(
	const float x,
	const float period) {

	const simd::f32x4 vx = simd::f32x4(x, simd::flag_zero());
	const simd::f32x4 vperiod = simd::f32x4(period, simd::flag_zero());
	const simd::u32x4 mask = vx >= vperiod;
	return x - simd::mask(vperiod, mask)[0];
}


static inline int32_t
reset_at_period(
	const int32_t x,
	const int32_t period) {

	const simd::s32x4 vx = simd::s32x4(x, simd::flag_zero());
	const simd::s32x4 vperiod = simd::s32x4(period, simd::flag_zero());
	const simd::u32x4 mask = vx < vperiod;
	return simd::mask(vx, mask)[0];
}

////////////////////////////////////////////////////////////////////////////////
// scene support
////////////////////////////////////////////////////////////////////////////////

class Scene {
protected:
	// scene offset in model space
	float offset_x;
	float offset_y;
	float offset_z;

	// scene orientation
	float azim; // around z
	float decl; // around x
	float roll; // around y

	// scene camera position
	float cam_x;
	float cam_y;
	float cam_z;

public:
	Scene()
	: offset_x(0.f)
	, offset_y(0.f)
	, offset_z(0.f)
	, azim(0.f)
	, decl(0.f)
	, roll(0.f)
	, cam_x(0.f)
	, cam_y(0.f)
	, cam_z(0.f) {
	}

	virtual bool init(Timeslice& scene) = 0;
	virtual bool frame(Timeslice& scene, const float dt) = 0;

	// scene offset in model space
	float get_offset_x() const {
		return offset_x;
	}

	float get_offset_y() const {
		return offset_y;
	}

	float get_offset_z() const {
		return offset_z;
	}

	// scene orientation
	float get_azim() const {
		return azim;
	}

	float get_decl() const {
		return decl;
	}

	float get_roll() const {
		return roll;
	}

	// scene camera position
	float get_cam_x() const {
		return cam_x;
	}

	float get_cam_y() const {
		return cam_y;
	}

	float get_cam_z() const {
		return cam_z;
	}
};

////////////////////////////////////////////////////////////////////////////////
// Scene1: Deathstar Treadmill
////////////////////////////////////////////////////////////////////////////////

class Scene1 : public virtual Scene {

	// scene camera properties
	float accum_x;
	float accum_y;

	enum {
		grid_rows = 20,
		grid_cols = 20,
		dist_unit = 1
	};

	float accum_time;
	float generation;

	Array< Voxel > content;
	BBox contentBox;

	bool update(
		Timeslice& scene,
		const float generation);

	void camera(
		const float dt);

public:
	Scene1() : contentBox(BBox::flag_noinit()) {}

	// virtual from Scene
	bool init(
		Timeslice& scene);

	// virtual from Scene
	bool frame(
		Timeslice& scene,
		const float dt);
};


bool Scene1::init(
	Timeslice& scene) {

	accum_x = 0.f;
	accum_y = 0.f;
	accum_time = 0.f;
	generation = grid_rows;

	if (!content.setCapacity(grid_rows * grid_cols))
		return false;

	const float unit = dist_unit;
	const float alt = unit * .5f;
	contentBox = BBox();

	for (int y = 0; y < grid_rows; ++y)
		for (int x = 0; x < grid_cols; ++x) {
			const BBox box(
				vect3(x * unit,        y * unit,        0.f),
				vect3(x * unit + unit, y * unit + unit, alt * (rand() % 4 + 1)),
				BBox::flag_direct());

			contentBox.grow(box);
			content.addElement(Voxel(box.get_min(), box.get_max()));
		}

	return scene.set_payload_array(content, contentBox);
}


inline bool Scene1::update(
	Timeslice& scene,
	const float generation) {

	const float unit = dist_unit;
	const float alt = unit * .5f;
	size_t index = 0;
	contentBox = BBox();

	for (index = 0; index < (grid_rows - 1) * grid_cols; ++index) {
		contentBox.grow(content.getElement(index + grid_cols).get_bbox());
		content.getMutable(index) = content.getElement(index + grid_cols);
	}

	const float y = generation;

	for (int x = 0; x < grid_cols; ++x, ++index) {
		const BBox box(
			vect3(x * unit,        y * unit,        0.f),
			vect3(x * unit + unit, y * unit + unit, alt * (rand() % 4 + 1)),
			BBox::flag_direct());

		contentBox.grow(box);
		content.getMutable(index) = Voxel(box.get_min(), box.get_max());
	}

	return scene.set_payload_array(content, contentBox);
}


inline void Scene1::camera(
	const float dt) {

	const float period_x = 3.f; // seconds
	const float period_y = 2.f; // seconds

	const float deviate_x = 1.f / 32.f; // distance
	const float deviate_y = 1.f / 32.f; // distance

	const float roll_factor = 1 / 32.f; // of pi

	accum_x = wrap_at_period(accum_x + dt, period_x);
	accum_y = wrap_at_period(accum_y + dt, period_y);

	simd::f32x4 sin_xy;
	simd::f32x4 cos_x;
	simd::sincos(simd::f32x4(
		accum_x / period_x * float(M_PI * 2.0),
		accum_y / period_y * float(M_PI * 2.0), 0.f, 0.f),
		sin_xy, cos_x);

	cam_x = sin_xy[0] * deviate_x;
	cam_y = sin_xy[1] * deviate_y;
	roll  = cos_x[0] * float(M_PI * roll_factor);
}


bool Scene1::frame(
	Timeslice& scene,
	const float dt) {

	camera(dt);

	const float update_period = .25f;

	offset_y -= dist_unit * dt / update_period;
	accum_time += dt;

	if (accum_time < update_period)
		return scene.set_payload_array(content, contentBox);

	accum_time -= update_period;

	return update(scene, generation++);
}

////////////////////////////////////////////////////////////////////////////////
// Scene2: Sine Floater
////////////////////////////////////////////////////////////////////////////////

class Scene2 : virtual public Scene {

	// scene camera properties
	float accum_x;
	float accum_y;

	enum {
		grid_rows = 20,
		grid_cols = 20,
		dist_unit = 1
	};

	float accum_time;

	Array< Voxel > content;

	bool update(
		Timeslice& scene,
		const float dt);

	void camera(
		const float dt);

public:
	// virtual from Scene
	bool init(
		Timeslice& scene);

	// virtual from Scene
	bool frame(
		Timeslice& scene,
		const float dt);
};


bool Scene2::init(
	Timeslice& scene) {

	accum_x = 0.f;
	accum_y = 0.f;
	accum_time = 0.f;

	if (!content.setCapacity(grid_rows * grid_cols))
		return false;

	const float unit = dist_unit;
	BBox contentBox;

	for (int y = 0; y < grid_rows; ++y)
		for (int x = 0; x < grid_cols; ++x) {
			const BBox box(
				vect3(x * unit,        y * unit,        0.f),
				vect3(x * unit + unit, y * unit + unit, 1.f),
				BBox::flag_direct());

			contentBox.grow(box);
			content.addElement(Voxel(box.get_min(), box.get_max()));

		}

	return scene.set_payload_array(content, contentBox);
}


inline bool Scene2::update(
	Timeslice& scene,
	const float dt) {

	const float period = 2.f; // seconds

	accum_time = wrap_at_period(accum_time + dt, period);

	const float time_factor = simd::sin(simd::f32x4(accum_time / period * float(M_PI * 2.0), simd::flag_zero()))[0];

	const float unit = dist_unit;
	const float alt = unit * .5f;
	size_t index = 0;
	BBox contentBox;

	for (int y = 0; y < grid_rows; ++y)
		for (int x = 0; x < grid_cols; ++x, ++index) {
			const simd::f32x4 sin_xy = simd::sin(simd::f32x4(x * alt, y * alt, 0.f, 0.f));
			const BBox box(
				vect3(x * unit,        y * unit,        0.f),
				vect3(x * unit + unit, y * unit + unit, 1.f + time_factor * unit * (sin_xy[0] * sin_xy[1])),
				BBox::flag_direct());

			contentBox.grow(box);
			content.getMutable(index) = Voxel(box.get_min(), box.get_max());
		}

	return scene.set_payload_array(content, contentBox);
}


inline void Scene2::camera(
	const float dt) {

	const float period_x = 3.f; // seconds
	const float period_y = 2.f; // seconds

	const float deviate_x = 1.f / 32.f; // distance
	const float deviate_y = 1.f / 32.f; // distance

	const float roll_factor = 1 / 32.f; // of pi

	accum_x = wrap_at_period(accum_x + dt, period_x);
	accum_y = wrap_at_period(accum_y + dt, period_y);

	simd::f32x4 sin_xy;
	simd::f32x4 cos_x;
	simd::sincos(simd::f32x4(
		accum_x / period_x * float(M_PI * 2.0),
		accum_y / period_y * float(M_PI * 2.0), 0.f, 0.f),
		sin_xy, cos_x);

	cam_x = sin_xy[0] * deviate_x;
	cam_y = sin_xy[1] * deviate_y;
	roll  = cos_x[0] * float(M_PI * roll_factor);
}


bool Scene2::frame(
	Timeslice& scene,
	const float dt) {

	camera(dt);

	return update(scene, dt);
}

////////////////////////////////////////////////////////////////////////////////
// Scene3: Serpents
////////////////////////////////////////////////////////////////////////////////

class Scene3 : virtual public Scene {

	// scene camera properties
	float accum_x;
	float accum_y;

	enum {
		queue_length = 96,
		main_radius = 8,
	};

	float accum_time;

	Array< Voxel > content;

	bool update(
		Timeslice& scene,
		const float dt);

public:
	// virtual from Scene
	bool init(
		Timeslice& scene);

	// virtual from Scene
	bool frame(
		Timeslice& scene,
		const float dt);
};


bool Scene3::init(
	Timeslice& scene) {

	offset_x = 10.f; // matching the center of scene_1
	offset_y = 10.f; // matching the center of scene_1

	accum_x = 0.f;
	accum_y = 0.f;
	accum_time = 0.f;

	if (!content.setCapacity(queue_length * 2 + 1))
		return false;

	const float radius = main_radius;
	BBox contentBox;

	for (int i = 0; i < queue_length; ++i) {
		const float angle = float(M_PI * 2.0) * (i / float(queue_length));
		const float sin_i = simd::sin(simd::f32x4(i / float(queue_length) * float(M_PI * 16.0), simd::flag_zero()))[0];

		{
			const matx3_rotate rot(angle, 0.f, 0.f, 1.f);
			const vect3 pos = vect3(radius + sin_i, 0.f, sin_i) * rot;
			const BBox box(
				pos - vect3(.5f),
				pos + vect3(.5f),
				BBox::flag_direct());

			contentBox.grow(box);
			content.addElement(Voxel(box.get_min(), box.get_max()));
		}
		{
			const matx3_rotate rot(-angle, 0.f, 0.f, 1.f);
			const vect3 pos = vect3(radius * .5f + sin_i, 0.f, sin_i) * rot;
			const BBox box(
				pos - vect3(.25f),
				pos + vect3(.25f),
				BBox::flag_direct());

			contentBox.grow(box);
			content.addElement(Voxel(box.get_min(), box.get_max()));
		}
	}

	const BBox box(
		vect3(-main_radius, -main_radius, -.25f),
		vect3(+main_radius, +main_radius, +.25f),
		BBox::flag_direct());

	contentBox.grow(box);
	content.addElement(Voxel(box.get_min(), box.get_max()));

	return scene.set_payload_array(content, contentBox);
}


inline bool Scene3::update(
	Timeslice& scene,
	const float dt) {

	const float period = 32.f; // seconds

	accum_time = wrap_at_period(accum_time + dt, period);

	const float radius = main_radius;
	size_t index = 0;
	BBox contentBox;

	for (int i = 0; i < queue_length; ++i, index += 2) {
		const float angle = float(M_PI * 2.0) * (i / float(queue_length)) + accum_time * float(M_PI * 2.0) / period;
		const simd::f32x4 sin_iz = simd::sin(simd::f32x4(
			i / float(queue_length) * float(M_PI * 16.0),
			i / float(queue_length) * float(M_PI * 16.0) + accum_time * float(M_PI * 32.0) / period, 0.f, 0.f));
		const float sin_i = sin_iz[0];
		const float sin_z = sin_iz[1];

		{
			const matx3_rotate rot(angle, 0.f, 0.f, 1.f);
			const vect3 pos = vect3(radius + sin_i, 0.f, sin_z) * rot;
			const BBox box(
				pos - vect3(.5f),
				pos + vect3(.5f),
				BBox::flag_direct());

			contentBox.grow(box);
			content.getMutable(index + 0) = Voxel(box.get_min(), box.get_max());
		}
		{
			const matx3_rotate rot(-angle, 0.f, 0.f, 1.f);
			const vect3 pos = vect3(radius * .5f + sin_i, 0.f, sin_z) * rot;
			const BBox box(
				pos - vect3(.35f),
				pos + vect3(.35f),
				BBox::flag_direct());

			contentBox.grow(box);
			content.getMutable(index + 1) = Voxel(box.get_min(), box.get_max());
		}
	}

	// account for the last element which is never updated
	contentBox.grow(content.getElement(index).get_bbox());

	return scene.set_payload_array(content, contentBox);
}


bool Scene3::frame(
	Timeslice& scene,
	const float dt) {

	return update(scene, dt);
}


enum {
	scene_1,
	scene_2,
	scene_3,

	scene_count
};

////////////////////////////////////////////////////////////////////////////////
// the global control state
////////////////////////////////////////////////////////////////////////////////

namespace c {
	static const float beat_period = 1.714288; // seconds

	static size_t scene_selector;

	// scene properties
	static float decl = 0.f;
	static float azim = 0.f;

	// camera properties
	static float pos_x = 0.f;
	static float pos_y = 0.f;
	static float pos_z = 0.f;

	// accumulators
	static float accum_time = 0.f;
	static float accum_beat = 0.f;
	static float accum_beat_2 = 0.f;

	// view properties
	static float contrast_middle = .5f;
	static float contrast_k = 1.f;
	static float blur_split = -1.f;

} // namespace c

////////////////////////////////////////////////////////////////////////////////
// scripting support
////////////////////////////////////////////////////////////////////////////////

class Action {
protected:
	float lifespan;

public:
	// start the action; return true if action is alive after the start
	virtual bool start(
		const float delay,
		const float duration) = 0;

	// perform action at tick (ie. at frame); return true if action still alive
	virtual bool frame(
		const float dt) = 0;
};

////////////////////////////////////////////////////////////////////////////////
// ActionSetScene1: establish scene_1
////////////////////////////////////////////////////////////////////////////////

class ActionSetScene1 : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionSetScene1::start(
	const float,
	const float) {

	c::scene_selector = scene_1;

	c::pos_x = 0.f;
	c::pos_y = .25f;
	c::pos_z = .875f;

	c::decl = M_PI / -2.0;

	return false;
}


bool ActionSetScene1::frame(
	const float) {

	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionSetScene2: establish scene_2
////////////////////////////////////////////////////////////////////////////////

class ActionSetScene2 : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionSetScene2::start(
	const float,
	const float) {

	c::scene_selector = scene_2;

	c::pos_x = 0.f;
	c::pos_y = .25f;
	c::pos_z = 1.f;

	c::decl = M_PI / -2.0;

	return false;
}


bool ActionSetScene2::frame(
	const float) {

	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionSetScene3: establish scene_3
////////////////////////////////////////////////////////////////////////////////

class ActionSetScene3 : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionSetScene3::start(
	const float,
	const float) {

	c::scene_selector = scene_3;

	c::pos_x = 0.f;
	c::pos_y = 0.f;
	c::pos_z = 1.125f;

	c::decl = M_PI / -2.125;

	return false;
}


bool ActionSetScene3::frame(
	const float) {

	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionViewBlur: blur the view
////////////////////////////////////////////////////////////////////////////////

class ActionViewBlur : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionViewBlur::start(
	const float,
	const float) {

	c::blur_split = -1.f;

	return false;
}


bool ActionViewBlur::frame(
	const float) {

	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionViewUnblur: de-blur the view
////////////////////////////////////////////////////////////////////////////////

class ActionViewUnblur : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionViewUnblur::start(
	const float,
	const float) {

	c::blur_split = 1.f;

	return false;
}


bool ActionViewUnblur::frame(
	const float) {

	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionViewBlurDt: blur the view non-instantaneously
////////////////////////////////////////////////////////////////////////////////

class ActionViewBlurDt : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionViewBlurDt::start(
	const float delay,
	const float) {

	c::blur_split -= (2.0 / c::beat_period) * delay;

	if (c::blur_split > -1.f)
		return true;

	c::blur_split = -1.f;
	return false;
}


bool ActionViewBlurDt::frame(
	const float dt) {

	c::blur_split -= (2.0 / c::beat_period) * dt;

	if (c::blur_split > -1.f)
		return true;

	c::blur_split = -1.f;
	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionViewUnblurDt: de-blur the view non-instantaneously
////////////////////////////////////////////////////////////////////////////////

class ActionViewUnblurDt : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionViewUnblurDt::start(
	const float delay,
	const float) {

	c::blur_split += (2.0 / c::beat_period) * delay;

	if (c::blur_split < 1.f)
		return true;

	c::blur_split = 1.f;
	return false;
}


bool ActionViewUnblurDt::frame(
	const float dt) {

	c::blur_split += (2.0 / c::beat_period) * dt;

	if (c::blur_split < 1.f)
		return true;

	c::blur_split = 1.f;
	return false;
}

////////////////////////////////////////////////////////////////////////////////
// ActionViewSplit: split view in sync to the beat
////////////////////////////////////////////////////////////////////////////////

class ActionViewSplit : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionViewSplit::start(
	const float delay,
	const float duration) {

	if (delay >= duration)
		return false;

	c::blur_split = simd::sin(simd::f32x4(c::accum_beat * float(M_PI * 2.0 / c::beat_period), simd::flag_zero()))[0] * .25f;

	lifespan = duration - delay;
	return true;
}


bool ActionViewSplit::frame(
	const float dt) {

	if (dt >= lifespan)
		return false;

	c::blur_split = simd::sin(simd::f32x4(c::accum_beat * float(M_PI * 2.0 / c::beat_period), simd::flag_zero()))[0] * .25f;

	lifespan -= dt;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// ActionContrastBeat: pulse contrast to the beat
////////////////////////////////////////////////////////////////////////////////

class ActionContrastBeat : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionContrastBeat::start(
	const float delay,
	const float duration) {

	if (delay >= duration)
		return false;

	c::contrast_middle = .5f;
	c::contrast_k = 1.f + simd::pow(simd::sin(simd::abs(simd::f32x4(float(-M_PI_2) + float(M_PI) * c::accum_beat / c::beat_period, simd::flag_zero()))), simd::f32x4(64.f))[0];

	lifespan = duration - delay;
	return true;
}


bool ActionContrastBeat::frame(
	const float dt) {

	if (dt >= lifespan) {
		c::contrast_middle = .5f;
		c::contrast_k = 1.f;
		return false;
	}

	c::contrast_k = 1.f + simd::pow(simd::sin(simd::abs(simd::f32x4(float(-M_PI_2) + float(M_PI) * c::accum_beat / c::beat_period, simd::flag_zero()))), simd::f32x4(64.f))[0];

	lifespan -= dt;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// ActionCameraSnake: snake-style camera azimuth (singe quadrant only)
////////////////////////////////////////////////////////////////////////////////

class ActionCameraSnake : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionCameraSnake::start(
	const float delay,
	const float duration) {

	if (delay >= duration)
		return false;

	c::azim += simd::sin(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * float(M_PI / 4.0) * delay;

	lifespan = duration - delay;
	return true;
}


bool ActionCameraSnake::frame(
	const float dt) {

	if (dt >= lifespan)
		return false;

	c::azim += simd::sin(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * float(M_PI / 4.0) * dt;

	lifespan -= dt;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// ActionCameraBounce: snake-style camera azimuth (full range)
////////////////////////////////////////////////////////////////////////////////

class ActionCameraBounce : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionCameraBounce::start(
	const float delay,
	const float duration) {

	if (delay >= duration)
		return false;

	c::azim += simd::cos(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * float(M_PI / 4.0) * delay;

	lifespan = duration - delay;
	return true;
}


bool ActionCameraBounce::frame(
	const float dt) {

	if (dt >= lifespan)
		return false;

	c::azim += simd::cos(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * float(M_PI / 4.0) * dt;

	lifespan -= dt;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// ActionCameraBnF: camera position back'n'forth
////////////////////////////////////////////////////////////////////////////////

class ActionCameraBnF : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionCameraBnF::start(
	const float delay,
	const float duration) {

	if (delay >= duration)
		return false;

	c::pos_z += simd::cos(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * delay;

	lifespan = duration - delay;
	return true;
}


bool ActionCameraBnF::frame(
	const float dt) {

	if (dt >= lifespan)
		return false;

	c::pos_z += simd::cos(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * dt;

	lifespan -= dt;
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// ActionCameraLean: camera leaning forth then back
////////////////////////////////////////////////////////////////////////////////

class ActionCameraLean : virtual public Action {
public:
	// virtual from Action
	bool start(const float, const float);

	// virtual from Action
	bool frame(const float);
};


bool ActionCameraLean::start(
	const float delay,
	const float duration) {

	if (delay >= duration)
		return false;

	c::decl += simd::sin(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * delay;

	lifespan = duration - delay;
	return true;
}


bool ActionCameraLean::frame(
	const float dt) {

	if (dt >= lifespan)
		return false;

	c::decl += simd::sin(simd::f32x4(c::accum_beat_2 * float(M_PI / c::beat_period), simd::flag_zero()))[0] * dt;

	lifespan -= dt;
	return true;
}


cli_param param;
extern "C" unsigned input;
unsigned input;

namespace { // anonymous

Scene1 scene1;
Scene2 scene2;
Scene3 scene3;

Scene* const scene[] = {
	&scene1,
	&scene2,
	&scene3
};

ActionSetScene1    actionSetScene1;
ActionSetScene2    actionSetScene2;
ActionSetScene3    actionSetScene3;
ActionContrastBeat actionContrastBeat;
ActionViewBlur     actionViewBlur;
ActionViewUnblur   actionViewUnblur;
ActionViewBlurDt   actionViewBlurDt;
ActionViewUnblurDt actionViewUnblurDt;
ActionViewSplit    actionViewSplit;
ActionCameraSnake  actionCameraSnake;
ActionCameraBounce actionCameraBounce;
ActionCameraBnF    actionCameraBnF;
ActionCameraLean   actionCameraLean;

// master track of the application (entries sorted by start time)
const struct {
	const float start;    // seconds
	const float duration; // seconds
	Action& action;
}
track[] = {
	{   0.f,        0.f,      actionSetScene1 },
	{   0.f,        9.428584, actionViewSplit },
	{   0.f,       68.571342, actionContrastBeat },
	{   9.428584,   0.f,      actionViewBlurDt },
	{  68.571342,   0.f,      actionSetScene2 },
	{  68.571342,   0.f,      actionViewUnblur },
	{  68.571342,  27.426758, actionCameraBounce },
	{  82.285824,  6.8551240, actionCameraLean },
	{  95.998100,   0.f,      actionSetScene3 },
	{  95.998100,  39.430652, actionCameraSnake },
	{ 109.714432,  25.714320, actionCameraBnF },
	{ 133.714464,   0.f,      actionViewBlurDt },
	{ 135.428752,   0.f,      actionSetScene1 },
	{ 136.285896,  56.571504, actionContrastBeat },
	{ 198.f,        0.f,      actionViewUnblurDt }
};

const size_t octet_w = 2;
const size_t octet_h = 4096;

const size_t mem_size_octet = octet_w * octet_h * sizeof(cl_ushort4);
const size_t octet_count = mem_size_octet / sizeof(cl_ushort8);

const size_t leaf_w = 4;
const size_t leaf_h = 4096;

const size_t mem_size_leaf = leaf_w * leaf_h * sizeof(cl_ushort4);
const size_t leaf_count = mem_size_leaf / sizeof(cl_ushort16);

const size_t voxel_w = 2;
const size_t voxel_h = 4096;

const size_t mem_size_voxel = voxel_w * voxel_h * sizeof(cl_float4);
const size_t voxel_count = mem_size_voxel / sizeof(cl_float8);

const size_t carb_w = 1;
const size_t carb_h = 6;

const size_t mem_size_carb = carb_w * carb_h * sizeof(cl_float4);
const size_t carb_count = mem_size_carb / sizeof(cl_float4);

uint64_t tlast;
size_t frame_idx;

cl_ushort8* octet_map_buffer[n_buffering];
cl_ushort16* leaf_map_buffer[n_buffering];
cl_float8* voxel_map_buffer[n_buffering];
cl_float4* carb_map_buffer[n_buffering];
cl_uchar* image_map_buffer[n_buffering];

Array< Timeslice > timeline;

cl_command_queue persist_queue;
cl_kernel persist_kernel;

cl_event event_data_set_ready[n_buffering][4];
cl_event event_kernel_complete[n_buffering];

cl_mem src_a_d[n_buffering];
cl_mem src_b_d[n_buffering];
cl_mem src_c_d[n_buffering];
cl_mem src_d_d[n_buffering];
cl_mem dst_d[n_buffering];

const cl_uint work_dim = 2;
size_t global_ws[work_dim];
size_t local_ws[work_dim];

size_t track_cursor;
Action* action[8];
size_t action_count;

simd::f32x4 bbox_min;
simd::f32x4 bbox_max;
simd::f32x4 centre;
simd::f32x4 extent;
float max_extent;

} // namespace anonymous

int initFrame(void)
{
	using testbed::scoped_linkage_ptr;
	using testbed::scoped_ptr;
	using testbed::scoped_functor;
	using testbed::generic_free;
	using testbed::deinit_resources_t;

	using clutil::ocl_ver;
	using clutil::cldevice_version;
	using clutil::cldevice_type_single_string;
	using clutil::reportCLError;
	using clutil::reportCLCaps;

	const size_t image_w = param.image_w;
	const size_t image_h = param.image_h;

	if (!testbed::util::reportGLCaps() || !testbed::monv::init_resources(image_w, image_h)) {
		stream::cerr << "failed to initialise GL resources; bailing out\n";
		return 1;
	}

	// ensure GPU resources we just initialized get deinitialized at early-out errors
	scoped_ptr< deinit_resources_t, scoped_functor > release_monv(testbed::monv::deinit_resources);

	const bool report_caps                  = bool(param.flags & cli_param::BIT_REPORT_CAPS);
	const bool discard_platform_version     = bool(param.flags & cli_param::BIT_DISCARD_PLATFORM_VERSION);
	const bool discard_device_version       = bool(param.flags & cli_param::BIT_DISCARD_DEVICE_VERSION);
	const KERNEL_PARAM_TYPE kern_param_type = param.kern_param_type;
	const unsigned platform_idx             = param.platform_idx;
	unsigned device_idx                     = param.device_idx;

	if (report_caps) {
		const int result_caps = reportCLCaps(discard_platform_version, discard_device_version);

		if (0 != result_caps)
			return result_caps;
	}

	cl_uint num_platforms = 0;
	cl_int success = clGetPlatformIDs(num_platforms, 0, &num_platforms);

	if (reportCLError(success)) {
		stream::cerr << "failure at clGetPlatformIDs; terminate\n";
		return -1;
	}

	if (platform_idx >= unsigned(num_platforms)) {
		stream::cerr << "requested platform index is out of bounds\n";
		return 1;
	}

	const scoped_ptr< cl_platform_id, generic_free > platform(
		reinterpret_cast< cl_platform_id* >(std::malloc(sizeof(cl_platform_id) * num_platforms)));

	success = clGetPlatformIDs(num_platforms, platform(), 0);

	if (reportCLError(success)) {
		stream::cerr << "failure at clGetPlatformIDs; terminate\n";
		return -1;
	}

	scoped_ptr< cl_device_id, generic_free > device;

	if (-1U != device_idx) {
		cl_uint num_devices = 0;

		success = clGetDeviceIDs(platform()[platform_idx],
			CL_DEVICE_TYPE_ALL, num_devices, 0, &num_devices);

		if (reportCLError(success)) {
			stream::cerr << "failure at clGetDeviceIDs; terminate\n";
			return -1;
		}

		if (device_idx >= unsigned(num_devices)) {
			stream::cerr << "requested device index is out of bound\n";
			return 1;
		}

		scoped_ptr< cl_device_id, generic_free > proto_device(
			reinterpret_cast< cl_device_id* >(std::malloc(sizeof(cl_device_id) * num_devices)));

		success = clGetDeviceIDs(platform()[platform_idx],
			CL_DEVICE_TYPE_ALL, num_devices, proto_device(), 0);

		if (reportCLError(success)) {
			stream::cerr << "failure at clGetDeviceIDs; terminate\n";
			return -1;
		}

		device.swap(proto_device);
	}
	else {
		const cl_device_type devtype = CL_DEVICE_TYPE_GPU;

		scoped_ptr< cl_device_id, generic_free > proto_device(
			reinterpret_cast< cl_device_id* >(std::malloc(sizeof(cl_device_id))));

		success = clGetDeviceIDs(platform()[platform_idx],
			devtype, 1, proto_device(), 0);

		if (reportCLError(success)) {
			stream::cerr << "error getting device of type " <<
				cldevice_type_single_string(devtype) << '\n';
			return -1;
		}

		device_idx = 0;
		device.swap(proto_device);
	}

	ocl_ver platform_version, device_version;

	if (!clplatform_version(platform()[platform_idx], platform_version) ||
		!cldevice_version(device()[device_idx], device_version)) {

		stream::cerr << "error getting platform/device version\n";
		return -1;
	}

	if (!clwrap::init(platform_version)) {
		stream::cerr << "error initializing API wrappers\n";
		return -1;
	}

	const cl_context_properties prop[] = {
		CL_CONTEXT_PLATFORM, cl_context_properties(platform()[platform_idx]),
		0
	};

	cl_context context = clCreateContext(prop, 1, device() + device_idx, 0, 0, &success);

	if (reportCLError(success)) {
		stream::cerr << "error creating context\n";
		return -1;
	}

	const scoped_ptr< cl_context, scoped_functor > release_context(&context);

	cl_command_queue queue = clCreateCommandQueue(context,
		device()[device_idx], CL_QUEUE_PROFILING_ENABLE, &success);

	if (reportCLError(success)) {
		stream::cerr << "error creating command queue\n";
		return -1;
	}

	scoped_ptr< cl_command_queue, scoped_functor > release_queue(&queue);

	// octet map element:
	// struct Octet {
	//     OctetId child[8]; // OctetId := ushort
	// }
	// represent in image as: ushort4 child[2]
	scoped_ptr< cl_ushort8, generic_free > octet_map(
		reinterpret_cast< cl_ushort8* >(std::malloc(mem_size_octet * n_buffering)));

	if (0 == octet_map()) {
		stream::cerr << "error allocating octet_map\n";
		return -1;
	}

	for (size_t i = 0; i < sizeof(octet_map_buffer) / sizeof(octet_map_buffer[0]); ++i)
		octet_map_buffer[i] = octet_map() + i * octet_count;

	// leaf map element:
	// struct Leaf {
	//     PayloadId start[8]; // PayloadId := ushort
	//     PayloadId count[8]; // PayloadId := ushort
	// }
	// represent in image as: ushort4 start[2], ushort4 count[2]
	scoped_ptr< cl_ushort16, generic_free > leaf_map(
		reinterpret_cast< cl_ushort16* >(std::malloc(mem_size_leaf * n_buffering)));

	if (0 == leaf_map()) {
		stream::cerr << "error allocating leaf_map\n";
		return -1;
	}

	for (size_t i = 0; i < sizeof(leaf_map_buffer) / sizeof(leaf_map_buffer[0]); ++i)
		leaf_map_buffer[i] = leaf_map() + i * leaf_count;

	// voxel map element:
	// struct BBox {
	//     float min[3];
	//     uint32 min_cookie;
	//     float max[3];
	//     uint32 max_cookie;
	// }
	// represent in image as: float4 min, float4 max
	scoped_ptr< cl_float8, generic_free > voxel_map(
		reinterpret_cast< cl_float8* >(std::malloc(mem_size_voxel * n_buffering)));

	if (0 == voxel_map()) {
		stream::cerr << "error allocating voxel_map\n";
		return -1;
	}

	for (size_t i = 0; i < sizeof(voxel_map_buffer) / sizeof(voxel_map_buffer[0]); ++i)
		voxel_map_buffer[i] = voxel_map() + i * voxel_count;

	// CARB (camera and root bbox) constant buffer element:
	// float4
	scoped_ptr< cl_float4, generic_free > carb_map(
		reinterpret_cast< cl_float4* >(std::malloc(mem_size_carb * n_buffering)));

	if (0 == carb_map()) {
		stream::cerr << "error allocating carb_map\n";
		return -1;
	}

	for (size_t i = 0; i < sizeof(carb_map_buffer) / sizeof(carb_map_buffer[0]); ++i)
		carb_map_buffer[i] = carb_map() + i * carb_count;

	// create cl mem objects ///////////////////////////////////////////////////
	const cl_image_format src_image_format_ushort4 = {
		CL_RGBA,				// cl_channel_order image_channel_order;
		CL_UNSIGNED_INT16		// cl_channel_type image_channel_data_type;
	};

	const cl_image_format src_image_format_float4 = {
		CL_RGBA,				// cl_channel_order image_channel_order;
		CL_FLOAT				// cl_channel_type image_channel_data_type;
	};

	const cl_image_desc src_image_desc_octet = {
		CL_MEM_OBJECT_IMAGE2D,	// cl_mem_object_type image_type;
		octet_w,				// size_t image_width;
		octet_h,				// size_t image_height;
		size_t(0),				// size_t image_depth;
		size_t(0),				// size_t image_array_size;
		size_t(0),				// size_t image_row_pitch;
		size_t(0),				// size_t image_slice_pitch;
		cl_uint(0),				// cl_uint num_mip_levels;
		cl_uint(0),				// cl_uint num_samples;
		cl_mem(0)				// cl_mem buffer;
	};

	const cl_image_desc src_image_desc_leaf = {
		CL_MEM_OBJECT_IMAGE2D,	// cl_mem_object_type image_type;
		leaf_w,					// size_t image_width;
		leaf_h,					// size_t image_height;
		size_t(0),				// size_t image_depth;
		size_t(0),				// size_t image_array_size;
		size_t(0),				// size_t image_row_pitch;
		size_t(0),				// size_t image_slice_pitch;
		cl_uint(0),				// cl_uint num_mip_levels;
		cl_uint(0),				// cl_uint num_samples;
		cl_mem(0)				// cl_mem buffer;
	};

	const cl_image_desc src_image_desc_voxel = {
		CL_MEM_OBJECT_IMAGE2D,	// cl_mem_object_type image_type;
		voxel_w,				// size_t image_width;
		voxel_h,				// size_t image_height;
		size_t(0),				// size_t image_depth;
		size_t(0),				// size_t image_array_size;
		size_t(0),				// size_t image_row_pitch;
		size_t(0),				// size_t image_slice_pitch;
		cl_uint(0),				// cl_uint num_mip_levels;
		cl_uint(0),				// cl_uint num_samples;
		cl_mem(0)				// cl_mem buffer;
	};

	cl_mem_flags src_flags = CL_MEM_READ_ONLY;

	if (device_version >= ocl_ver(1, 2))
		src_flags |= CL_MEM_HOST_WRITE_ONLY;

	// octet map source ///////////////////////////////////////////////////////////
	scoped_ptr< cl_mem, scoped_functor > release_src_a(src_a_d);

	switch (kern_param_type) {
	case PARAM_TYPE_BUFFER:
		for (size_t i = 0; i < sizeof(src_a_d) / sizeof(src_a_d[0]); ++i) {
			src_a_d[i] = clCreateBuffer(context, src_flags, mem_size_octet, 0, &success);

			if (reportCLError(success)) {
				stream::cerr << "error creating buffer for src_a\n";
				return -1;
			}
		}
		break;

	case PARAM_TYPE_IMAGE:
		for (size_t i = 0; i < sizeof(src_a_d) / sizeof(src_a_d[0]); ++i) {
			src_a_d[i] = clwrap::clCreateImage(context,
				src_flags, &src_image_format_ushort4, &src_image_desc_octet, 0, &success);

			if (reportCLError(success)) {
				stream::cerr << "error creating image for src_a\n";
				return -1;
			}
		}
		break;
	}

	// leaf map source ////////////////////////////////////////////////////////////
	scoped_ptr< cl_mem, scoped_functor > release_src_b(src_b_d);

	switch (kern_param_type) {
	case PARAM_TYPE_BUFFER:
		for (size_t i = 0; i < sizeof(src_b_d) / sizeof(src_b_d[0]); ++i) {
			src_b_d[i] = clCreateBuffer(context, src_flags, mem_size_leaf, 0, &success);

			if (reportCLError(success)) {
				stream::cerr << "error creating buffer for src_b\n";
				return -1;
			}
		}
		break;

	case PARAM_TYPE_IMAGE:
		for (size_t i = 0; i < sizeof(src_b_d) / sizeof(src_b_d[0]); ++i) {
			src_b_d[i] = clwrap::clCreateImage(context,
				src_flags, &src_image_format_ushort4, &src_image_desc_leaf, 0, &success);

			if (reportCLError(success)) {
				stream::cerr << "error creating image for src_b\n";
				return -1;
			}
		}
		break;
	}

	// voxel map source ///////////////////////////////////////////////////////////
	scoped_ptr< cl_mem, scoped_functor > release_src_c(src_c_d);

	switch (kern_param_type) {
	case PARAM_TYPE_BUFFER:
		for (size_t i = 0; i < sizeof(src_c_d) / sizeof(src_c_d[0]); ++i) {
			src_c_d[i] = clCreateBuffer(context, src_flags, mem_size_voxel, 0, &success);

			if (reportCLError(success)) {
				stream::cerr << "error creating buffer for src_c\n";
				return -1;
			}
		}
		break;

	case PARAM_TYPE_IMAGE:
		for (size_t i = 0; i < sizeof(src_c_d) / sizeof(src_c_d[0]); ++i) {
			src_c_d[i] = clwrap::clCreateImage(context,
				src_flags, &src_image_format_float4, &src_image_desc_voxel, 0, &success);

			if (reportCLError(success)) {
				stream::cerr << "error creating image for src_c\n";
				return -1;
			}
		}
		break;
	}

	// cam mem source /////////////////////////////////////////////////////////////
	scoped_ptr< cl_mem, scoped_functor > release_src_d(src_d_d);

	for (size_t i = 0; i < sizeof(src_d_d) / sizeof(src_d_d[0]); ++i) {
		src_d_d[i] = clCreateBuffer(context, src_flags, mem_size_carb, 0, &success);

		if (reportCLError(success)) {
			stream::cerr << "error creating buffer for src_d\n";
			return -1;
		}
	}

	// output map /////////////////////////////////////////////////////////////////
	const size_t mem_size_image = image_w * image_h * sizeof(cl_uchar);
	const size_t pixel_count = image_w * image_h;

	scoped_ptr< cl_uchar, generic_free > image_map(
		reinterpret_cast< cl_uchar* >(std::malloc(mem_size_image * n_buffering)));

	if (0 == image_map()) {
		stream::cerr << "error allocating image_map\n";
		return -1;
	}

	for (size_t i = 0; i < sizeof(image_map_buffer) / sizeof(image_map_buffer[0]); ++i)
		image_map_buffer[i] = image_map() + i * pixel_count;

	cl_mem_flags dst_flags = CL_MEM_WRITE_ONLY;

	if (device_version >= ocl_ver(1, 2))
		dst_flags |= CL_MEM_HOST_READ_ONLY;

	scoped_ptr< cl_mem, scoped_functor > release_dst(dst_d);

	for (size_t i = 0; i < sizeof(dst_d) / sizeof(dst_d[0]); ++i) {
		dst_d[i] = clCreateBuffer(context, dst_flags, mem_size_image, 0, &success);

		if (reportCLError(success)) {
			stream::cerr << "error creating buffer for result\n";
			return -1;
		}
	}

	// kernel sources /////////////////////////////////////////////////////////////
	const char source_prologue[] =
		"struct BBox {\n"
		"	float3 min;\n"
		"	float3 max;\n"
		"};\n"
		"\n"
		"struct Ray {\n"
		"	float3 origin;\n"
		"	float3 rcpdir;\n"
		"};\n"
		"\n"
		"struct Octet {\n"
		"	ushort8 child;\n"
		"};\n"
		"\n"
		"struct Leaf {\n"
		"	ushort8 start;\n"
		"	ushort8 count;\n"
		"};\n"
		"\n"
		"struct Voxel {\n"
		"	float4 min;\n"
		"	float4 max;\n"
		"};\n"
		"\n"
		"struct ChildIndex {\n"
		"	float8 distance;\n"
		"	uint8 index;\n"
		"};\n"
		"\n"
		"inline float intersect(\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray)\n"
		"{\n"
		"	const float3 t0 = (bbox->min - ray->origin) * ray->rcpdir;\n"
		"	const float3 t1 = (bbox->max - ray->origin) * ray->rcpdir;\n"
		"\n"
		"	const float3 axial_min = fmin(t0, t1);\n"
		"	const float3 axial_max = fmax(t0, t1);\n"
		"\n"
		"	const float min = fmax(fmax(axial_min.x, axial_min.y), axial_min.z);\n"
		"	const float max = fmin(fmin(axial_max.x, axial_max.y), axial_max.z);\n"
		"\n"
		"	const int2 msk = isless((float2)(0.f, min), (float2)(min, max));\n"
		"	return select(INFINITY, min, msk.x & msk.y);\n"
		"}\n"
		"\n"
		"inline void intersect8(\n"
		"	const float8 bbox_min_x,\n"
		"	const float8 bbox_min_y,\n"
		"	const float8 bbox_min_z,\n"
		"	const float8 bbox_max_x,\n"
		"	const float8 bbox_max_y,\n"
		"	const float8 bbox_max_z,\n"
		"	const struct Ray* const ray,\n"
		"	float8* const t,\n"
		"	int8* const r)\n"
		"{\n"
		"	const float3 ray_origin = ray->origin;\n"
		"	const float3 ray_rcpdir = ray->rcpdir;\n"
		"\n"
		"	const float8 tmin_x = (bbox_min_x - ray_origin.xxxxxxxx) * ray_rcpdir.xxxxxxxx;\n"
		"	const float8 tmax_x = (bbox_max_x - ray_origin.xxxxxxxx) * ray_rcpdir.xxxxxxxx;\n"
		"	const float8 tmin_y = (bbox_min_y - ray_origin.yyyyyyyy) * ray_rcpdir.yyyyyyyy;\n"
		"	const float8 tmax_y = (bbox_max_y - ray_origin.yyyyyyyy) * ray_rcpdir.yyyyyyyy;\n"
		"	const float8 tmin_z = (bbox_min_z - ray_origin.zzzzzzzz) * ray_rcpdir.zzzzzzzz;\n"
		"	const float8 tmax_z = (bbox_max_z - ray_origin.zzzzzzzz) * ray_rcpdir.zzzzzzzz;\n"
		"\n"
		"	const float8 x_min = fmin(tmin_x, tmax_x);\n"
		"	const float8 x_max = fmax(tmin_x, tmax_x);\n"
		"	const float8 y_min = fmin(tmin_y, tmax_y);\n"
		"	const float8 y_max = fmax(tmin_y, tmax_y);\n"
		"	const float8 z_min = fmin(tmin_z, tmax_z);\n"
		"	const float8 z_max = fmax(tmin_z, tmax_z);\n"
		"\n"
		"	const float8 min = fmax(fmax(x_min, y_min), z_min);\n"
		"	const float8 max = fmin(fmin(x_max, y_max), z_max);\n"
		"	*t = max;\n"
		"\n"
		"	const int8 msk = isless(min, max) & isless((float8)(0.f), max);\n"
		"	*r = msk;\n"
		"}\n"
		"\n"
		"uint octlf_intersect_wide(\n"
		"	const struct Leaf octet,\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray,\n"
		"	struct ChildIndex* const child_index)\n"
		"{\n"
		"	const float3 par_min = bbox->min;\n"
		"	const float3 par_max = bbox->max;\n"
		"	const float3 par_mid = (par_min + par_max) * 0.5f;\n"
		"\n"
		"	const float8 bbox_min_x = (float8)( par_min.x, par_mid.x, par_min.x, par_mid.x, par_min.x, par_mid.x, par_min.x, par_mid.x );\n"
		"	const float8 bbox_min_y = (float8)( par_min.yy, par_mid.yy, par_min.yy, par_mid.yy );\n"
		"	const float8 bbox_min_z = (float8)( par_min.zzzz, par_mid.zzzz );\n"
		"	const float8 bbox_max_x = (float8)( par_mid.x, par_max.x, par_mid.x, par_max.x, par_mid.x, par_max.x, par_mid.x, par_max.x );\n"
		"	const float8 bbox_max_y = (float8)( par_mid.yy, par_max.yy, par_mid.yy, par_max.yy );\n"
		"	const float8 bbox_max_z = (float8)( par_mid.zzzz, par_max.zzzz );\n"
		"\n"
		"	float8 t;\n"
		"	int8 r;\n"
		"	intersect8(bbox_min_x, bbox_min_y, bbox_min_z, bbox_max_x, bbox_max_y, bbox_max_z, ray, &t, &r);\n"
		"	const int8 occupancy = convert_int8((ushort8)(0) != octet.count);\n"
		"	r &= occupancy;\n"
		"\n"
		"	int count = 0;\n"
		"	count -= r.s0;\n"
		"	count -= r.s1;\n"
		"	count -= r.s2;\n"
		"	count -= r.s3;\n"
		"	count -= r.s4;\n"
		"	count -= r.s5;\n"
		"	count -= r.s6;\n"
		"	count -= r.s7;\n"
		"	t = select((float8)(INFINITY), t, r);\n"
		"\n"
		"	const float4 r0_A = (float4)(t.s0, t.s3, t.s4, t.s7);\n"
		"	const float4 r0_B = (float4)(t.s1, t.s2, t.s5, t.s6);\n"
		"	const int4 r0x_A = (int4)(0, 3, 4, 7);\n"
		"	const int4 r0x_B = (int4)(1, 2, 5, 6);\n"
		"	const int4 m0 = islessequal(r0_A, r0_B);\n"
		"	const float4 r0_min = fmin(r0_A, r0_B);\n"
		"	const float4 r0_max = fmax(r0_A, r0_B);\n"
		"	const int4 r0x_min = select(r0x_B, r0x_A, m0);\n"
		"	const int4 r0x_max = select(r0x_A, r0x_B, m0);\n"
		"\n"
		"	const float4 r1_A = (float4)(r0_min.s0, r0_max.s0, r0_max.s3, r0_min.s3);\n"
		"	const float4 r1_B = (float4)(r0_max.s1, r0_min.s1, r0_min.s2, r0_max.s2);\n"
		"	const int4 r1x_A = (int4)(r0x_min.s0, r0x_max.s0, r0x_max.s3, r0x_min.s3);\n"
		"	const int4 r1x_B = (int4)(r0x_max.s1, r0x_min.s1, r0x_min.s2, r0x_max.s2);\n"
		"	const int4 m1 = islessequal(r1_A, r1_B);\n"
		"	const float4 r1_min = fmin(r1_A, r1_B);\n"
		"	const float4 r1_max = fmax(r1_A, r1_B);\n"
		"	const int4 r1x_min = select(r1x_B, r1x_A, m1);\n"
		"	const int4 r1x_max = select(r1x_A, r1x_B, m1);\n"
		"\n"
		"	const float4 r2_A = (float4)(r1_min.s0, r1_max.s0, r1_max.s3, r1_min.s3);\n"
		"	const float4 r2_B = (float4)(r1_min.s1, r1_max.s1, r1_max.s2, r1_min.s2);\n"
		"	const int4 r2x_A = (int4)(r1x_min.s0, r1x_max.s0, r1x_max.s3, r1x_min.s3);\n"
		"	const int4 r2x_B = (int4)(r1x_min.s1, r1x_max.s1, r1x_max.s2, r1x_min.s2);\n"
		"	const int4 m2 = islessequal(r2_A, r2_B);\n"
		"	const float4 r2_min = fmin(r2_A, r2_B);\n"
		"	const float4 r2_max = fmax(r2_A, r2_B);\n"
		"	const int4 r2x_min = select(r2x_B, r2x_A, m2);\n"
		"	const int4 r2x_max = select(r2x_A, r2x_B, m2);\n"
		"\n"
		"	const float4 r3_A = (float4)(r2_min.s0, r2_max.s0, r2_min.s1, r2_max.s1);\n"
		"	const float4 r3_B = (float4)(r2_max.s2, r2_min.s2, r2_max.s3, r2_min.s3);\n"
		"	const int4 r3x_A = (int4)(r2x_min.s0, r2x_max.s0, r2x_min.s1, r2x_max.s1);\n"
		"	const int4 r3x_B = (int4)(r2x_max.s2, r2x_min.s2, r2x_max.s3, r2x_min.s3);\n"
		"	const int4 m3 = islessequal(r3_A, r3_B);\n"
		"	const float4 r3_min = fmin(r3_A, r3_B);\n"
		"	const float4 r3_max = fmax(r3_A, r3_B);\n"
		"	const int4 r3x_min = select(r3x_B, r3x_A, m3);\n"
		"	const int4 r3x_max = select(r3x_A, r3x_B, m3);\n"
		"\n"
		"	const float4 r4_A = (float4)(r3_min.s0, r3_min.s1, r3_max.s0, r3_max.s1);\n"
		"	const float4 r4_B = (float4)(r3_min.s2, r3_min.s3, r3_max.s2, r3_max.s3);\n"
		"	const int4 r4x_A = (int4)(r3x_min.s0, r3x_min.s1, r3x_max.s0, r3x_max.s1);\n"
		"	const int4 r4x_B = (int4)(r3x_min.s2, r3x_min.s3, r3x_max.s2, r3x_max.s3);\n"
		"	const int4 m4 = islessequal(r4_A, r4_B);\n"
		"	const float4 r4_min = fmin(r4_A, r4_B);\n"
		"	const float4 r4_max = fmax(r4_A, r4_B);\n"
		"	const int4 r4x_min = select(r4x_B, r4x_A, m4);\n"
		"	const int4 r4x_max = select(r4x_A, r4x_B, m4);\n"
		"\n"
		"	const float4 r5_A = (float4)(r4_min.s0, r4_max.s0, r4_min.s2, r4_max.s2);\n"
		"	const float4 r5_B = (float4)(r4_min.s1, r4_max.s1, r4_min.s3, r4_max.s3);\n"
		"	const int4 r5x_A = (int4)(r4x_min.s0, r4x_max.s0, r4x_min.s2, r4x_max.s2);\n"
		"	const int4 r5x_B = (int4)(r4x_min.s1, r4x_max.s1, r4x_min.s3, r4x_max.s3);\n"
		"	const int4 m5 = islessequal(r5_A, r5_B);\n"
		"	const float4 r5_min = fmin(r5_A, r5_B);\n"
		"	const float4 r5_max = fmax(r5_A, r5_B);\n"
		"	const int4 r5x_min = select(r5x_B, r5x_A, m5);\n"
		"	const int4 r5x_max = select(r5x_A, r5x_B, m5);\n"
		"\n"
		"	child_index->distance = (float8)(r5_min.s0, r5_max.s0, r5_min.s1, r5_max.s1, r5_min.s2, r5_max.s2, r5_min.s3, r5_max.s3);\n"
		"	child_index->index = (uint8)(r5x_min.s0, r5x_max.s0, r5x_min.s1, r5x_max.s1, r5x_min.s2, r5x_max.s2, r5x_min.s3, r5x_max.s3);\n"
		"	return as_uint(count);\n"
		"}\n"
		"\n"
		"uint octet_intersect_wide(\n"
		"	const struct Octet octet,\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray,\n"
		"	struct ChildIndex* const child_index,\n"
		"	struct BBox child_bbox[8])\n"
		"{\n"
		"	const float3 par_min = bbox->min;\n"
		"	const float3 par_max = bbox->max;\n"
		"	const float3 par_mid = (par_min + par_max) * 0.5f;\n"
		"\n"
		"	const float8 bbox_min_x = (float8)( par_min.x, par_mid.x, par_min.x, par_mid.x, par_min.x, par_mid.x, par_min.x, par_mid.x );\n"
		"	const float8 bbox_min_y = (float8)( par_min.yy, par_mid.yy, par_min.yy, par_mid.yy );\n"
		"	const float8 bbox_min_z = (float8)( par_min.zzzz, par_mid.zzzz );\n"
		"	const float8 bbox_max_x = (float8)( par_mid.x, par_max.x, par_mid.x, par_max.x, par_mid.x, par_max.x, par_mid.x, par_max.x );\n"
		"	const float8 bbox_max_y = (float8)( par_mid.yy, par_max.yy, par_mid.yy, par_max.yy );\n"
		"	const float8 bbox_max_z = (float8)( par_mid.zzzz, par_max.zzzz );\n"
		"\n"
		"	child_bbox[0] = (struct BBox){ (float3)( bbox_min_x.s0, bbox_min_y.s0, bbox_min_z.s0 ), (float3)( bbox_max_x.s0, bbox_max_y.s0, bbox_max_z.s0 ) };\n"
		"	child_bbox[1] = (struct BBox){ (float3)( bbox_min_x.s1, bbox_min_y.s1, bbox_min_z.s1 ), (float3)( bbox_max_x.s1, bbox_max_y.s1, bbox_max_z.s1 ) };\n"
		"	child_bbox[2] = (struct BBox){ (float3)( bbox_min_x.s2, bbox_min_y.s2, bbox_min_z.s2 ), (float3)( bbox_max_x.s2, bbox_max_y.s2, bbox_max_z.s2 ) };\n"
		"	child_bbox[3] = (struct BBox){ (float3)( bbox_min_x.s3, bbox_min_y.s3, bbox_min_z.s3 ), (float3)( bbox_max_x.s3, bbox_max_y.s3, bbox_max_z.s3 ) };\n"
		"	child_bbox[4] = (struct BBox){ (float3)( bbox_min_x.s4, bbox_min_y.s4, bbox_min_z.s4 ), (float3)( bbox_max_x.s4, bbox_max_y.s4, bbox_max_z.s4 ) };\n"
		"	child_bbox[5] = (struct BBox){ (float3)( bbox_min_x.s5, bbox_min_y.s5, bbox_min_z.s5 ), (float3)( bbox_max_x.s5, bbox_max_y.s5, bbox_max_z.s5 ) };\n"
		"	child_bbox[6] = (struct BBox){ (float3)( bbox_min_x.s6, bbox_min_y.s6, bbox_min_z.s6 ), (float3)( bbox_max_x.s6, bbox_max_y.s6, bbox_max_z.s6 ) };\n"
		"	child_bbox[7] = (struct BBox){ (float3)( bbox_min_x.s7, bbox_min_y.s7, bbox_min_z.s7 ), (float3)( bbox_max_x.s7, bbox_max_y.s7, bbox_max_z.s7 ) };\n"
		"\n"
		"	float8 t;\n"
		"	int8 r;\n"
		"	intersect8(bbox_min_x, bbox_min_y, bbox_min_z, bbox_max_x, bbox_max_y, bbox_max_z, ray, &t, &r);\n"
		"	const int8 occupancy = convert_int8((ushort8)(-1) != octet.child);\n"
		"	r &= occupancy;\n"
		"\n"
		"	int count = 0;\n"
		"	count -= r.s0;\n"
		"	count -= r.s1;\n"
		"	count -= r.s2;\n"
		"	count -= r.s3;\n"
		"	count -= r.s4;\n"
		"	count -= r.s5;\n"
		"	count -= r.s6;\n"
		"	count -= r.s7;\n"
		"	t = select((float8)(INFINITY), t, r);\n"
		"\n"
		"	const float4 r0_A = (float4)(t.s0, t.s3, t.s4, t.s7);\n"
		"	const float4 r0_B = (float4)(t.s1, t.s2, t.s5, t.s6);\n"
		"	const int4 r0x_A = (int4)(0, 3, 4, 7);\n"
		"	const int4 r0x_B = (int4)(1, 2, 5, 6);\n"
		"	const int4 m0 = islessequal(r0_A, r0_B);\n"
		"	const float4 r0_min = fmin(r0_A, r0_B);\n"
		"	const float4 r0_max = fmax(r0_A, r0_B);\n"
		"	const int4 r0x_min = select(r0x_B, r0x_A, m0);\n"
		"	const int4 r0x_max = select(r0x_A, r0x_B, m0);\n"
		"\n"
		"	const float4 r1_A = (float4)(r0_min.s0, r0_max.s0, r0_max.s3, r0_min.s3);\n"
		"	const float4 r1_B = (float4)(r0_max.s1, r0_min.s1, r0_min.s2, r0_max.s2);\n"
		"	const int4 r1x_A = (int4)(r0x_min.s0, r0x_max.s0, r0x_max.s3, r0x_min.s3);\n"
		"	const int4 r1x_B = (int4)(r0x_max.s1, r0x_min.s1, r0x_min.s2, r0x_max.s2);\n"
		"	const int4 m1 = islessequal(r1_A, r1_B);\n"
		"	const float4 r1_min = fmin(r1_A, r1_B);\n"
		"	const float4 r1_max = fmax(r1_A, r1_B);\n"
		"	const int4 r1x_min = select(r1x_B, r1x_A, m1);\n"
		"	const int4 r1x_max = select(r1x_A, r1x_B, m1);\n"
		"\n"
		"	const float4 r2_A = (float4)(r1_min.s0, r1_max.s0, r1_max.s3, r1_min.s3);\n"
		"	const float4 r2_B = (float4)(r1_min.s1, r1_max.s1, r1_max.s2, r1_min.s2);\n"
		"	const int4 r2x_A = (int4)(r1x_min.s0, r1x_max.s0, r1x_max.s3, r1x_min.s3);\n"
		"	const int4 r2x_B = (int4)(r1x_min.s1, r1x_max.s1, r1x_max.s2, r1x_min.s2);\n"
		"	const int4 m2 = islessequal(r2_A, r2_B);\n"
		"	const float4 r2_min = fmin(r2_A, r2_B);\n"
		"	const float4 r2_max = fmax(r2_A, r2_B);\n"
		"	const int4 r2x_min = select(r2x_B, r2x_A, m2);\n"
		"	const int4 r2x_max = select(r2x_A, r2x_B, m2);\n"
		"\n"
		"	const float4 r3_A = (float4)(r2_min.s0, r2_max.s0, r2_min.s1, r2_max.s1);\n"
		"	const float4 r3_B = (float4)(r2_max.s2, r2_min.s2, r2_max.s3, r2_min.s3);\n"
		"	const int4 r3x_A = (int4)(r2x_min.s0, r2x_max.s0, r2x_min.s1, r2x_max.s1);\n"
		"	const int4 r3x_B = (int4)(r2x_max.s2, r2x_min.s2, r2x_max.s3, r2x_min.s3);\n"
		"	const int4 m3 = islessequal(r3_A, r3_B);\n"
		"	const float4 r3_min = fmin(r3_A, r3_B);\n"
		"	const float4 r3_max = fmax(r3_A, r3_B);\n"
		"	const int4 r3x_min = select(r3x_B, r3x_A, m3);\n"
		"	const int4 r3x_max = select(r3x_A, r3x_B, m3);\n"
		"\n"
		"	const float4 r4_A = (float4)(r3_min.s0, r3_min.s1, r3_max.s0, r3_max.s1);\n"
		"	const float4 r4_B = (float4)(r3_min.s2, r3_min.s3, r3_max.s2, r3_max.s3);\n"
		"	const int4 r4x_A = (int4)(r3x_min.s0, r3x_min.s1, r3x_max.s0, r3x_max.s1);\n"
		"	const int4 r4x_B = (int4)(r3x_min.s2, r3x_min.s3, r3x_max.s2, r3x_max.s3);\n"
		"	const int4 m4 = islessequal(r4_A, r4_B);\n"
		"	const float4 r4_min = fmin(r4_A, r4_B);\n"
		"	const float4 r4_max = fmax(r4_A, r4_B);\n"
		"	const int4 r4x_min = select(r4x_B, r4x_A, m4);\n"
		"	const int4 r4x_max = select(r4x_A, r4x_B, m4);\n"
		"\n"
		"	const float4 r5_A = (float4)(r4_min.s0, r4_max.s0, r4_min.s2, r4_max.s2);\n"
		"	const float4 r5_B = (float4)(r4_min.s1, r4_max.s1, r4_min.s3, r4_max.s3);\n"
		"	const int4 r5x_A = (int4)(r4x_min.s0, r4x_max.s0, r4x_min.s2, r4x_max.s2);\n"
		"	const int4 r5x_B = (int4)(r4x_min.s1, r4x_max.s1, r4x_min.s3, r4x_max.s3);\n"
		"	const int4 m5 = islessequal(r5_A, r5_B);\n"
		"	const float4 r5_min = fmin(r5_A, r5_B);\n"
		"	const float4 r5_max = fmax(r5_A, r5_B);\n"
		"	const int4 r5x_min = select(r5x_B, r5x_A, m5);\n"
		"	const int4 r5x_max = select(r5x_A, r5x_B, m5);\n"
		"\n"
		"	child_index->distance = (float8)(r5_min.s0, r5_max.s0, r5_min.s1, r5_max.s1, r5_min.s2, r5_max.s2, r5_min.s3, r5_max.s3);\n"
		"	child_index->index = (uint8)(r5x_min.s0, r5x_max.s0, r5x_min.s1, r5x_max.s1, r5x_min.s2, r5x_max.s2, r5x_min.s3, r5x_max.s3);\n"
		"	return as_uint(count);\n"
		"}\n";
	const char source_buffer[] =
		"inline struct Octet get_octet(\n"
		"	__global const ushort4* const octet,\n"
		"	const uint idx)\n"
		"{\n"
		"	return (struct Octet){\n"
		"		(ushort8)(\n"
		"			octet[idx * 2 + 0],\n"
		"			octet[idx * 2 + 1]\n"
		"		)\n"
		"	};\n"
		"}\n"
		"inline struct Leaf get_leaf(\n"
		"	__global const ushort4* const leaf,\n"
		"	const uint idx)\n"
		"{\n"
		"	return (struct Leaf){\n"
		"		(ushort8)(\n"
		"			leaf[idx * 4 + 0],\n"
		"			leaf[idx * 4 + 1]\n"
		"		),\n"
		"		(ushort8)(\n"
		"			leaf[idx * 4 + 2],\n"
		"			leaf[idx * 4 + 3]\n"
		"		)\n"
		"	};\n"
		"}\n"
		"inline struct Voxel get_voxel(\n"
		"	__global const float4* const voxel,\n"
		"	const uint idx)\n"
		"{\n"
		"	return (struct Voxel){\n"
		"		voxel[idx * 2 + 0],\n"
		"		voxel[idx * 2 + 1]\n"
		"	};\n"
		"}\n"
		"\n"
		"uint traverself(\n"
		"	const struct Leaf leaf,\n"
		"	__global const float4* const voxel,\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray)\n"
		"{\n"
		"	struct ChildIndex child_index;\n"
		"\n"
		"	const uint hitCount = octlf_intersect_wide(\n"
		"		leaf,\n"
		"		bbox,\n"
		"		ray,\n"
		"		&child_index);\n"
		"\n"
		"	const uint8 leaf_start = convert_uint8(leaf.start);\n"
		"	const uint leafStart[8] = {\n"
		"		leaf_start.s0,\n"
		"		leaf_start.s1,\n"
		"		leaf_start.s2,\n"
		"		leaf_start.s3,\n"
		"		leaf_start.s4,\n"
		"		leaf_start.s5,\n"
		"		leaf_start.s6,\n"
		"		leaf_start.s7\n"
		"	};\n"
		"	const uint8 leaf_count = convert_uint8(leaf.count);\n"
		"	const uint leafCount[8] = {\n"
		"		leaf_count.s0,\n"
		"		leaf_count.s1,\n"
		"		leaf_count.s2,\n"
		"		leaf_count.s3,\n"
		"		leaf_count.s4,\n"
		"		leaf_count.s5,\n"
		"		leaf_count.s6,\n"
		"		leaf_count.s7\n"
		"	};\n"
		"	const float distance[8] = {\n"
		"		child_index.distance.s0,\n"
		"		child_index.distance.s1,\n"
		"		child_index.distance.s2,\n"
		"		child_index.distance.s3,\n"
		"		child_index.distance.s4,\n"
		"		child_index.distance.s5,\n"
		"		child_index.distance.s6,\n"
		"		child_index.distance.s7\n"
		"	};\n"
		"	const uint index[8] = {\n"
		"		child_index.index.s0,\n"
		"		child_index.index.s1,\n"
		"		child_index.index.s2,\n"
		"		child_index.index.s3,\n"
		"		child_index.index.s4,\n"
		"		child_index.index.s5,\n"
		"		child_index.index.s6,\n"
		"		child_index.index.s7\n"
		"	};\n"
		"	for (uint i = 0; i < hitCount; ++i) {\n"
		"		const uint payload_start = leafStart[index[i]];\n"
		"		const uint payload_count = leafCount[index[i]];\n"
		"		float nearest_dist = distance[i];\n"
		"		uint voxel_id = -1U;\n"
		"\n"
		"		for (uint j = payload_start; j < payload_start + payload_count; ++j) {\n"
		"			const struct Voxel payload = get_voxel(voxel, j);\n"
		"			const struct BBox payload_bbox = (struct BBox){ payload.min.xyz, payload.max.xyz };\n"
		"			const float dist = intersect(&payload_bbox, ray);\n"
		"\n"
		"			if (dist < nearest_dist) {\n"
		"				nearest_dist = dist;\n"
		"				voxel_id = as_uint(payload.min.w);\n"
		"			}\n"
		"		}\n"
		"\n"
		"		if (-1U != voxel_id)\n"
		"			return voxel_id;\n"
		"	}\n"
		"	return -1U;\n"
		"}\n"
		"\n"
		"uint traverse(\n"
		"	const struct Octet octet,\n"
		"	__global const ushort4* const leaf,\n"
		"	__global const float4* const voxel,\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray)\n"
		"{\n"
		"	struct ChildIndex child_index;\n"
		"	struct BBox child_bbox[8];\n"
		"\n"
		"	const uint hitCount = octet_intersect_wide(\n"
		"		octet,\n"
		"		bbox,\n"
		"		ray,\n"
		"		&child_index,\n"
		"		child_bbox);\n"
		"\n"
		"	const uint8 octet_child = convert_uint8(octet.child);\n"
		"	const uint octetChild[8] = {\n"
		"		octet_child.s0,\n"
		"		octet_child.s1,\n"
		"		octet_child.s2,\n"
		"		octet_child.s3,\n"
		"		octet_child.s4,\n"
		"		octet_child.s5,\n"
		"		octet_child.s6,\n"
		"		octet_child.s7\n"
		"	};\n"
		"	const uint index[8] = {\n"
		"		child_index.index.s0,\n"
		"		child_index.index.s1,\n"
		"		child_index.index.s2,\n"
		"		child_index.index.s3,\n"
		"		child_index.index.s4,\n"
		"		child_index.index.s5,\n"
		"		child_index.index.s6,\n"
		"		child_index.index.s7\n"
		"	};\n"
		"	for (uint i = 0; i < hitCount; ++i) {\n"
		"		const uint hit = traverself(get_leaf(leaf, octetChild[index[i]]), voxel, child_bbox + index[i], ray);\n"
		"		if (-1U != hit)\n"
		"			return hit;\n"
		"	}\n"
		"	return -1U;\n"
		"}\n"
		"\n"
		"__kernel __attribute__ ((vec_type_hint(float4)))\n"
		"void monokernel(\n"
		"	__global const ushort4* const src_a,\n"
		"	__global const ushort4* const src_b,\n"
		"	__global const float4* const src_c,\n"
#if OCL_QUIRK_0001
		"	__constant float* const src_d,\n"
#else
		"	__constant float4* const src_d,\n"
#endif
		"	__global uchar* const dst)\n"
		"{\n";

	const char source_image[] =
		"__constant sampler_t sampler_a = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n"
		"__constant sampler_t sampler_b = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n"
		"__constant sampler_t sampler_c = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n"
		"\n"
		"inline struct Octet get_octet(\n"
		"	__read_only image2d_t octet,\n"
		"	const uint idx)\n"
		"{\n"
		"	return (struct Octet){\n"
		"		(ushort8)(\n"
		"			convert_ushort4(read_imageui(octet, sampler_a, (int2)(0, idx))),\n"
		"			convert_ushort4(read_imageui(octet, sampler_a, (int2)(1, idx)))\n"
		"		)\n"
		"	};\n"
		"}\n"
		"inline struct Leaf get_leaf(\n"
		"	__read_only image2d_t leaf,\n"
		"	const uint idx)\n"
		"{\n"
		"	return (struct Leaf){\n"
		"		(ushort8)(\n"
		"			convert_ushort4(read_imageui(leaf, sampler_b, (int2)(0, idx))),\n"
		"			convert_ushort4(read_imageui(leaf, sampler_b, (int2)(1, idx)))\n"
		"		),\n"
		"		(ushort8)(\n"
		"			convert_ushort4(read_imageui(leaf, sampler_b, (int2)(2, idx))),\n"
		"			convert_ushort4(read_imageui(leaf, sampler_b, (int2)(3, idx)))\n"
		"		)\n"
		"	};\n"
		"}\n"
		"inline struct Voxel get_voxel(\n"
		"	__read_only image2d_t voxel,\n"
		"	const uint idx)\n"
		"{\n"
		"	return (struct Voxel){\n"
		"		read_imagef(voxel, sampler_c, (int2)(0, idx)),\n"
		"		read_imagef(voxel, sampler_c, (int2)(1, idx))\n"
		"	};\n"
		"}\n"
		"\n"
		"uint traverself(\n"
		"	const struct Leaf leaf,\n"
		"	__read_only image2d_t voxel,\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray)\n"
		"{\n"
		"	struct ChildIndex child_index;\n"
		"\n"
		"	const uint hitCount = octlf_intersect_wide(\n"
		"		leaf,\n"
		"		bbox,\n"
		"		ray,\n"
		"		&child_index);\n"
		"\n"
		"	const uint8 leaf_start = convert_uint8(leaf.start);\n"
		"	const uint leafStart[8] = {\n"
		"		leaf_start.s0,\n"
		"		leaf_start.s1,\n"
		"		leaf_start.s2,\n"
		"		leaf_start.s3,\n"
		"		leaf_start.s4,\n"
		"		leaf_start.s5,\n"
		"		leaf_start.s6,\n"
		"		leaf_start.s7\n"
		"	};\n"
		"	const uint8 leaf_count = convert_uint8(leaf.count);\n"
		"	const uint leafCount[8] = {\n"
		"		leaf_count.s0,\n"
		"		leaf_count.s1,\n"
		"		leaf_count.s2,\n"
		"		leaf_count.s3,\n"
		"		leaf_count.s4,\n"
		"		leaf_count.s5,\n"
		"		leaf_count.s6,\n"
		"		leaf_count.s7\n"
		"	};\n"
		"	const float distance[8] = {\n"
		"		child_index.distance.s0,\n"
		"		child_index.distance.s1,\n"
		"		child_index.distance.s2,\n"
		"		child_index.distance.s3,\n"
		"		child_index.distance.s4,\n"
		"		child_index.distance.s5,\n"
		"		child_index.distance.s6,\n"
		"		child_index.distance.s7\n"
		"	};\n"
		"	const uint index[8] = {\n"
		"		child_index.index.s0,\n"
		"		child_index.index.s1,\n"
		"		child_index.index.s2,\n"
		"		child_index.index.s3,\n"
		"		child_index.index.s4,\n"
		"		child_index.index.s5,\n"
		"		child_index.index.s6,\n"
		"		child_index.index.s7\n"
		"	};\n"
		"	for (uint i = 0; i < hitCount; ++i) {\n"
		"		const uint payload_start = leafStart[index[i]];\n"
		"		const uint payload_count = leafCount[index[i]];\n"
		"		float nearest_dist = distance[i];\n"
		"		uint voxel_id = -1U;\n"
		"\n"
		"		for (uint j = payload_start; j < payload_start + payload_count; ++j) {\n"
		"			const struct Voxel payload = get_voxel(voxel, j);\n"
		"			const struct BBox payload_bbox = (struct BBox){ payload.min.xyz, payload.max.xyz };\n"
		"			const float dist = intersect(&payload_bbox, ray);\n"
		"\n"
		"			if (dist < nearest_dist) {\n"
		"				nearest_dist = dist;\n"
		"				voxel_id = as_uint(payload.min.w);\n"
		"			}\n"
		"		}\n"
		"\n"
		"		if (-1U != voxel_id)\n"
		"			return voxel_id;\n"
		"	}\n"
		"	return -1U;\n"
		"}\n"
		"\n"
		"uint traverse(\n"
		"	const struct Octet octet,\n"
		"	__read_only image2d_t leaf,\n"
		"	__read_only image2d_t voxel,\n"
		"	const struct BBox* const bbox,\n"
		"	const struct Ray* const ray)\n"
		"{\n"
		"	struct ChildIndex child_index;\n"
		"	struct BBox child_bbox[8];\n"
		"\n"
		"	const uint hitCount = octet_intersect_wide(\n"
		"		octet,\n"
		"		bbox,\n"
		"		ray,\n"
		"		&child_index,\n"
		"		child_bbox);\n"
		"\n"
		"	const uint8 octet_child = convert_uint8(octet.child);\n"
		"	const uint octetChild[8] = {\n"
		"		octet_child.s0,\n"
		"		octet_child.s1,\n"
		"		octet_child.s2,\n"
		"		octet_child.s3,\n"
		"		octet_child.s4,\n"
		"		octet_child.s5,\n"
		"		octet_child.s6,\n"
		"		octet_child.s7\n"
		"	};\n"
		"	const uint index[8] = {\n"
		"		child_index.index.s0,\n"
		"		child_index.index.s1,\n"
		"		child_index.index.s2,\n"
		"		child_index.index.s3,\n"
		"		child_index.index.s4,\n"
		"		child_index.index.s5,\n"
		"		child_index.index.s6,\n"
		"		child_index.index.s7\n"
		"	};\n"
		"	for (uint i = 0; i < hitCount; ++i) {\n"
		"		const uint hit = traverself(get_leaf(leaf, octetChild[index[i]]), voxel, child_bbox + index[i], ray);\n"
		"		if (-1U != hit)\n"
		"			return hit;\n"
		"	}\n"
		"	return -1U;\n"
		"}\n"
		"\n"
		"__kernel __attribute__ ((vec_type_hint(float4)))\n"
		"void monokernel(\n"
		"	__read_only image2d_t src_a,\n"
		"	__read_only image2d_t src_b,\n"
		"	__read_only image2d_t src_c,\n"
#if OCL_QUIRK_0001
		"	__constant float* const src_d,\n"
#else
		"	__constant float4* const src_d,\n"
#endif
		"	__global uchar* const dst)\n"
		"{\n";

	const char source_buffer_output[] =
		"	dst[get_global_size(0) * idy + idx] = (uchar)result;\n"
		"}";

	const char source_implicit_mad[] =
		"	const int idx = (int)get_global_id(0);\n"
		"	const int idy = (int)get_global_id(1);\n"
		"	const int dimx = (int)get_global_size(0);\n"
		"	const int dimy = (int)get_global_size(1);\n"
#if OCL_QUIRK_0001
		"	const float3 cam0 = (float3)(src_d[0], src_d[1], src_d[ 2]);\n"
		"	const float3 cam1 = (float3)(src_d[4], src_d[5], src_d[ 6]);\n"
		"	const float3 cam2 = (float3)(src_d[8], src_d[9], src_d[10]);\n"
		"	const float3 ray_origin = (float3)(src_d[12], src_d[13], src_d[14]);\n"
		"	const float3 bbox_min   = (float3)(src_d[16], src_d[17], src_d[18]);\n"
		"	const float3 bbox_max   = (float3)(src_d[20], src_d[21], src_d[22]);\n"
#else
		"	const float3 cam0 = src_d[0].xyz;\n"
		"	const float3 cam1 = src_d[1].xyz;\n"
		"	const float3 cam2 = src_d[2].xyz;\n"
		"	const float3 ray_origin = src_d[3].xyz;\n"
		"	const float3 bbox_min   = src_d[4].xyz;\n"
		"	const float3 bbox_max   = src_d[5].xyz;\n"
#endif
		"	const float3 ray_direction =\n"
		"		cam0 * ((idx * 2 - dimx) * (1.0f / dimx)) +\n"
		"		cam1 * ((idy * 2 - dimy) * (1.0f / dimy)) +\n"
		"		cam2;\n"
		"	const float3 ray_rcpdir = fmax(fmin(1.0f / ray_direction, MAXFLOAT), -MAXFLOAT);\n"
		"	const struct Ray ray = (struct Ray){ ray_origin, ray_rcpdir };\n"
		"	const struct BBox root_bbox = (struct BBox){ bbox_min, bbox_max };\n"
		"	uint result = traverse(get_octet(src_a, 0), src_b, src_c, &root_bbox, &ray);\n"
		"	result = -1U != result ? 16 + result * 8: 0;\n";

	const char* source[] = {
		source_prologue,
		source_buffer,
		source_implicit_mad,
		source_buffer_output
	};

	if (PARAM_TYPE_IMAGE == kern_param_type)
		source[1] = source_image;

	cl_program program = clCreateProgramWithSource(context,
		COUNT_OF(source), source, 0, &success);

	if (reportCLError(success)) {
		stream::cerr << "error creating program with source\n";
		return -1;
	}

	const scoped_ptr< cl_program, scoped_functor > release_program(&program);

	success = clBuildProgram(program, 1, device() + device_idx, "-cl-mad-enable", 0, 0);

	if (reportCLError(success) && CL_BUILD_PROGRAM_FAILURE != success) {
		stream::cerr << "error building program\n";
		return -1;
	}

	cl_build_status build_status = CL_BUILD_NONE;
	success = clGetProgramBuildInfo(program,
		device()[device_idx], CL_PROGRAM_BUILD_STATUS, sizeof(build_status), &build_status, 0);

	if (reportCLError(success)) {
		stream::cerr << "error getting build info (build_status)\n";
		return -1;
	}

	if (CL_BUILD_SUCCESS != build_status || OCL_KERNEL_BUILD_VERBOSE) {
		size_t log_len = 0;
		success = clGetProgramBuildInfo(program,
			device()[device_idx], CL_PROGRAM_BUILD_LOG, 0, 0, &log_len);

		if (reportCLError(success)) {
			stream::cerr << "error getting build info (log_len)\n";
			return -1;
		}

		const scoped_ptr< char, generic_free > build_log(
			reinterpret_cast< char* >(std::calloc(log_len, sizeof(char))));

		success = clGetProgramBuildInfo(program,
			device()[device_idx], CL_PROGRAM_BUILD_LOG, log_len, build_log(), 0);

		if (reportCLError(success)) {
			stream::cerr << "error getting build info (build_log)\n";
			return -1;
		}

		stream::cerr << build_log() << '\n';

		if (CL_BUILD_SUCCESS != build_status)
			return -1;
	}

	cl_kernel kernel = clCreateKernel(program, "monokernel", &success);

	if (reportCLError(success)) {
		stream::cerr << "error creating kernel\n";
		return -1;
	}

	scoped_ptr< cl_kernel, scoped_functor > release_kernel(&kernel);

	size_t local_ws_multiple = 0;
	success = clGetKernelWorkGroupInfo(kernel,
		device()[device_idx], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
		sizeof(local_ws_multiple), &local_ws_multiple, 0);

	if (reportCLError(success)) {
		stream::cerr << "error getting kernel preferred workgroup size multiple info\n";
		return -1;
	}

	stream::cout << "kernel preferred workgroup size multiple: " << local_ws_multiple << '\n';

	cl_uint max_dim = 0;
	success = clGetDeviceInfo(
		device()[device_idx], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
		sizeof(max_dim), &max_dim, 0);

	if (reportCLError(success)) {
		stream::cerr << "error getting device max work-item dimensions info\n";
		return -1;
	}

	size_t buffer_devinfo[8];

	if (max_dim > sizeof(buffer_devinfo) / sizeof(buffer_devinfo[0])) {
		stream::cerr << "device max work-item dimensions exceed expectations; bailing out\n";
		return -1;
	}

	success = clGetDeviceInfo(
		device()[device_idx], CL_DEVICE_MAX_WORK_ITEM_SIZES,
		sizeof(buffer_devinfo), buffer_devinfo, 0);

	if (reportCLError(success)) {
		stream::cerr << "error getting device max work-item sizes info\n";
		return -1;
	}

	stream::cout << "device max work-item sizes:";
	for (cl_uint i = 0; i < max_dim; ++i)
		stream::cout << ' ' << buffer_devinfo[i];
	stream::cout << '\n';

	const size_t work_item_size_0 = buffer_devinfo[0];
	const size_t work_item_size_1 = buffer_devinfo[1];

	success = clGetDeviceInfo(
		device()[device_idx], CL_DEVICE_MAX_WORK_GROUP_SIZE,
		sizeof(buffer_devinfo), buffer_devinfo, 0);

	if (reportCLError(success)) {
		stream::cerr << "error getting device max work-group size info\n";
		return -1;
	}

	const size_t max_local_ws = buffer_devinfo[0];
	const size_t latency_hiding_factor = local_ws_multiple * 2 > max_local_ws ? 1 : 2;
	const size_t combined_item_size_0 = local_ws_multiple * latency_hiding_factor;

	global_ws[0] = image_w;
	global_ws[1] = image_h;

	if (latency_hiding_factor > work_item_size_1) {
		if (combined_item_size_0 > work_item_size_0) {
			local_ws[0] = work_item_size_0;
			local_ws[1] = 1;
		}
		else {
			local_ws[0] = combined_item_size_0;
			local_ws[1] = 1;
		}
	}
	else {
		local_ws[0] = local_ws_multiple;
		local_ws[1] = latency_hiding_factor;
	}

	for (cl_uint i = 0; i < work_dim; ++i)
		if (0 != global_ws[i] % local_ws[i]) {
			stream::cerr << "ND range not a multiple of workgroup size; bailing out\n";
			return 1;
		}

	// prepare the playfield ///////////////////////////////////////////////////
	timeline.setCapacity(scene_count);
	timeline.addMultiElement(scene_count);

	// set initial external storage to the octree of the 1st scene
	// note: practically all buffers of the octree require 16-byte alignment; since this is
	// guaranteed by 64-bit malloc, we don't do anything WRT alignment here /32-bit caveat
	// note: all we need from the tree of each scene prior to frame loop is the root bbox;
	// let the trees of all scenes overwrite each other in the same buffer /practical cheat
	timeline.getMutable(scene_1).set_extrnal_storage(
		octet_count, octet_map_buffer[0],
		leaf_count, leaf_map_buffer[0],
		voxel_count, voxel_map_buffer[0]);

	// set initial external storage to the octree of the 2nd scene (same as storage for 1st and 3rd scenes)
	timeline.getMutable(scene_2).set_extrnal_storage(
		octet_count, octet_map_buffer[0],
		leaf_count, leaf_map_buffer[0],
		voxel_count, voxel_map_buffer[0]);

	// set initial external storage to the octree of the 3rd scene (same as storage for 1st and 2nd scenes)
	timeline.getMutable(scene_3).set_extrnal_storage(
		octet_count, octet_map_buffer[0],
		leaf_count, leaf_map_buffer[0],
		voxel_count, voxel_map_buffer[0]);

	if (!scene1.init(timeline.getMutable(scene_1)))
		return 1;

	if (!scene2.init(timeline.getMutable(scene_2)))
		return 2;

	if (!scene3.init(timeline.getMutable(scene_3)))
		return 3;

	track_cursor = 0;
	action_count = 0;

	// use first scene's initial world bbox to compute a normalization (pan_n_zoom) matrix
	const BBox& world_bbox = timeline.getElement(scene_1).get_root_bbox();

	bbox_min = world_bbox.get_min();
	bbox_max = world_bbox.get_max();
	centre = (bbox_max + bbox_min) * simd::f32x4(.5f);
	extent = (bbox_max - bbox_min) * simd::f32x4(.5f);
	max_extent = std::max(extent[0], std::max(extent[1], extent[2]));

	// frame initialization finished successfully -- dismiss any deinitters and persist non-locals
	octet_map.reset();
	leaf_map.reset();
	voxel_map.reset();
	carb_map.reset();
	image_map.reset();

	persist_queue = queue;
	persist_kernel = kernel;

	release_src_a.reset();
	release_src_b.reset();
	release_src_c.reset();
	release_src_d.reset();
	release_dst.reset();
	release_kernel.reset();
	release_queue.reset();
	release_monv.reset();
	return 0;
}

int deinitFrame(void)
{
	for (size_t i = 0; i < n_buffering; ++i) {
		clReleaseMemObject(src_a_d[i]);
		clReleaseMemObject(src_b_d[i]);
		clReleaseMemObject(src_c_d[i]);
		clReleaseMemObject(src_d_d[i]);
		clReleaseMemObject(dst_d[i]);
	}
	clReleaseKernel(persist_kernel);
	clReleaseCommandQueue(persist_queue);
	testbed::monv::deinit_resources();
	return 0;
}

int renderFrame(void)
{
	using clutil::reportCLError;

	const bool report_kernel_time           = param.flags & cli_param::BIT_REPORT_KERNEL_TIME;
	const KERNEL_PARAM_TYPE kern_param_type = param.kern_param_type;
	const size_t image_w                    = param.image_w;
	const size_t image_h                    = param.image_h;
	const size_t frames                     = param.frames;

	const size_t mem_size_image = image_w * image_h * sizeof(cl_uchar);
	const size_t frame = frame_idx++;
	const cl_command_queue queue = persist_queue;
	const cl_kernel kernel = persist_kernel;

	// have we produced enough frames?
	if (frame >= frames)
		return 1;

	const compile_assert< 2 == n_buffering > assert_double_buffering;

#if FRAME_RATE == 0
	const uint64_t tframe = timer_ns();

	if (0 == frame)
		tlast = tframe;

	const float dt = double(tframe - tlast) * 1e-9;
	tlast = tframe;

#else
	const float dt = 1.0 / FRAME_RATE;

#endif
	// upate run time (we aren't supposed to run long - fp32 should do) and beat time
	c::accum_time += dt;
	c::accum_beat   = wrap_at_period(c::accum_beat   + dt, c::beat_period);
	c::accum_beat_2 = wrap_at_period(c::accum_beat_2 + dt, c::beat_period * 2.0);

	// run all live actions, retiring the completed ones
	for (size_t i = 0; i < action_count; ++i)
		if (!action[i]->frame(dt))
			action[i--] = action[--action_count];

	// start any pending actions
	for (; track_cursor < sizeof(track) / sizeof(track[0]) && c::accum_time >= track[track_cursor].start; ++track_cursor)
		if (track[track_cursor].action.start(c::accum_time - track[track_cursor].start, track[track_cursor].duration)) {
			if (action_count == sizeof(action) / sizeof(action[0])) {
				stream::cerr << "error: too many pending actions\n";
				return 999;
			}

			action[action_count++] = &track[track_cursor].action;
		}

	// set proper external storage to the octree of the live scene
	// note: practically all buffers of the octree require 16-byte alignment; since this is
	// guaranteed by 64-bit malloc, we don't do anything WRT alignment here /32-bit caveat
	timeline.getMutable(c::scene_selector).set_extrnal_storage(
		octet_count, octet_map_buffer[frame & 1],
		leaf_count, leaf_map_buffer[frame & 1],
		voxel_count, voxel_map_buffer[frame & 1]);

	// run the live scene
	if (!scene[c::scene_selector]->frame(timeline.getMutable(c::scene_selector), dt))
		stream::cerr << "failure building frame " << frame << '\n';

	// produce camera for the new frame;
	// collapse S * T and T * S operators as follows:
	//
	//	s	0	0	0		1	0	0	0		s	0	0	0
	//	0	s	0	0	*	0	1	0	0	=	0	s	0	0
	//	0	0	s	0		0	0	1	0		0	0	s	0
	//	0	0	0	1		x	y	z	1		x	y	z	1
	//
	//	1	0	0	0		s	0	0	0		s	0	0	0
	//	0	1	0	0	*	0	s	0	0	=	0	s	0	0
	//	0	0	1	0		0	0	s	0		0	0	s	0
	//	x	y	z	1		0	0	0	1		sx	sy	sz	1

	// forward: pan * zoom * rot * eyep
	// inverse: (eyep)-1 * rotT * (zoom)-1 * (pan)-1

	const matx4 rot =
		matx4_rotate(scene[c::scene_selector]->get_roll(),           0.f, 1.f, 0.f) *
		matx4_rotate(scene[c::scene_selector]->get_azim() + c::azim, 0.f, 0.f, 1.f) *
		matx4_rotate(scene[c::scene_selector]->get_decl() + c::decl, 1.f, 0.f, 0.f);

	const matx4 eyep(
		1.f, 0.f, 0.f, 0.f,
		0.f, 1.f, 0.f, 0.f,
		0.f, 0.f, 1.f, 0.f,
		scene[c::scene_selector]->get_cam_x() + c::pos_x,
		scene[c::scene_selector]->get_cam_y() + c::pos_y,
		scene[c::scene_selector]->get_cam_z() + c::pos_z, 1.f);

	const matx4 zoom_n_pan(
		max_extent, 0.f, 0.f, 0.f,
		0.f, max_extent, 0.f, 0.f,
		0.f, 0.f, max_extent, 0.f,
		centre[0] - scene[c::scene_selector]->get_offset_x(),
		centre[1] - scene[c::scene_selector]->get_offset_y(),
		centre[2] - scene[c::scene_selector]->get_offset_z(), 1.f);

	const matx4 mv_inv = eyep * transpose(rot) * zoom_n_pan;

	vect3 (& carb)[carb_count] = *reinterpret_cast< vect3 (*)[carb_count] >(carb_map_buffer[frame & 1]);
	// camera
	carb[0] = vect3(mv_inv[0][0], mv_inv[0][1], mv_inv[0][2]);
	carb[1] = vect3(mv_inv[1][0], mv_inv[1][1], mv_inv[1][2]) * vect3(float(image_h) / image_w);
	carb[2] = vect3(mv_inv[2][0], mv_inv[2][1], mv_inv[2][2]) * vect3(-1);
	carb[3] = vect3(mv_inv[3][0], mv_inv[3][1], mv_inv[3][2]);
	// root bbox
	carb[4] = timeline.getElement(c::scene_selector).get_root_bbox().get_min();
	carb[5] = timeline.getElement(c::scene_selector).get_root_bbox().get_max();

	cl_int success;

	// fill-in new frame data set
	if (PARAM_TYPE_IMAGE == kern_param_type) {
		const size_t origin[] = { 0, 0, 0 };

		const size_t region_src_a[] = { octet_w, octet_h, 1 };
		success = clEnqueueWriteImage(queue, src_a_d[frame & 1], CL_FALSE, origin, region_src_a, sizeof(*octet_map_buffer[0]), 0, octet_map_buffer[frame & 1],
			0, 0, &event_data_set_ready[frame & 1][0]);

		if (reportCLError(success)) {
			stream::cerr << "error enqueuing image write for image octet\n";
			return -1;
		}

		const size_t region_src_b[] = { leaf_w, leaf_h, 1 };
		success = clEnqueueWriteImage(queue, src_b_d[frame & 1], CL_FALSE, origin, region_src_b, sizeof(*leaf_map_buffer[0]), 0, leaf_map_buffer[frame & 1],
			0, 0, &event_data_set_ready[frame & 1][1]);

		if (reportCLError(success)) {
			stream::cerr << "error enqueuing image write for image leaf\n";
			return -1;
		}

		const size_t region_src_c[] = { voxel_w, voxel_h, 1 };
		success = clEnqueueWriteImage(queue, src_c_d[frame & 1], CL_FALSE, origin, region_src_c, sizeof(*voxel_map_buffer[0]), 0, voxel_map_buffer[frame & 1],
			0, 0, &event_data_set_ready[frame & 1][2]);

		if (reportCLError(success)) {
			stream::cerr << "error enqueuing image write for image voxel\n";
			return -1;
		}
	}
	else {
		success = clEnqueueWriteBuffer(queue, src_a_d[frame & 1], CL_FALSE, 0, mem_size_octet, octet_map_buffer[frame & 1],
			0, 0, &event_data_set_ready[frame & 1][0]);

		if (reportCLError(success)) {
			stream::cerr << "error enqueuing buffer write for buffer octet\n";
			return -1;
		}

		success = clEnqueueWriteBuffer(queue, src_b_d[frame & 1], CL_FALSE, 0, mem_size_leaf, leaf_map_buffer[frame & 1],
			0, 0, &event_data_set_ready[frame & 1][1]);

		if (reportCLError(success)) {
			stream::cerr << "error enqueuing buffer write for buffer leaf\n";
			return -1;
		}

		success = clEnqueueWriteBuffer(queue, src_c_d[frame & 1], CL_FALSE, 0, mem_size_voxel, voxel_map_buffer[frame & 1],
			0, 0, &event_data_set_ready[frame & 1][2]);

		if (reportCLError(success)) {
			stream::cerr << "error enqueuing buffer write for buffer voxel\n";
			return -1;
		}
	}

	success = clEnqueueWriteBuffer(queue, src_d_d[frame & 1], CL_FALSE, 0, mem_size_carb, carb_map_buffer[frame & 1],
		0, 0, &event_data_set_ready[frame & 1][3]);

	if (reportCLError(success)) {
		stream::cerr << "error enqueuing buffer write for buffer cam\n";
		return -1;
	}

	// show result from one frame ago
	if (0 != frame) {
		const size_t ready_frame = frame + 1;

		success = clEnqueueReadBuffer(queue, dst_d[ready_frame & 1], CL_TRUE, 0, mem_size_image, image_map_buffer[ready_frame & 1],
			1, &event_kernel_complete[ready_frame & 1], 0);

		if (reportCLError(success)) {
			stream::cerr << "error waiting for ocl done\n";
			return -1;
		}

		testbed::monv::render(image_map_buffer[ready_frame & 1]);

		// profile the frame
		if (report_kernel_time) {
			cl_ulong t0, t1;

			success = clGetEventProfilingInfo(
				event_kernel_complete[ready_frame & 1], CL_PROFILING_COMMAND_START, sizeof(t0), &t0, 0);

			if (reportCLError(success)) {
				stream::cerr << "error getting profiling info (command_start)\n";
				return -1;
			}

			success = clGetEventProfilingInfo(
				event_kernel_complete[ready_frame & 1], CL_PROFILING_COMMAND_END, sizeof(t1), &t1, 0);

			if (reportCLError(success)) {
				stream::cerr << "error getting profiling info (command_end)\n";
				return -1;
			}

			stream::cout << "elapsed time (ns): " << t1 - t0 << '\n';
		}
	}

	// launch kernel on new frame data set
	success = clSetKernelArg(kernel, 0, sizeof(src_a_d[0]), &src_a_d[frame & 1]);

	if (reportCLError(success)) {
		stream::cerr << "error setting kernel arg 0\n";
		return -1;
	}

	success = clSetKernelArg(kernel, 1, sizeof(src_b_d[0]), &src_b_d[frame & 1]);

	if (reportCLError(success)) {
		stream::cerr << "error setting kernel arg 1\n";
		return -1;
	}

	success = clSetKernelArg(kernel, 2, sizeof(src_c_d[0]), &src_c_d[frame & 1]);

	if (reportCLError(success)) {
		stream::cerr << "error setting kernel arg 2\n";
		return -1;
	}

	success = clSetKernelArg(kernel, 3, sizeof(src_d_d[0]), &src_d_d[frame & 1]);

	if (reportCLError(success)) {
		stream::cerr << "error setting kernel arg 3\n";
		return -1;
	}

	success = clSetKernelArg(kernel, 4, sizeof(dst_d[0]), &dst_d[frame & 1]);

	if (reportCLError(success)) {
		stream::cerr << "error setting kernel arg 4\n";
		return -1;
	}

	success = clEnqueueNDRangeKernel(queue, kernel, work_dim, 0, global_ws, local_ws,
		sizeof(event_data_set_ready[0]) / sizeof(event_data_set_ready[0][0]), event_data_set_ready[frame & 1], &event_kernel_complete[frame & 1]);

	if (reportCLError(success)) {
		stream::cerr << "error enqueuing kernel\n";
		return -1;
	}

	return 0;
}
