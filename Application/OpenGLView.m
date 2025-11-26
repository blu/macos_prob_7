/*
 Copyright (C) 2015 Apple Inc. All Rights Reserved.
 See LICENSE.txt for this sample’s licensing information
 
 Abstract:
 OpenGL view subclass.
 */

#import "OpenGLView.h"
#import "OpenGLRenderer.h"
#import "param.h"

#if CGL_VERSION_1_3 == 0
#error missing CGL_VERSION_1_3
#endif

@interface OpenGLView ()
{
	CVDisplayLinkRef _displayLink;

	OpenGLRenderer* _renderer;
}
@end

@implementation OpenGLView

- (CVReturn) getFrameForTime:(const CVTimeStamp*)outputTime
{
	// There is no autorelease pool when this method is called
	// because it will be called from a background thread.
	// It's important to create one or app can leak objects.
	@autoreleasepool {
		[self drawView];
	}
	return kCVReturnSuccess;
}

// This is the renderer output callback function
static CVReturn MyDisplayLinkCallback(CVDisplayLinkRef displayLink,
									  const CVTimeStamp* now,
									  const CVTimeStamp* outputTime,
									  CVOptionFlags flagsIn,
									  CVOptionFlags* flagsOut, 
									  void* displayLinkContext)
{
    CVReturn result = [(__bridge OpenGLView*)displayLinkContext getFrameForTime:outputTime];
    return result;
}

- (instancetype) init
{
	self = [super init];
	if (self) {
		NSOpenGLPixelFormatAttribute attrs[256] = {
			NSOpenGLPFADoubleBuffer,
			NSOpenGLPFADepthSize, 24,
			// Must specify the 3.2 Core Profile to use OpenGL 3.2
			NSOpenGLPFAOpenGLProfile,
			NSOpenGLProfileVersion3_2Core,
			0
		};

		const size_t attrsMax = sizeof(attrs) / sizeof(attrs[0]);
		size_t extraAttrs = 0;
		while (extraAttrs < attrsMax && attrs[extraAttrs])
			++extraAttrs;

		if (param.fsaa) {
			if (extraAttrs + 5 < attrsMax) {
				attrs[extraAttrs++] = NSOpenGLPFAMultisample;
				attrs[extraAttrs++] = NSOpenGLPFASampleBuffers;
				attrs[extraAttrs++] = (NSOpenGLPixelFormatAttribute) 1;
				attrs[extraAttrs++] = NSOpenGLPFASamples;
				attrs[extraAttrs++] = (NSOpenGLPixelFormatAttribute) param.fsaa;
			}
			else {
				NSLog(@"error: NSOpenGLPixelFormatAttribute array capacity exceeded");
			}
		}

		const unsigned pixelBitness = param.bitness[0] + param.bitness[1] + param.bitness[2] + param.bitness[3];
		if (pixelBitness) {
			if (extraAttrs + 2 < attrsMax) {
				attrs[extraAttrs++] = NSOpenGLPFAColorSize;
				attrs[extraAttrs++] = (NSOpenGLPixelFormatAttribute) pixelBitness;
			}
			else {
				NSLog(@"error: NSOpenGLPixelFormatAttribute array capacity exceeded");
			}
		}

		// terminate the attrib array
		attrs[extraAttrs] = 0;

		NSOpenGLPixelFormat *pf = [[NSOpenGLPixelFormat alloc] initWithAttributes:attrs];

		if (!pf) {
			NSLog(@"error: no OpenGL pixel format");
		}

		GLint attrFeedback;
		[pf getValues:&attrFeedback forAttribute:NSOpenGLPFAColorSize forVirtualScreen:0];

		NSLog(@"OpenGL context of NSOpenGLPFAColorSize %u", (uint32_t) attrFeedback);

		NSOpenGLContext* context = [[NSOpenGLContext alloc] initWithFormat:pf shareContext:nil];

#if defined(DEBUG)
		// When we're using a CoreProfile context, crash if we call a legacy OpenGL function
		// This will make it much more obvious where and when such a function call is made so
		// that we can remove such calls.
		// Without this we'd simply get GL_INVALID_OPERATION error for calling legacy functions
		// but it would be more difficult to see where that function was called.
		CGLEnable([context CGLContextObj], kCGLCECrashOnRemovedFunctions);
#endif

		[self setPixelFormat:pf];

		[self setOpenGLContext:context];

		// Opt-In to Retina resolution
		[self setWantsBestResolutionOpenGLSurface:YES];
	}
	return self;
}

- (void) prepareOpenGL
{
	[super prepareOpenGL];

	// Make all the OpenGL calls to setup rendering  
	//  and build the necessary rendering objects
	[self initGL];

	// Create a display link capable of being used with all active displays
	CVDisplayLinkCreateWithActiveCGDisplays(&_displayLink);

	// Set the renderer output callback function
	CVDisplayLinkSetOutputCallback(_displayLink, &MyDisplayLinkCallback, (__bridge void*)self);

	// Set the display link for the current renderer
	CGLContextObj cglContext = [[self openGLContext] CGLContextObj];
	CGLPixelFormatObj cglPixelFormat = [[self pixelFormat] CGLPixelFormatObj];
	CVDisplayLinkSetCurrentCGDisplayFromOpenGLContext(_displayLink, cglContext, cglPixelFormat);

	// Activate the display link
	CVDisplayLinkStart(_displayLink);
}

- (void) initGL
{
	// The reshape function may have changed the thread to which our OpenGL
	// context is attached before prepareOpenGL and initGL are called.  So call
	// makeCurrentContext to ensure that our OpenGL context current to this 
	// thread (i.e. makeCurrentContext directs all OpenGL calls on this thread
	// to [self openGLContext])
	[[self openGLContext] makeCurrentContext];

	// Synchronize buffer swaps with vertical refresh rate
	GLint swapInt = 1;
	[[self openGLContext] setValues:&swapInt forParameter:NSOpenGLCPSwapInterval];

	// Init our renderer.  Use 0 for the defaultFBO which is appropriate for
	// OSX (but not iOS since iOS apps must create their own FBO)
	_renderer = [[OpenGLRenderer alloc] initWithDefaultFBO:0];
}

- (void)reshape
{
	[super reshape];

	// We draw on a secondary thread through the display link. However, when
	// resizing the view, -drawRect is called on the main thread.
	// Add a mutex around to avoid the threads accessing the context
	// simultaneously when resizing.
	CGLLockContext([[self openGLContext] CGLContextObj]);

	// Get the view size in Points
	NSRect viewRectPoints = [self bounds];

	// Rendering at retina resolutions will reduce aliasing, but at the potential
	// cost of framerate and battery life due to the GPU needing to render more
	// pixels.

	// Any calculations the renderer does which use pixel dimentions, must be
	// in "retina" space.  [NSView convertRectToBacking] converts point sizes
	// to pixel sizes.  Thus the renderer gets the size in pixels, not points,
	// so that it can set it's viewport and perform and other pixel based
	// calculations appropriately.
	// viewRectPixels will be larger than viewRectPoints for retina displays.
	// viewRectPixels will be the same as viewRectPoints for non-retina displays
	NSRect viewRectPixels = [self convertRectToBacking:viewRectPoints];

	if (viewRectPixels.size.width > 0 && viewRectPixels.size.height > 0) {
		// Set the new dimensions in our renderer
		[_renderer resizeWithWidth:viewRectPixels.size.width
						 AndHeight:viewRectPixels.size.height];
	}

	CGLUnlockContext([[self openGLContext] CGLContextObj]);
}

- (void)renewGState
{
	// Called whenever graphics state updated (such as window resize)
	
	// OpenGL rendering is not synchronous with other rendering on the OSX.
	// Therefore, call disableScreenUpdatesUntilFlush so the window server
	// doesn't render non-OpenGL content in the window asynchronously from
	// OpenGL content, which could cause flickering.  (non-OpenGL content
	// includes the title bar and drawing done by the app with other APIs)
	[[self window] disableScreenUpdatesUntilFlush];

	[super renewGState];
}

- (void) drawRect: (NSRect) theRect
{
	// Called during resize operations
	
	// Avoid flickering during resize by drawiing	
	[self drawView];
}

- (void) drawView
{
	[[self openGLContext] makeCurrentContext];

	// We draw on a secondary thread through the display link
	// When resizing the view, -reshape is called automatically on the main
	// thread. Add a mutex around to avoid the threads accessing the context
	// simultaneously when resizing
	CGLLockContext([[self openGLContext] CGLContextObj]);

	[_renderer render];

	CGLFlushDrawable([[self openGLContext] CGLContextObj]);
	CGLUnlockContext([[self openGLContext] CGLContextObj]);
}

- (void) deinit
{
	// Stop the display link BEFORE releasing anything in the view
	// otherwise the display link thread may call into the view and crash
	// when it encounters something that has been release
	CVDisplayLinkStop(_displayLink);

	CVDisplayLinkRelease(_displayLink);

	_renderer = nil;
}

@end
