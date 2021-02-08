/*
 Copyright (C) 2015 Apple Inc. All Rights Reserved.
 See LICENSE.txt for this sample’s licensing information
 
 Abstract:
 The OpenGLRenderer class creates and draws objects.
  Most of the code is OS independent.
 */

#import <AppKit/NSApplication.h>
#import "OpenGLRenderer.h"
#import "../param.h"

#ifndef NULL
#define NULL 0
#endif

#define BUFFER_OFFSET(i) ((char *)NULL + (i))

static inline const char * GetGLErrorString(GLenum error)
{
	const char *str;
	switch( error ) {
		case GL_NO_ERROR:
			str = "GL_NO_ERROR";
			break;
		case GL_INVALID_ENUM:
			str = "GL_INVALID_ENUM";
			break;
		case GL_INVALID_VALUE:
			str = "GL_INVALID_VALUE";
			break;
		case GL_INVALID_OPERATION:
			str = "GL_INVALID_OPERATION";
			break;
#if defined __gl_h_ || defined __gl3_h_
		case GL_OUT_OF_MEMORY:
			str = "GL_OUT_OF_MEMORY";
			break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			str = "GL_INVALID_FRAMEBUFFER_OPERATION";
			break;
#endif
#if defined __gl_h_
		case GL_STACK_OVERFLOW:
			str = "GL_STACK_OVERFLOW";
			break;
		case GL_STACK_UNDERFLOW:
			str = "GL_STACK_UNDERFLOW";
			break;
		case GL_TABLE_TOO_LARGE:
			str = "GL_TABLE_TOO_LARGE";
			break;
#endif
		default:
			str = "(ERROR: Unknown Error Enum)";
			break;
	}
	return str;
}

#define GetGLError()                                    \
{                                                       \
	GLenum err = glGetError();                          \
	while (err != GL_NO_ERROR) {                        \
		NSLog(@"GLError %s set in File:%s Line:%d\n",   \
		GetGLErrorString(err), __FILE__, __LINE__);     \
		err = glGetError();                             \
	}                                                   \
}

@interface OpenGLRenderer ()
{
    GLuint _defaultFBOName;
}
@end

@implementation OpenGLRenderer


- (void) resizeWithWidth:(GLuint)width AndHeight:(GLuint)height
{
	NSLog(@"OpenGLRenderer resize %u %u", width, height);
	glViewport(0, 0, width, height);
}


- (void) render
{
	if (renderFrame()) {
		[[NSApplication sharedApplication] terminate:nil];
	}
}


- (id) initWithDefaultFBO: (GLuint) defaultFBOName
{
	if((self = [super init]))
	{
		NSLog(@"%s %s", glGetString(GL_RENDERER), glGetString(GL_VERSION));

		_defaultFBOName = defaultFBOName;

		if (initFrame())
		{
			[[NSApplication sharedApplication] terminate:nil];
		}

		// Check for errors to make sure all of our setup went ok
		GetGLError();
	}

	return self;
}


- (void) dealloc
{
	deinitFrame();
}

@end
