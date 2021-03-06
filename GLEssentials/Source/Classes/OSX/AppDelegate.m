/*
 Copyright (C) 2015 Apple Inc. All Rights Reserved.
 See LICENSE.txt for this sample’s licensing information
 
 Abstract:
 The Application Delegate.
 */

#import "AppDelegate.h"
#import "GLEssentialsWindowController.h"
#import "GLEssentialsGLView.h"
#import "../../param.h"

@interface AppDelegate ()
{
	GLEssentialsWindowController *_controller;
}
@end

@implementation AppDelegate

- (instancetype) init {
	if (self = [super init]) {
#if SUPPORT_RETINA_RESOLUTION
		const CGFloat res_w = param.image_w / 2;
		const CGFloat res_h = param.image_h / 2;
#else
		const CGFloat res_w = param.image_w;
		const CGFloat res_h = param.image_h;
#endif
		GLEssentialsGLView *view = [[GLEssentialsGLView alloc] init];
		NSWindow *window = [[NSWindow alloc] initWithContentRect: NSMakeRect(0, 0, res_w, res_h)
		                                               styleMask: NSWindowStyleMaskTitled
		                                                 backing: NSBackingStoreBuffered
		                                                   defer: NO];
		[window setContentView: view];
		_controller = [[GLEssentialsWindowController alloc] initWithWindow: window];
	}
	return self;
}

- (void)applicationWillFinishLaunching:(NSNotification *)aNotification {
	NSWindow *window = [_controller window];
	window.title = NSProcessInfo.processInfo.processName;
	[window cascadeTopLeftFromPoint: NSMakePoint(20,20)];
	[window makeKeyAndOrderFront: self];
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
	// Insert code here to initialize your application
}

- (void)applicationWillTerminate:(NSNotification *)aNotification {
	// Insert code here to tear down your application
}

@end
