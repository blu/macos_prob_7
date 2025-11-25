/*
 Copyright (C) 2015 Apple Inc. All Rights Reserved.
 See LICENSE.txt for this sample’s licensing information
 
 Abstract:
 The Application Delegate.
 */

#import "AppDelegate.h"
#import "OpenGLView.h"
#import "param.h"

@interface AppDelegate ()
{
	NSWindowController *_controller;
}
@end

@implementation AppDelegate

- (instancetype) init
{
	self = [super init];
	if (self) {
		NSScreen *mainScreen = [NSScreen mainScreen];
		const size_t retina = mainScreen.backingScaleFactor == 2.f ? 1 : 0;
		const NSRect rect = NSMakeRect(0, 0, param.image_w >> retina, param.image_h >> retina);
		OpenGLView *view = [[OpenGLView alloc] init];
		NSWindow *window = [[NSWindow alloc] initWithContentRect:rect
		                                               styleMask:NSWindowStyleMaskTitled
		                                                 backing:NSBackingStoreBuffered
		                                                   defer:NO
		                                                  screen:mainScreen];

		window.title = NSProcessInfo.processInfo.processName;
		window.contentView = view;

		_controller = [[NSWindowController alloc] initWithWindow:window];
	}

	return self;
}

- (void)applicationWillFinishLaunching:(NSNotification *)aNotification
{
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
	[_controller.window makeKeyAndOrderFront:self];
}

- (void)applicationWillTerminate:(NSNotification *)aNotification
{
	[_controller.window.contentView deinit];
}

@end
