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
		NSScreen *mainScreen = [NSScreen mainScreen];
		const size_t retina = mainScreen.backingScaleFactor == 2.f ? 1 : 0;
		GLEssentialsGLView *view = [[GLEssentialsGLView alloc] init];
		NSWindow *window = [[NSWindow alloc] initWithContentRect: NSMakeRect(0, 0, param.image_w >> retina, param.image_h >> retina)
		                                               styleMask: NSWindowStyleMaskTitled
		                                                 backing: NSBackingStoreBuffered
		                                                   defer: NO
		                                                  screen: mainScreen];
		window.title = NSProcessInfo.processInfo.processName;
		[window setContentView: view];
		_controller = [[GLEssentialsWindowController alloc] initWithWindow: window];
	}
	return self;
}

- (void)applicationWillFinishLaunching:(NSNotification *)aNotification {
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
	NSWindow *window = [_controller window];
	[window makeKeyAndOrderFront: self];
}

- (void)applicationWillTerminate:(NSNotification *)aNotification {
}

@end
