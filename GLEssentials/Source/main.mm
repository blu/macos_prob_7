/*
 Copyright (C) 2015 Apple Inc. All Rights Reserved.
 See LICENSE.txt for this sample’s licensing information
 
 Abstract:
 Standard AppKit entry point.
 */

#import <Cocoa/Cocoa.h>
#import "AppDelegate.h"
#import "param.h"
#import "stream.hpp"

int main(int argc, char * argv[])
{
	// set up cin, cout and cerr substitute streams
	stream::cin.open(stdin);
	stream::cout.open(stdout);
	stream::cerr.open(stderr);

	param.platform_idx = 0;
	param.device_idx = -1U;
	param.flags = 0;
	param.kern_param_type = PARAM_TYPE_BUFFER;
	param.fsaa = 0;
	param.image_w = 512;
	param.image_h = 512;
	param.bitness[0] = 0;
	param.bitness[1] = 0;
	param.bitness[2] = 0;
	param.bitness[3] = 0;
	param.frames = -1U;

	// read render setup from CLI
	const int result_cli = parseCLI(argc, argv, &param);

	if (0 != result_cli)
		return result_cli;

	@autoreleasepool {
		NSApplication *application = [NSApplication sharedApplication];
		[application activateIgnoringOtherApps: YES];

		// provide app menu with one item: quit
		NSMenuItem *item = [[NSMenuItem alloc] init];
		[item setSubmenu: [[NSMenu alloc] init]];
		[item.submenu addItem: [[NSMenuItem alloc] initWithTitle: [@"Quit " stringByAppendingString: NSProcessInfo.processInfo.processName] action:@selector(terminate:) keyEquivalent:@"q"]];
		[application setMainMenu: [[NSMenu alloc] init]];
		[application.mainMenu addItem: item];

		AppDelegate *appDelegate = [[AppDelegate alloc] init];
		[application setDelegate: appDelegate];
		[application run];
	}
	return EXIT_SUCCESS;
}
