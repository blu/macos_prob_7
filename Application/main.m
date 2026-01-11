/*
See the LICENSE_apple.txt file for this sample’s licensing information.

Abstract:
Application entry point for macOS.
*/

#import <Cocoa/Cocoa.h>
#import "AppDelegate.h"
#import "param.h"

int main(int argc, const char * argv[])
{
	param.platform_idx = 0;
	param.device_idx = -1U;
	param.flags = 0;
	param.kern_param_type = PARAM_TYPE_BUFFER;
	param.image_w = 512;
	param.image_h = 512;
	param.image_hz = 60;
	param.frames = -1U;
	param.workgroup_size = 0;

	// read render setup from CLI
	const int result_cli = parseCLI(argc, argv, &param);

	if (0 != result_cli)
		return result_cli;

	@autoreleasepool {
		NSApplication *application = [NSApplication sharedApplication];
		[application setActivationPolicy:NSApplicationActivationPolicyRegular];

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
