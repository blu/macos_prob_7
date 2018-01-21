/*
 Copyright (C) 2015 Apple Inc. All Rights Reserved.
 See LICENSE.txt for this sample’s licensing information
 
 Abstract:
 Window controller subclass.
 */

#import "GLEssentialsWindowController.h"
#import "GLEssentialsFullscreenWindow.h"
#import "input.h"

extern unsigned input;

@interface GLEssentialsWindowController ()
{
    // Fullscreen window
    GLEssentialsFullscreenWindow *_fullscreenWindow;

    // Non-Fullscreen window (also the initial window)
    NSWindow* _standardWindow;
}
@end

@implementation GLEssentialsWindowController

- (instancetype)initWithWindow:(NSWindow *)window
{
    self = [super initWithWindow:window];

	if (self)
	{
		// Initialize to nil since it indicates app is not fullscreen
		_fullscreenWindow = nil;
    }

	return self;
}

- (void) goFullscreen
{
	// If app is already fullscreen...
	if(_fullscreenWindow)
	{
		//...don't do anything
		return;
	}

	// Allocate a new fullscreen window
	_fullscreenWindow = [[GLEssentialsFullscreenWindow alloc] init];

	// Resize the view to screensize
	NSRect viewRect = [_fullscreenWindow frame];

	// Set the view to the size of the fullscreen window
	[self.window.contentView setFrameSize: viewRect.size];

	// Set the view in the fullscreen window
	[_fullscreenWindow setContentView:self.window.contentView];

	_standardWindow = [self window];

	// Hide non-fullscreen window so it doesn't show up when switching out
	// of this app (i.e. with CMD-TAB)
	[_standardWindow orderOut:self];

	// Set controller to the fullscreen window so that all input will go to
	// this controller (self)
	[self setWindow:_fullscreenWindow];

	// Show the window and make it the key window for input
	[_fullscreenWindow makeKeyAndOrderFront:self];

}

- (void) goWindow
{
	// If controller doesn't have a full screen window...
	if(_fullscreenWindow == nil)
	{
		//...app is already windowed so don't do anything
		return;
	}

	// Get the rectangle of the original window
	NSRect viewRect = [_standardWindow frame];
	
	// Set the view rect to the new size
	[self.window.contentView setFrame:viewRect];

	// Set controller to the standard window so that all input will go to
	// this controller (self)
	[self setWindow:_standardWindow];

	// Set the content of the orginal window to the view
	[[self window] setContentView:_fullscreenWindow.contentView];

	// Show the window and make it the key window for input
	[[self window] makeKeyAndOrderFront:self];

	// Ensure we set fullscreen Window to nil so our checks for 
	// windowed vs. fullscreen mode elsewhere are correct
	_fullscreenWindow = nil;
}


- (void) keyDown:(NSEvent *)event
{
	unichar c = [[event charactersIgnoringModifiers] characterAtIndex:0];
	unsigned input_mask = 0;

	switch (c)
	{
		// Have f key toggle fullscreen
		case 'f':
			if(_fullscreenWindow == nil)
			{
				[self goFullscreen];
			}
			else
			{
				[self goWindow];
			}
			return;
		case 27:
			input_mask = INPUT_MASK_ESC;
			break;
		case ' ':
			input_mask = INPUT_MASK_ACTION;
			break;
		case '1':
			input_mask = INPUT_MASK_OPTION_1;
			break;
		case '2':
			input_mask = INPUT_MASK_OPTION_2;
			break;
		case '3':
			input_mask = INPUT_MASK_OPTION_3;
			break;
		case '4':
			input_mask = INPUT_MASK_OPTION_4;
			break;
		case 'a':
		case 'A':
			input_mask = INPUT_MASK_ALT_LEFT;
			break;
		case 'd':
		case 'D':
			input_mask = INPUT_MASK_ALT_RIGHT;
			break;
		case 'i':
		case 'I':
			input_mask = INPUT_MASK_UP;
			break;
		case 'j':
		case 'J':
			input_mask = INPUT_MASK_LEFT;
			break;
		case 'l':
		case 'L':
			input_mask = INPUT_MASK_RIGHT;
			break;
		case 'm':
		case 'M':
			input_mask = INPUT_MASK_DOWN;
			break;
		case 'w':
		case 'W':
			input_mask = INPUT_MASK_ALT_UP;
			break;
		case 'z':
		case 'Z':
			input_mask = INPUT_MASK_ALT_DOWN;
			break;
	}

	if (input_mask)
	{
		input |= input_mask;
		return;
	}

	// Allow other character to be handled (or not and beep)
	[super keyDown:event];
}

- (void) keyUp:(NSEvent *)event
{
	unichar c = [[event charactersIgnoringModifiers] characterAtIndex:0];
	unsigned input_mask = 0;

	switch (c)
	{
		// We have consumed keyDown and handled f key -- do nothing here
		case 'f':
			return;
		case 27:
			input_mask = INPUT_MASK_ESC;
			break;
		case ' ':
			input_mask = INPUT_MASK_ACTION;
			break;
		case '1':
			input_mask = INPUT_MASK_OPTION_1;
			break;
		case '2':
			input_mask = INPUT_MASK_OPTION_2;
			break;
		case '3':
			input_mask = INPUT_MASK_OPTION_3;
			break;
		case '4':
			input_mask = INPUT_MASK_OPTION_4;
			break;
		case 'a':
		case 'A':
			input_mask = INPUT_MASK_ALT_LEFT;
			break;
		case 'd':
		case 'D':
			input_mask = INPUT_MASK_ALT_RIGHT;
			break;
		case 'i':
		case 'I':
			input_mask = INPUT_MASK_UP;
			break;
		case 'j':
		case 'J':
			input_mask = INPUT_MASK_LEFT;
			break;
		case 'l':
		case 'L':
			input_mask = INPUT_MASK_RIGHT;
			break;
		case 'm':
		case 'M':
			input_mask = INPUT_MASK_DOWN;
			break;
		case 'w':
		case 'W':
			input_mask = INPUT_MASK_ALT_UP;
			break;
		case 'z':
		case 'Z':
			input_mask = INPUT_MASK_ALT_DOWN;
			break;
	}

	if (input_mask)
	{
		input &= ~input_mask;
		return;
	}

	// Allow other character to be handled (or not and beep)
	[super keyUp:event];
}

@end
