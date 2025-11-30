/*
See the LICENSE_apple.txt file for this sample’s licensing information.

Abstract:
Implementation of a platform independent renderer class, which performs Metal setup and per frame rendering
*/

#import "MetalRenderer.h"
#import "param.h"

// Main class performing the rendering
@implementation MetalRenderer
{
	id<MTLDevice> _device;
	id<MTLCommandQueue> _commandQueue;
}

void back_to_caller(void *texture, void *bytes, size_t perRow)
{
	[(__bridge id<MTLTexture>)texture replaceRegion:MTLRegionMake2D(0, 0, param.image_w, param.image_h)
	                                    mipmapLevel:0
	                                      withBytes:bytes
	                                    bytesPerRow:perRow];
}

- (nonnull instancetype)initWithMTLDevice:(nonnull id<MTLDevice>)device
{
	self = [super init];
	if (self) {
		if (content_init()) {
			[[NSApplication sharedApplication] terminate:nil];
		}

		_device = device;
		_commandQueue = [_device newCommandQueue];
	}

	return self;
}

// Called whenever the view needs to render a frame.
- (void)drawInMTKView:(nonnull MTKView *)view
{
	@autoreleasepool {
		id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];

		// Get the drawable that will be presented at the end of the frame
		id<CAMetalDrawable> drawable = view.currentDrawable;

		// Get the drawable's texture to render content into
		id<MTLTexture> texture = drawable.texture;
		if (content_frame((__bridge void *)texture)) {
			[[NSApplication sharedApplication] terminate:nil];
		}

		// Request that the drawable texture be presented by the windowing system once drawing is done
		[commandBuffer presentDrawable:drawable];
		[commandBuffer commit];
	}
}

// Called whenever view changes orientation or is resized
- (void)mtkView:(nonnull MTKView *)view drawableSizeWillChange:(CGSize)size
{
}

- (void) dealloc
{
	content_deinit();
}

@end
