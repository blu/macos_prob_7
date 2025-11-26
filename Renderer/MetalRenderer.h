/*
See the LICENSE_apple.txt file for this sample’s licensing information.

Abstract:
Header for a platform independent renderer class, which performs Metal setup and per frame rendering.
*/

#import <MetalKit/MetalKit.h>
#import <Foundation/Foundation.h>

@interface MetalRenderer : NSObject<MTKViewDelegate>

- (nonnull instancetype)initWithMTLDevice:(nonnull id<MTLDevice>)device;
- (void) dealloc;

@end
