#include <stdint.h>
#include <math/TVector.h>

namespace Colors {

uint32_t Red   = 0xff0000;
uint32_t Green = 0x00ff00;
uint32_t Blue  = 0x0000ff;

uint32_t BuPu[9] = { 0xf7fcfd, 0xe0ecf4, 0xbfd3e6,
                     0x9ebcda, 0x8c96c6, 0x8c6bb1,
                     0x88419d, 0x810f7c, 0x4d004b };
const size_t BuPuLength = sizeof(BuPu) / sizeof(uint32_t);
uint32_t Own1[9] = { 0xefdd10 /* dark yellow */, 0xa75e07 /* brown */, 
                     0x6a9aff /* light blue  */, 0xa75eff /* lila  */,
                     0x63a338 /* dark green  */ };
const size_t Own1Length = sizeof(Own1) / sizeof(uint32_t);
uint32_t SatSpec[6] = { 0xff0000, 0xffff00, 0x00ff00,
                        0x00ffff, 0x0000ff, 0xff00ff };
const size_t SatSpecLength = sizeof(SatSpec) / sizeof(uint32_t);
                     
uint32_t getRed( uint32_t color ) {
    return ( color & 0xff0000 ) >> 16;
}
uint32_t getGreen( uint32_t color ) {
    return ( color & 0x00ff00 ) >> 8;
}
uint32_t getBlue( uint32_t color ) {
    return ( color & 0x0000ff ) >> 0;
}

Vec<int,3> getColorVector( const uint32_t color ) {
    return Vec<int,3>( getRed(color), getGreen(color), getBlue(color) );
}

} /* namespace Colors */
