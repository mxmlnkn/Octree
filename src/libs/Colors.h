#include <stdint.h>

namespace Colors {

uint32_t Red   = 0xff0000;
uint32_t Green = 0x00ff00;
uint32_t Blue  = 0x0000ff;

uint32_t BuPu[9] = { 0xf7fcfd, 0xe0ecf4, 0xbfd3e6,
                     0x9ebcda, 0x8c96c6, 0x8c6bb1,
                     0x88419d, 0x810f7c, 0x4d004b };

uint32_t getRed( uint32_t color ) {
    return ( color & 0xff0000 ) >> 16;
}
uint32_t getGreen( uint32_t color ) {
    return ( color & 0x00ff00 ) >> 16;
}
uint32_t getBlue( uint32_t color ) {
    return ( color & 0x0000ff ) >> 16;
}

}
