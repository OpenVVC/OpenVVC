#ifndef FG_COMPUTE_BLOCK_AVG_H
#define FG_COMPUTE_BLOCK_AVG_H

int16_t asm_fg_compute_block_avg_sse4(int16_t *dstSampleBlk8, uint32_t widthComp, uint16_t *pNumSamples,
                      uint8_t ySize, uint8_t xSize, uint8_t bitDepth);

#endif